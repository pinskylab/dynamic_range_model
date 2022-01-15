#############
# load packages and data
#############

set.seed(42)
library(tidyverse)
library(tidybayes)
library(Cairo)
library(here)
library(magrittr)
library(rstan)
library(Matrix)
library(ggridges)
library(rstanarm)
library(geosphere)
library(ggridges)
library(purrr)

funs <- list.files("functions")
sapply(funs, function(x) source(file.path("functions",x)))

rstan_options(javascript=FALSE, auto_write =TRUE)

dat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv"))
# dat %<>% filter(length >17)
dat_test <- read_csv(here("processed-data","flounder_catch_at_length_fall_testing.csv"))


#############
# make model decisions
#############
trim_to_abundant_patches=FALSE
if(trim_to_abundant_patches==TRUE){
  focal_patch_n = 7 # adjust accordingly
}
do_dirichlet = 1
eval_l_comps = 0 # evaluate length composition data? 0=no, 1=yes
T_dep_mortality = 0 # CURRENTLY NOT REALLY WORKING
T_dep_recruitment = 1 # think carefully before making more than one of the temperature dependencies true
spawner_recruit_relationship = 0
run_forecast=1
time_varying_f = TRUE

if(time_varying_f==TRUE){
# the f-at-age data starts in 1982; fill in the previous years with the earliest year of data
dat_f_age_prep <- read_csv(here("processed-data","summer_flounder_F_by_age.csv")) %>%
  rename_with(str_to_lower)
f_early <- expand_grid(year=seq(1972, 1981, 1), age=unique(dat_f_age_prep$age)) %>% 
  left_join(dat_f_age_prep %>% filter(year==1982) %>% select(age, f)) 
dat_f_age_prep <- bind_rows(dat_f_age_prep, f_early)
}

make_data_plots <- FALSE

if(make_data_plots==TRUE){
  # how much variation is there in length frequency over time?
  gglength <- dat %>% 
    mutate(lat_floor = floor(lat)) %>% 
    filter(between(lat_floor, 37,41)) %>% 
    group_by(length, year, lat_floor) %>% 
    summarise(total = sum(number_at_length)) %>% 
    mutate(year = as.character(year)) %>% 
    ggplot(aes(x=length, y=year)) +
    geom_density_ridges(aes(height=..density..,
                            weight=total),    
                        scale= 4,
                        stat="density") +
    scale_y_discrete(expand = c(0, 0))  + 
    facet_wrap(~lat_floor) + 
    coord_flip()
  # gglength  
  ggsave(gglength, filename=here("results","summer_flounder_length_freq.png"))
  
  # by patch -- 
  ggpatchlength <- dat %>% 
    mutate(lat_floor = floor(lat)) %>% 
    group_by(length, lat_floor, year) %>% 
    summarise(total = sum(number_at_length)) %>% 
    mutate(year = as.character(year)) %>% 
    ggplot(aes(x=length, y=year)) +
    geom_density_ridges(aes(height=..density..,
                            weight=total),    
                        scale= 4,
                        stat="density") +
    scale_y_discrete(expand = c(0, 0)) +
    facet_wrap(~lat_floor, scales = "free_y")
  # ggpatchlength  
  ggsave(ggpatchlength, filename=here("results","summer_flounder_length_freq_by_patch.png"), height=8, width=5, dpi=160)
  # not sure this doesn't have any data in the other patches, maybe they are missing years? 
}

##########
# prep data for fitting
##########

if(trim_to_abundant_patches==TRUE){
  # reshape fish data 
  use_patches <- dat %>% 
    mutate(lat_floor = floor(lat)) %>% 
    group_by(lat_floor) %>% 
    summarise(total = sum(number_at_length)) %>% 
    arrange(-total) %>% 
    slice(1:7) # CHECK THIS MANUALLY TO BE SURE IT'S SANE, AND PATCHES ARE CONTIGUOUS
  
  patches <- sort(unique(use_patches$lat_floor))
  np = length(patches) 
}

if(trim_to_abundant_patches==FALSE){
  # reshape fish data 
  use_patches <- dat %>% 
    mutate(lat_floor = floor(lat)) %>% 
    group_by(lat_floor) %>% 
    summarise(total = sum(number_at_length))
  
  patches <- sort(unique(use_patches$lat_floor))
  np = length(patches) 
}

dat_train_lengths <- dat %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(length, year, lat_floor) %>% 
  summarise(sum_num_at_length = sum(number_at_length)) %>% 
  filter(lat_floor %in% patches)%>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_test_lengths <- dat_test %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(length, year, lat_floor) %>% 
  summarise(sum_num_at_length = sum(number_at_length)) %>% 
  filter(lat_floor %in% patches)%>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_train_dens <- dat %>% 
  mutate(lat_floor = floor(lat)) %>% 
  filter(lat_floor %in% patches) %>% 
  group_by(haulid) %>% 
  mutate(dens = sum(number_at_length)) %>% # get total no. fish in each haul, of any size
  group_by(year, lat_floor) %>% 
  summarise(mean_dens = mean(dens)) %>%  # get mean density (all sizes) / haul for the patch*year combo 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_test_dens <- dat_test %>% 
  mutate(lat_floor = floor(lat)) %>% 
  filter(lat_floor %in% patches) %>% 
  group_by(haulid) %>% 
  mutate(dens = sum(number_at_length)) %>% # get total no. fish in each haul, of any size
  group_by(year, lat_floor) %>% 
  summarise(mean_dens = mean(dens)) %>%  # get mean density (all sizes) / haul for the patch*year combo 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

# get time dimension
years <- sort(unique(dat_train_lengths$year)) 
years_proj <- sort(unique(dat_test_lengths$year))
ny <- length(years)
ny_proj <- length(years_proj)

# get patch area 
patchdat <- dat %>% 
  select(lat, lon) %>%
  distinct() %>%
  mutate(lat_floor = floor(lat)) %>% 
  filter(lat_floor %in% use_patches$lat_floor) %>% 
  group_by(lat_floor) %>% 
  summarise(max_lon=max(lon),
            min_lon=min(lon)) %>% 
  rowwise() %>% 
  mutate(lon_dist = distGeo(p1=c(max_lon, lat_floor+0.5), p2=c(min_lon, lat_floor+0.5))/1000, # get distance between the furthest longitudes in km, at the midpoint of the lat band 
         patch_area_km2 = lon_dist * 111) %>%  # 1 degree latitude = 111 km 
  select(lat_floor, patch_area_km2) %>% 
  filter(lat_floor %in% patches) %>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))
meanpatcharea <- mean(patchdat$patch_area_km2)

dat_train_sbt <- dat %>%   
  mutate(lat_floor = floor(lat)) %>% 
  group_by(lat_floor, year) %>% 
  summarise(sbt = mean(btemp, na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(lat_floor %in% patches)%>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_test_sbt <- dat_test %>%   
  mutate(lat_floor = floor(lat)) %>% 
  group_by(lat_floor, year) %>% 
  summarise(sbt = mean(btemp, na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(lat_floor %in% patches)%>% 
  mutate(patch = as.integer(as.factor(lat_floor))) %>% 
  filter(!is.na(sbt))

# how many rows should each df have?
nrow(dat_train_sbt)==(np*ny)
nrow(dat_test_sbt) == (np*ny_proj) 

# some SBT data are missing: lat 35, 36, and 37 in 2008
# THIS IS VERY HARD CODED! DANGER
library(lme4)
sbt_lme <- lmer(btemp ~ lat + (1|year), data=bind_rows(dat, dat_test))
sbt_test_fill <- data.frame(
  lat_floor = c(35, 36, 37),
  year = rep(2008, 3),
  patch = c(1,2,3),
  sbt = c(predict(sbt_lme, newdata=data.frame(lat=35.5, year=2008)),
          predict(sbt_lme, newdata=data.frame(lat=36.5, year=2008)),
          predict(sbt_lme, newdata=data.frame(lat=35.5, year=2008))
  ))
dat_test_sbt <- bind_rows(dat_test_sbt, sbt_test_fill)

# set fixed parameters from stock assessment
loo = 83.6
k = 0.14
m = 0.25
age_at_maturity = 2
t0=-.2
cv= 0.2 # guess
min_age = 0
max_age = 20

if(time_varying_f==TRUE){
# now that we have the max_age, fill in f for years above 7 (since the f for age=7 is really for 7+)
older_ages <- expand_grid(age=seq(max(dat_f_age_prep$age)+1, max_age, 1), year= unique(dat_f_age_prep$year)) %>% 
  left_join(dat_f_age_prep %>% filter(age==max(age)) %>% select(year, f))
dat_f_age_prep %<>% bind_rows(older_ages)
}

# make length to age conversions
length_at_age_key <-
  generate_length_at_age_key(
    min_age = min_age,
    max_age = max_age,
    cv = cv,
    linf = loo,
    k = k,
    t0 = t0,
    time_step = 1,
    linf_buffer = 1.5
  )

length_at_age_key %>% 
  filter(age > 5) %>% 
  ggplot(aes(age, length_bin, fill = p_bin)) + 
  geom_hline(aes(yintercept = loo)) +
  geom_tile() + 
  scale_fill_viridis_c()

l_at_a_mat <- length_at_age_key %>% 
  select(age, length_bin, p_bin) %>% 
  pivot_wider(names_from = length_bin, values_from = p_bin) %>% 
  ungroup() %>% 
  select(-age) %>% 
  as.matrix()

# prep f data
if(time_varying_f==TRUE){
  dat_f_age <- dat_f_age_prep %>% 
    filter(year %in% years) 
  
  dat_f_age_proj <- dat_f_age_prep %>% 
    filter(year %in% years_proj) %>%
    bind_rows(dat_f_age %>% filter(year==max(year))) # need final year of training data to initialize projection
} else {
  f_prep=0.334
}

lbins <- unique(length_at_age_key$length_bin)
# lbins <- sort(unique(dat_train_lengths$length))
n_lbins <- length(lbins) 

n_ages <- nrow(l_at_a_mat)

# now that we have n_ages, calculate weight at age
# THIS IS VERY HARD CODED, AND ALSO NOT PROBABILISTIC
# NEED TO FIX
wt_at_age_prep <- read_csv(here("processed-data","summer_flounder_wt_at_age.csv")) %>% 
  filter(!Age %in% seq(7, 10, 1)) %>% 
  mutate(Age = gsub("over7",7,Age),
         Age = as.numeric(Age),
         Age = Age + 1) %>% # start at 1 not 0
  group_by(Age) %>% 
  summarise(wt = mean(Wt)) %>% # average over all years 
  ungroup() %>% 
  arrange(Age)
# fill in plus ages
wt_at_age_add <- data.frame(Age = seq(max(wt_at_age_prep+1),n_ages,1), wt = slice_tail(wt_at_age_prep)$wt)
wt_at_age <- rbind(wt_at_age_prep, wt_at_age_add)$wt

# now that years are defined above, convert them into indices in the datasets
# be sure all these dataframes have exactly the same year range! 

dat_train_dens$year = as.integer(as.factor(dat_train_dens$year))
dat_test_dens$year = as.integer(as.factor(dat_test_dens$year))
dat_train_lengths$year = as.integer(as.factor(dat_train_lengths$year))
#dat_test_lengths$year = as.integer(as.factor(dat_test_lengths$year))
dat_test_sbt$year= as.integer(as.factor(dat_test_sbt$year))
dat_train_sbt$year= as.integer(as.factor(dat_train_sbt$year))
if(time_varying_f==TRUE){
  dat_f_age$year = as.integer(as.factor(dat_f_age$year))
  dat_f_age_proj$year = as.integer(as.factor(dat_f_age_proj$year))
}

# make matrices/arrays from dfs
len <- array(0, dim = c(np, n_lbins, ny)) 
for(p in 1:np){
  for(l in 1:n_lbins){
    for(y in 1:ny){
      tmp <- dat_train_lengths %>% filter(patch==p, round(length)==lbins[l], year==y) 
      if (nrow(tmp) > 0){
        len[p,l,y] <- tmp$sum_num_at_length
      }
    }
  }
}

plot(len[4,,20])

dens <- array(NA, dim=c(np, ny))
for(p in 1:np){
  for(y in 1:ny){
    tmp2 <- dat_train_dens %>% filter(patch==p, year==y) 
    #   left_join(patchdat, by = c("lat_floor","patch"))%>% 
    #   mutate(mean_dens = mean_dens * patch_area_km2)
    # dens[p,y] <- tmp2$mean_dens
    dens[p,y] <- tmp2$mean_dens *meanpatcharea
  }
}

sbt <- array(NA, dim=c(np,ny))
for(p in 1:np){
  for(y in 1:ny){
    tmp3 <- dat_train_sbt %>% filter(patch==p, year==y) 
    sbt[p,y] <- tmp3$sbt
  }
}

sbt_proj <- array(NA, dim=c(np,ny_proj))
for(p in 1:np){
  for(y in 1:ny_proj){
    tmp6 <- dat_test_sbt %>% filter(patch==p, year==y) 
    sbt_proj[p,y] <- tmp6$sbt
  }
}


f <- array(NA, dim=c(n_ages,ny))
for(a in min_age:max_age){
  for(y in 1:ny){
    if(time_varying_f==TRUE){
      tmp4 <- dat_f_age %>% filter(age==a, year==y) 
      f[a+1,y] <- tmp4$f # add 1 because matrix indexing starts at 1 not 0
    } else{
      f[a+1,y] <- f_prep
    }
  }
}
f_proj <- array(NA, dim=c(n_ages,(ny_proj+1)))
for(a in min_age:max_age){
  for(y in 1:(ny_proj+1)){
    if(time_varying_f==TRUE){
      tmp5 <- dat_f_age_proj %>% filter(age==a, year==y) 
      f_proj[a+1,y] <- tmp5$f # add 1 because matrix indexing starts at 1 not 0
    } else{
      f_proj[a+1,y] <-f_prep
    }
  }
}

a <- seq(min_age, max_age)

check <- a %*% l_at_a_mat
######
# fit model
######

stan_data <- list(
  np=np,
  n_ages=n_ages,
  ny_train=ny,
  ny_proj=ny_proj,
  n_lbins=n_lbins,
  n_p_l_y = len,
  abund_p_y = dens,
  sbt = sbt,
  sbt_proj=sbt_proj,
  m=m,
  f=f,
  f_proj=f_proj,
  k=k,
  loo=loo,
  t0=t0,
  cv=cv,
  length_50_sel_guess=20, # THIS IS A RANDOM GUESS, I can't find help in the stock assessment
  n_lbins = n_lbins, 
  age_sel = 0,
  bin_mids=lbins+0.5, # also not sure if this is the right way to calculate the midpoints
  sel_100 = 3, # not sure if this should be 2 or 3. it's age 2, but it's the third age category because we start at 0, which I think Stan will classify as 3...?
  age_at_maturity = age_at_maturity,
  l_at_a_key = l_at_a_mat,
  wt_at_age = wt_at_age,
  do_dirichlet = do_dirichlet,
  eval_l_comps = eval_l_comps, # evaluate length composition data? 0=no, 1=yes
  T_dep_mortality = T_dep_mortality, # CURRENTLY NOT REALLY WORKING
  T_dep_recruitment = T_dep_recruitment, # think carefully before making more than one of the temperature dependencies true
  spawner_recruit_relationship = spawner_recruit_relationship, 
  run_forecast=run_forecast
)

warmups <- 1000
total_iterations <- 2000
max_treedepth <-  10
n_chains <-  1
n_cores <- 1

stan_model_fit <- stan(file = here::here("src","process_sdm_stock_recruit.stan"), # check that it's the right model!
                       data = stan_data,
                       chains = n_chains,
                       warmup = warmups,
                       #     init = list(list(log_mean_recruits = log(1000),
                       #                      theta_d = 1,
                       #                     ssb0=1000000)),
                       iter = total_iterations,
                       cores = n_cores,
                       refresh = 250,
                       control = list(max_treedepth = max_treedepth,
                                      adapt_delta = 0.85)
)

# a = rstan::extract(stan_model_fit, "theta_d")
# write_rds(stan_model_fit,"sigh.rds")
# hist(a$sigma_obs)
# rstanarm::launch_shinystan(stan_model_fit)


# plot important parameters 
plot(stan_model_fit, pars=c('sigma_r','sigma_obs','d','width','Topt','alpha','beta_obs','theta_d'))
plot(stan_model_fit, pars=c('sigma_r','sigma_obs','d','alpha','beta_obs','theta_d'))


# assess abundance fits

hist(extract(stan_model_fit, "mean_recruits")$mean_recruits)

abund_p_y <- dat_train_dens %>%
  mutate(abundance = mean_dens * meanpatcharea) 
# left_join(patchdat, by = c("lat_floor", "patch")) %>% 
# group_by(patch, year) %>% 
# summarise(abundance = sum(mean_dens *meanpatcharea)) %>% 
# ungroup()


abund_p_y_hat <- tidybayes::spread_draws(stan_model_fit, dens_p_y_hat[patch,year])

abundance_v_time <- abund_p_y_hat %>% 
  ggplot(aes(year, dens_p_y_hat)) + 
  stat_lineribbon() + 
  geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch, scales = "free_y") +
  labs(x="Year",y="Abundance") + 
  scale_fill_brewer()
abundance_v_time
ggsave(abundance_v_time, filename=here("results","density_v_time_no_length_comps.png"), width=7, height=4)

# assess length comp fits

n_p_l_y_hat <- tidybayes::gather_draws(stan_model_fit, n_p_l_y_hat[year,patch,length], n = 500)

# neff <- tidybayes::gather_draws(stan_model_fit, n_eff[patch,year], n = 500)

#neff <- tidybayes::gather_draws(stan_model_fit, n_eff[patch,year], n = 500)

p = 5

dat_train_lengths <- dat_train_lengths %>% 
  group_by(patch, year) %>% 
  mutate(p_length = sum_num_at_length / sum(sum_num_at_length))

dat_train_lengths %>% 
  group_by(year,patch) %>% 
  summarise(n = sum(sum_num_at_length)) %>% 
  ggplot(aes(year, n, color =factor(patch))) + 
  geom_point()

n_p_l_y_hat %>% 
  ungroup() %>% 
  filter(patch == p) %>% 
  group_by(patch, year, .iteration) %>% 
  mutate(pvalue = .value / sum(.value)) %>% 
  ggplot(aes(length, pvalue)) + 
  stat_lineribbon() + 
  geom_point(data = dat_train_lengths %>% filter(patch == p), aes(length,p_length), color = "red", alpha = 0.2) +
  facet_wrap(~year, scales = "free_y")

# length frequency over time
l_freq_time <- n_p_l_y_hat %>% 
  ungroup() %>% 
  ggplot(aes(x=length, y=..density.., weight=.value)) + 
  geom_histogram(bins=50) +
  facet_grid(patch~year)

ggsave(l_freq_time, filename=here("results","length_freq_time_flounder.png"), scale=1.5, width=15, height=4)

# is there a temperature - recruitment relationship? 

## first need to figure out which lengths correspond to age 0 
length_at_age_key %>% 
  ggplot(aes(x=length_bin, y=p_bin)) +
  geom_line() + 
  facet_wrap(~age)

## let's call recruits <5cm

selectivity_at_bin <- gather_draws(stan_model_fit, selectivity_at_bin[n_lbins], n=500)

mean_sel <- selectivity_at_bin %>% 
  group_by(n_lbins) %>% 
  summarise(mean_sel = mean(.value)) %>% 
  rename(length = n_lbins)

# get proportion of recruits adjusted by selectivity
recruits <- n_p_l_y_hat %>% 
  left_join(mean_sel, by="length") %>% 
  mutate(recruit = ifelse(length <=5, "yes", "no"),
         val_post_sel = .value / mean_sel) %>% 
  group_by(recruit, patch, year, .draw) %>% 
  summarise(sumcount = sum(val_post_sel)) %>% 
  ungroup() %>% 
  group_by(patch, year, .draw) %>% 
  mutate(prop_recruit = sumcount / sum(sumcount)) %>% 
  filter(recruit == "yes")

recruits %>% group_by(patch, year) %>% 
  mutate(mean_prop_rec = mean(prop_recruit)) %>% 
  left_join(dat_train_sbt, by=c('patch','year')) %>% 
  ggplot(aes(x=sbt, y=mean_prop_rec, color=year, fill=year)) +
  geom_point()

# plot recruitment deviates
# note that the length of rec_dev is actually 34 not 35
gather_draws(stan_model_fit, rec_dev[ny]) %>% 
  ggplot(aes(x=ny, y=.value)) + 
  stat_lineribbon() + 
  scale_fill_brewer()

# plot raw
raw_v_time <- gather_draws(stan_model_fit, raw[ny]) %>%
  ggplot(aes(x=ny, y=.value)) + 
  stat_lineribbon() + 
  scale_fill_brewer()
ggsave(raw_v_time, filename=here("results","raw_v_time.png"))

# plot actual recruitment
n_p_a_y_hat <- gather_draws(stan_model_fit, n_p_a_y_hat[np, n_ages, ny])

n_p_a_y_hat %>%
  filter(n_ages==1) %>%
  ggplot(aes(x=ny, y=.value)) +
  stat_lineribbon() +
  scale_fill_brewer() +
  facet_wrap(~np)

# plot all ages over time
n_p_a_y_hat %>% 
  ggplot(aes(x=ny, y=.value)) +
  stat_lineribbon() +
  scale_fill_brewer() +
  facet_grid(np~n_ages)

proj_n_p_a_y_hat %>% 
  ggplot(aes(x=`(ny_proj + 1)`, y=.value)) +
  stat_lineribbon() +
  scale_fill_brewer() +
  facet_grid(np~n_ages)

# detection stats 
spread_draws(stan_model_fit, theta[patch,year]) %>% 
  ggplot(aes(x=year, y=theta)) + 
  stat_lineribbon() + 
  facet_wrap(~patch) +
  scale_fill_brewer()

beta_obs <- extract(stan_model_fit, "beta_obs")$beta_obs
hist(beta_obs)

# save model run if desired
saveRDS(stan_model_fit, here("results","summer_flounder_stan_fit.rds"))

###########
# calculate summary statistics and evaluate range shifts
###########

# centroid position by year 
dat_centroid <- abund_p_y %>% 
  group_by(year) %>% 
  summarise(centroid_lat = weighted.mean(x=patch, w=abundance))

# model fit centroid -- should eventually estimate in model for proper SE -- just exploring here
est_centroid <- abund_p_y_hat %>% 
  group_by(year, .iteration) %>%  # IS THIS SUPPOSED TO BE .ITERATION? CHECK WHEN MODEL IS RUN FOR LONGER 
  summarise(centroid_lat = weighted.mean(x=patch, w=dens_p_y_hat)) %>% 
  ungroup()

gg_centroid <- est_centroid %>% 
  ggplot(aes(year, centroid_lat)) + 
  stat_lineribbon() + 
  scale_fill_brewer() +
  geom_point(data = dat_centroid, aes(year, centroid_lat), color = "red") +
  theme(legend.position = "none") +
  labs(x="Year", y="Patch", title="Centroid Position", fill="Credible Interval")
gg_centroid
ggsave(gg_centroid, filename=here("results","centroid_v_time_no_length_comps.png"))

# centroid didn't shift at all!

# patch abundance fraction every year 
est_patch_abund <- abund_p_y_hat %>% 
  group_by(year, patch) %>% 
  summarise(abundance = mean(dens_p_y_hat))

observed_abundance_tile <- abund_p_y %>% 
  ggplot(aes(x=year, y=patch, fill=abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0, 36, 4)) +
  scale_y_continuous(breaks=seq(1, np, 1)) +
  labs(title="Observed", x="Year", y="Patch", fill="Abundance")


estimated_abundance_tile <- est_patch_abund %>% 
  ggplot(aes(x=year, y=patch, fill=abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0, 36, 4)) +
  scale_y_continuous(breaks=seq(1, np, 1)) +
  labs(title="Estimated")

ggsave(observed_abundance_tile, filename=here("results","observed_abundance_v_time_tileplot_no_length_comps.png"))
ggsave(estimated_abundance_tile, filename=here("results","estimated_abundance_v_time_tileplot_no_length_comps.png"))

# who's doing the colonizing?
dat_train_lengths %>% 
  group_by(patch, length) %>% 
  arrange(year) %>% 
  mutate(logratio = log(sum_num_at_length / lag(sum_num_at_length))) %>% 
  filter(logratio < Inf, logratio > -Inf) %>% 
  ungroup() %>% 
  ggplot(aes(x=year, y=logratio, group=length, color=length)) +
  geom_point() + 
  geom_line() +
  scale_color_viridis_c() +
  facet_wrap(~patch)

dat_train_lengths %>% 
  ggplot(aes(x=year, y=sum_num_at_length, fill=length)) + 
  geom_bar(stat="identity") +
  facet_wrap(~patch) +
  scale_fill_viridis_c() +
  theme_bw()


dat_train_lengths %>% 
  ggplot(aes(x=year, y=sum_num_at_length, fill=length)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_viridis_c() +
  theme_bw()  

# plot temperature difference from optimum over space and time
Topt <- extract(stan_model_fit, "Topt")$Topt
hist(Topt) # note, not normally distributed

dat_train_sbt %>%
  mutate(Tdiff = sbt - median(Topt)) %>%
  ggplot(aes(x=year, y=patch, fill=Tdiff)) +
  geom_tile() + 
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)

########
# evaluate forecast
########
proj_dens_p_y_hat <- gather_draws(stan_model_fit, proj_dens_p_y_hat[np, ny_proj])

proj_abund_p_y <- dat_test_dens %>%
  mutate(abundance = mean_dens * meanpatcharea)
# left_join(patchdat, by = c("lat_floor", "patch")) %>% 
# group_by(patch, year) %>% 
# summarise(abundance = sum(mean_dens *patch_area_km2)) %>% 
#  ungroup()

proj_observed_abundance_tile <- proj_abund_p_y %>% 
  mutate(Year = (year + min(years_proj) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>% 
  ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="Observed")


proj_est_patch_abund <- proj_dens_p_y_hat %>% 
  group_by(ny_proj, np) %>% 
  summarise(abundance = mean(.value))

proj_estimated_abundance_tile <- proj_est_patch_abund %>% 
  mutate(Year = (ny_proj + min(years_proj) - 1), Latitude = (np + min(patches) - 1), Abundance=abundance) %>% 
  ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="Estimated")
ggsave(proj_estimated_abundance_tile, filename=here("results","proj_estimated_abundance_v_time_tileplot.png"), scale=0.9)
ggsave(proj_observed_abundance_tile, filename=here("results","proj_observed_abundance_v_time_tileplot.png"), scale=0.9)

proj_abundance_v_time <- proj_dens_p_y_hat %>% 
  rename(patch=np) %>% 
  ggplot(aes(ny_proj, .value)) + 
  stat_lineribbon() + 
  geom_point(data = proj_abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch, scales = "free_y") +
  labs(x="Year",y="Abundance") + 
  scale_x_continuous(breaks=seq(0, 10, 2), limits=c(0, 10)) +
  theme(legend.position="none") +
  scale_fill_brewer()
ggsave(proj_abundance_v_time, filename=here("results","proj_density_v_time.png"), width=7, height=4)

# centroid position by year 
dat_centroid_proj <- proj_abund_p_y %>% 
  group_by(year) %>% 
  summarise(centroid_lat = weighted.mean(x=patch, w=abundance)) %>%
  mutate(Year = (year + min(years_proj) - 1), Latitude = (centroid_lat + min(patches) - 1))

# model fit centroid -- should eventually estimate in model for proper SE -- just exploring here
est_centroid_proj <- proj_dens_p_y_hat %>% 
  group_by(ny_proj, .iteration) %>%  # IS THIS SUPPOSED TO BE .ITERATION? CHECK WHEN MODEL IS RUN FOR LONGER 
  summarise(centroid_lat = weighted.mean(x=np, w=.value)) %>% 
  ungroup()

gg_centroid_proj_prep <- dat_centroid_proj %>% 
  ggplot(aes(Year, Latitude)) + 
  geom_point(color = "red") +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(37.6, 39.2, 0.2), limits=c(37.6, 39.2)) +
  theme(legend.position = "none") +
  labs(title="Centroid Position") 

gg_centroid_proj <- est_centroid_proj %>% 
  mutate(Year = (ny_proj + min(years_proj) - 1), Latitude = (centroid_lat + min(patches) - 1)) %>% 
  ggplot(aes(Year, Latitude)) + 
  stat_lineribbon() + 
  scale_fill_brewer() +
  geom_point(data = dat_centroid_proj, aes(Year, Latitude), color = "red") +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(37.6, 39.2, 0.2)) +
  theme(legend.position = "none") +
  labs(title="Centroid Position") 
ggsave(gg_centroid_proj_prep, filename=here("results","proj_centroid_v_time_prep.png"), scale=0.9)
ggsave(gg_centroid_proj, filename=here("results","proj_centroid_v_time.png"), scale=0.9)
