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

# the f-at-age data starts in 1982; fill in the previous years with the earliest year of data
dat_f_age_prep <- read_csv(here("processed-data","summer_flounder_F_by_age.csv")) %>%
  rename_with(str_to_lower)
f_early <- expand_grid(year=seq(1972, 1981, 1), age=unique(dat_f_age_prep$age)) %>% 
  left_join(dat_f_age_prep %>% filter(year==1982) %>% select(age, f)) 
dat_f_age_prep <- bind_rows(dat_f_age_prep, f_early)

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

# reshape fish data 
top_patches <- dat %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(lat_floor) %>% 
  summarise(total = sum(number_at_length)) %>% 
  arrange(-total) %>% 
  slice(1:7) # CHECK THIS MANUALLY TO BE SURE IT'S SANE, AND PATCHES ARE CONTIGUOUS

dat_train_lengths <- dat %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(length, year, lat_floor) %>% 
  summarise(sum_num_at_length = sum(number_at_length)) %>% 
  filter(lat_floor %in% top_patches$lat_floor)%>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_test_lengths <- dat_test %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(length, year, lat_floor) %>% 
  summarise(sum_num_at_length = sum(number_at_length)) %>% 
  filter(lat_floor %in% top_patches$lat_floor)%>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_train_dens <- dat %>% 
  mutate(lat_floor = floor(lat)) %>% 
  filter(lat_floor %in% top_patches$lat_floor) %>% 
  group_by(haulid) %>% 
  mutate(dens = sum(number_at_length)) %>% # get total no. fish in each haul, of any size
  group_by(year, lat_floor) %>% 
  summarise(mean_dens = mean(dens)) %>%  # get mean density (all sizes) / haul for the patch*year combo 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_test_dens <- dat_test %>% 
  mutate(lat_floor = floor(lat)) %>% 
  filter(lat_floor %in% top_patches$lat_floor) %>% 
  group_by(haulid) %>% 
  mutate(dens = sum(number_at_length)) %>% # get total no. fish in each haul, of any size
  group_by(year, lat_floor) %>% 
  summarise(mean_dens = mean(dens)) %>%  # get mean density (all sizes) / haul for the patch*year combo 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

# get patch area 
patchdat <- dat %>% 
  select(lat, lon) %>%
  distinct() %>%
  mutate(lat_floor = floor(lat)) %>% 
  group_by(lat_floor) %>% 
  summarise(max_lon=max(lon),
            min_lon=min(lon)) %>% 
  rowwise() %>% 
  mutate(lon_dist = distGeo(p1=c(max_lon, lat_floor+0.5), p2=c(min_lon, lat_floor+0.5))/1000, # get distance between the furthest longitudes in km, at the midpoint of the lat band 
         patch_area_km2 = lon_dist * 111) %>%  # 1 degree latitude = 111 km 
  select(lat_floor, patch_area_km2) %>% 
  filter(lat_floor %in% top_patches$lat_floor) %>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_train_sbt <- dat %>%   
  mutate(lat_floor = floor(lat)) %>% 
  group_by(lat_floor, year) %>% 
  summarise(sbt = mean(btemp, na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(lat_floor %in% top_patches$lat_floor)%>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_test_sbt <- dat_test %>%   
  mutate(lat_floor = floor(lat)) %>% 
  group_by(lat_floor, year) %>% 
  summarise(sbt = mean(btemp, na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(lat_floor %in% top_patches$lat_floor)%>% 
  mutate(patch = as.integer(as.factor(lat_floor))) %>% 
  filter(!is.na(sbt))

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
#f = 0.334
#z = exp(-m-f)
age_at_maturity = 2
t0=-.2
cv= 0.2 # guess
min_age = 0
max_age = 20

# now that we have the max_age, fill in f for years above 7 (since the f for age=7 is really for 7+)
older_ages <- expand_grid(age=seq(max(dat_f_age_prep$age)+1, max_age, 1), year= unique(dat_f_age_prep$year)) %>% 
  left_join(dat_f_age_prep %>% filter(age==max(age)) %>% select(year, f))

dat_f_age_prep %<>% bind_rows(older_ages)

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

# get time dimension
years <- sort(unique(dat_train_lengths$year)) 
years_proj <- sort(unique(dat_test_lengths$year))
ny <- length(years)
ny_proj <- length(years_proj)
dat_f_age <- dat_f_age_prep %>% 
  filter(year %in% years) 

dat_f_age_proj <- dat_f_age_prep %>% 
  filter(year %in% years_proj) %>%
  bind_rows(dat_f_age %>% filter(year==max(year))) # need final year of training data to initialize projection

#get other dimensions
patches <- sort(unique(dat_train_lengths$lat_floor))
patcharea <- patchdat %>% arrange(patch) %>% pull(patch_area_km2)
np = length(patches) 

lbins <- unique(length_at_age_key$length_bin)
# lbins <- sort(unique(dat_train_lengths$length))
n_lbins <- length(lbins) 

n_ages <- nrow(l_at_a_mat)

# now that years are defined above, convert them into indices in the datasets
# be sure all these dataframes have exactly the same year range! 

dat_train_dens$year = as.integer(as.factor(dat_train_dens$year))
dat_test_dens$year = as.integer(as.factor(dat_test_dens$year))
dat_train_lengths$year = as.integer(as.factor(dat_train_lengths$year))
#dat_test_lengths$year = as.integer(as.factor(dat_test_lengths$year))
dat_test_sbt$year= as.integer(as.factor(dat_test_sbt$year))
dat_train_sbt$year= as.integer(as.factor(dat_train_sbt$year))
dat_f_age$year = as.integer(as.factor(dat_f_age$year))
dat_f_age_proj$year = as.integer(as.factor(dat_f_age_proj$year))

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
    tmp2 <- dat_train_dens %>% filter(patch==p, year==y) %>% 
      left_join(patchdat, by = c("lat_floor","patch"))%>% 
      mutate(mean_dens = mean_dens * patch_area_km2)
    dens[p,y] <- tmp2$mean_dens
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
    tmp4 <- dat_f_age %>% filter(age==a, year==y) 
    f[a+1,y] <- tmp4$f # add 1 because matrix indexing starts at 1 not 0
  }
}
f_proj <- array(NA, dim=c(n_ages,(ny_proj+1)))
for(a in min_age:max_age){
  for(y in 1:(ny_proj+1)){
    tmp5 <- dat_f_age_proj %>% filter(age==a, year==y) 
    f_proj[a+1,y] <- tmp5$f # add 1 because matrix indexing starts at 1 not 0
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
  patcharea = patcharea,
  l_at_a_key = l_at_a_mat,
  do_dirichlet = 1,
  eval_l_comps = 0, # evaluate length composition data? 0=no, 1=yes
  T_dep_mortality = 0, # 0=off, 1=on
  T_dep_recruitment = 1 # think carefully before making more than one of the temperature dependencies true
  
)

# only run 1 iteration with Fixed_param
warmups <- 0
total_iterations <- 1
max_treedepth <-  10
n_chains <-  1
n_cores <- 1
stan_model_fit <- stan(file = here::here("src","process_sdm.stan"), # check that it's the right model!
                       data = stan_data,
                       algorithm = "Fixed_param",
                       chains = n_chains,
                       warmup = warmups,
                       init = list(list(log_mean_recruits = log(1000000),
                                        theta_d = 1,
                                        p_length_50_sel = 0.25)),
                       iter = total_iterations,
                       cores = n_cores,
                       refresh = 250,
                       control = list(max_treedepth = max_treedepth,
                                      adapt_delta = 0.85)
)


# ok, so the two things are are fitting to are abund_p_y and n_p_l_y
# When we fit a model, we generate the "model" version of these, somewhat confusingly named now
# and dens_p_y_hat and n_p_l_y_hat. 
# So, the idea now is to take the values of dens_p_y_hat and n_p_l_y_hat, and then convert them into data in place of
# abund_p_y and n_p_l_y, and then fit that model

dim(stan_data$abund_p_y)

str(stan_data$abund_p_y)

# note taht parameters are fixed at their starting guesses

theta_d <-  extract(stan_model_fit, "theta_d")

theta_d$theta_d

log_mean_recruits <-  extract(stan_model_fit, "log_mean_recruits")

log_mean_recruits$log_mean_recruits

p_length_50_sel <-  extract(stan_model_fit, "p_length_50_sel")

p_length_50_sel$p_length_50_sel

# now need to wrangle this into a patch x year matrix
new_abund_p_y <-  extract(stan_model_fit, "dens_p_y_hat")
# because life is annoying, this is a list with dimensions iteration, patch, year. because we only used 1 interation (warmup = 0, iter = 1), we don'tn eed to worry about iterations, so just need to collapse this list back down to something more usable. 
# 
# Could use tidybayes for this but sticking with base for now just to show how it's really working

new_abund_p_y <- new_abund_p_y$dens_p_y_hat[1,,] # first iteration, all patches, all years. Note that this would work even if you had more than 1 iteration since with Fixed_params every iteration is the same

plot(new_abund_p_y[4,])

dim(new_abund_p_y) <- dim(stan_data$abund_p_y) # set the dimensions to be the same

lines(new_abund_p_y[4,]) # OK that checks out


# repeat for length comps

dim(stan_data$n_p_l_y)


new_n_p_l_y <-  extract(stan_model_fit, "n_p_l_y_hat") # this one will be a pain

dim(new_n_p_l_y$n_p_l_y_hat)

new_n_p_l_y <- new_n_p_l_y$n_p_l_y_hat[1,,,] # first iteration, all years, all patches, all lengths

dim(new_n_p_l_y)

# sigh. can't just do dim(new_n_p_l_y) <- dim(stan_data$n_p_l_y) since that doesn't preserve the order of things correctly


temp <- array(NA, dim = dim(stan_data$n_p_l_y))

# going to be lazy and just do this like this for now since i can't be bothered to think of a more elegant strategy. 

for( i in 1:stan_data$ny_train){
  
  
  temp[,,i] <- as.integer(round(new_n_p_l_y[i,,] * 5e9)) # oh right, these have to be integers, and for some reason these are crazy small numbers, so just multipldying by something crazy big to make the rounding not set things to 0
  
}

plot(round(new_n_p_l_y[5,6,] * 5e9)) #So, need to convert this to the right dimentions

lines(temp[6,,5]) 
# yay. 

new_n_p_l_y <- temp

fixed_param_stan_data <- stan_data

fixed_param_stan_data$n_p_l_y <- new_n_p_l_y

fixed_param_stan_data$abund_p_y <- new_abund_p_y

# now fit to the fixed_param generated thingy

warmups <- 5000
total_iterations <- 10000
max_treedepth <-  10
n_chains <-  1
n_cores <- 1

fixed_param_stan_model_fit <- stan(file = here::here("src","process_sdm.stan"), # check that it's the right model!
                       data = fixed_param_stan_data,
                       chains = n_chains,
                       warmup = warmups,
                       init = list(list(log_mean_recruits = log(1000000),
                                        theta_d = 1,
                                        p_length_50_sel = 0.25)),
                       iter = total_iterations,
                       cores = n_cores,
                       refresh = 250,
                       control = list(max_treedepth = max_treedepth,
                                      adapt_delta = 0.85)
)

p_length_50_sel <-  extract(fixed_param_stan_model_fit, "p_length_50_sel")

p_length_50_sel$p_length_50_sel

hist(p_length_50_sel$p_length_50_sel)

# now go through and check the fits, if everything is working right the model should fit the generated data perfectly, and should reciver all the parameters etc. If not, it means there are weird pathologies in the model, likely related to multiple stable states