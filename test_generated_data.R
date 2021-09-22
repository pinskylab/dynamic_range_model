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

funs <- list.files("functions")
sapply(funs, function(x) source(file.path("functions",x)))

rstan_options(javascript=FALSE, auto_write =TRUE)

dat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv"))

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

dat_train_dens <- dat %>% 
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

# set fixed parameters from stock assessment
loo = 83.6
k = 0.14
m = 0.25
f = 0.334
z = exp(-m-f)
age_at_maturity = 2
t0=-.2
cv= 0.2 # guess
min_age = 0
max_age = 20


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
ny <- length(years)
ny_proj <- 10

#get other dimensions
patches <- sort(unique(dat_train_lengths$lat_floor))
np = length(patches) 

lbins <- unique(length_at_age_key$length_bin)
# lbins <- sort(unique(dat_train_lengths$length))
n_lbins <- length(lbins) 

n_ages <- nrow(l_at_a_mat)

# now that years are defined above, convert them into indices in the datasets
dat_train_dens$year = as.integer(as.factor(dat_train_dens$year))
dat_train_lengths$year = as.integer(as.factor(dat_train_lengths$year))
dat_train_sbt$year= as.integer(as.factor(dat_train_sbt$year))

# make matrices/arrays from dfs
len <- array(0, dim = c(np, n_lbins, ny)) 
# for(p in 1:np){
#   for(l in 1:n_lbins){
#     for(y in 1:ny){
#       tmp <- dat_train_lengths %>% filter(patch==p, round(length)==lbins[l], year==y) 
#       if (nrow(tmp) > 0){
#         len[p,l,y] <- tmp$sum_num_at_length
#       }
#     }
#   }
# }

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

a <- seq(min_age, max_age)

check <- a %*% l_at_a_mat
######
# fit model
######

stan_data <- list(
  np=np,
  n_ages=n_ages,
  ny_train=ny,
  n_lbins=n_lbins,
  n_p_l_y = len,
  abund_p_y = dens,
  sbt = sbt,
  m=m,
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
  patcharea = patchdat$patch_area_km2,
  l_at_a_key = l_at_a_mat,
  do_dirichlet = 1
)

warmups <- 1500
total_iterations <- 3000
max_treedepth <-  10
n_chains <-  1
n_cores <- 1
stan_model_fit <- stan(file = here::here("src","process_sdm.stan"), # check that it's the right model!
                       data = stan_data,
                       chains = n_chains,
                       warmup = warmups,
                        iter = total_iterations,
                       cores = n_cores,
                       refresh = 250,
                       control = list(max_treedepth = max_treedepth,
                                      adapt_delta = 0.85)
)

# a = extract(stan_model_fit, "sigma_obs")

# hist(a$sigma_obs)
rstanarm::launch_shinystan(stan_model_fit)

# assess abundance fits

abund_p_y <- dat_train_dens %>%
  left_join(patchdat, by = c("lat_floor", "patch")) %>% 
  group_by(patch, year) %>% 
  summarise(abundance = sum(mean_dens *patch_area_km2)) %>% 
  ungroup()


abund_p_y_hat <- tidybayes::spread_draws(stan_model_fit, dens_p_y_hat[patch,year])

abund_p_y_hat %>% 
  ggplot(aes(year, dens_p_y_hat)) + 
  stat_lineribbon() + 
  geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch, scales = "free_y")

# assess length comp fits

n_p_l_y_hat <- tidybayes::gather_draws(stan_model_fit, n_p_l_y_hat[year,patch,length], n = 500)

dat_train_lengths

p = 5

dat_train_lengths <- dat_train_lengths %>% 
  group_by(patch, year) %>% 
  mutate(p_length = sum_num_at_length / sum(sum_num_at_length))

n_p_l_y_hat %>% 
  ungroup() %>% 
  filter(patch == p) %>% 
  group_by(patch, year, .iteration) %>% 
  mutate(pvalue = .value / sum(.value)) %>% 
  ggplot(aes(length, pvalue)) + 
  stat_lineribbon() + 
  geom_point(data = dat_train_lengths %>% filter(patch == p), aes(length,p_length), color = "red", alpha = 0.2) +
  facet_wrap(~year, scales = "free_y")




stop()
# now fit the generated length comps --------------------------------------

tmp_len <- rstan::extract(stan_model_fit, "n_p_l_y_hat")$n_p_l_y_hat

# make matrices/arrays from dfs
len2 <- array(0, dim = c(np, n_lbins, ny)) 
for(p in 1:np){
  for(y in 1:ny){
    tmp <- round(colMeans(tmp_len[,y,p,]))
    len2[p,,y] <- tmp
  }
}


plot(len2[5,,1])

stan_data_2 <- list(
  np=np,
  n_ages=n_ages,
  ny_train=ny,
  n_lbins=n_lbins,
  n_p_l_y = len2,
  abund_p_y = dens,
  sbt = sbt,
  m=m,
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
  patcharea = patchdat$patch_area_km2,
  l_at_a_key = l_at_a_mat,
  do_dirichlet = 1
)

stan_model_fit_2 <- stan(file = here::here("src","process_sdm.stan"), # check that it's the right model!
                         data = stan_data_2,
                         chains = n_chains,
                         warmup = warmups,
                         iter = total_iterations,
                         cores = n_cores,
                         refresh = 250,
                         control = list(max_treedepth = max_treedepth,
                                        adapt_delta = 0.85)
)

tmp <- rstan::extract(stan_model_fit_2, "sigma_obs")$sigma_obs

tmp2 <- rstan::extract(stan_model_fit_2, "length_50_sel")$length_50_sel


abund_p_y <- dat_train_dens %>%
  left_join(patchdat, by = c("lat_floor", "patch")) %>% 
  group_by(patch, year) %>% 
  summarise(abundance = sum(mean_dens *patch_area_km2)) %>% 
  ungroup()


abund_p_y_hat_2 <- tidybayes::spread_draws(stan_model_fit_2, dens_p_y_hat[patch,year])

abund_p_y_hat_2 %>% 
  ggplot(aes(year, dens_p_y_hat)) + 
  stat_lineribbon() +
  geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch, scales = "free_y")

# assess length comp fits

n_p_l_y_hat_2 <- tidybayes::gather_draws(stan_model_fit_2, n_p_l_y_hat[year,patch,length], n = 500)

n_p_l_y_2 <- reshape2::melt(len2)

names(n_p_l_y_2) <- c("patch","length","year","number")

n_p_l_y_2 <- as_tibble(n_p_l_y_2)


p = 5

n_p_l_y_2 <- n_p_l_y_2 %>% 
  group_by(patch, year) %>% 
  mutate(p_length = number / sum(number))

n_p_l_y_hat_2 %>% 
  ungroup() %>% 
  filter(patch == p) %>% 
  group_by(patch, year, .iteration) %>% 
  mutate(pvalue = .value / sum(.value)) %>% 
  ggplot(aes(length, pvalue)) + 
  stat_lineribbon() + 
  geom_point(data = n_p_l_y_2 %>% filter(patch == p), aes(length,p_length), color = "red", alpha = 0.2) +
  facet_wrap(~year, scales = "free_y")

save(list = c("stan_model_fit","stan_model_fit_2"), file = "flounderfits.Rdata")



