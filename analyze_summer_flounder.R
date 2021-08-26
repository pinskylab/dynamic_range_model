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
library(spasm)
library(geosphere)
rstan_options(javascript=FALSE)

dat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv"))

# how much variation is there in length frequency over time?
library(ggridges)
gglength <- dat %>% 
  group_by(length, year) %>% 
  summarise(total = sum(number_at_length)) %>% 
  mutate(year = as.character(year)) %>% 
  ggplot(aes(x=length, y=year)) +
  geom_density_ridges(aes(height=..density..,
                          weight=total),    
                      scale= 4,
                      stat="density") +
  scale_y_discrete(expand = c(0, 0)) 
gglength  
ggsave(gglength, filename=here("results","summer_flounder_length_freq.png"), height=8, width=3, dpi=160)

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
ggpatchlength  
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
l0=0
cv=0.2 # guess

# get time dimension
years <- sort(unique(dat_train_lengths$year)) 
ny <- length(years)
ny_proj <- 10

#get other dimensions
patches <- sort(unique(dat_train_lengths$lat_floor))
np = length(patches) 
lbins <- sort(unique(dat_train_lengths$length))
n_lbins <- length(lbins) 

# now that years are defined above, convert them into indices in the datasets
dat_train_dens$year = as.integer(as.factor(dat_train_dens$year))
dat_train_lengths$year = as.integer(as.factor(dat_train_lengths$year))
dat_train_sbt$year= as.integer(as.factor(dat_train_sbt$year))

# make matrices/arrays from dfs
len <- array(NA, dim = c(np, n_lbins, ny)) 
for(p in 1:np){
  for(l in 1:n_lbins){
    for(y in 1:ny){
      tmp <- dat_train_lengths %>% filter(patch==p, length==lbins[l], year==y) 
      len[p,l,y] <- tmp$sum_num_at_length
    }
  }
}

dens <- array(NA, dim=c(np, ny))
for(p in 1:np){
  for(y in 1:ny){
    tmp2 <- dat_train_dens %>% filter(patch==p, year==y) %>% 
      left_join(patchdat)%>% 
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
  z=z,
  k=k,
  loo=loo,
  l0=l0,
  cv=cv,
  length_50_sel_guess=20, # THIS IS A RANDOM GUESS, I can't find help in the stock assessment
  n_lbins = n_lbins, 
  age_sel = 0,
  bin_mids=lbins+0.5, # also not sure if this is the right way to calculate the midpoints
  sel_100 = 3, # not sure if this should be 2 or 3. it's age 2, but it's the third age category because we start at 0, which I think Stan will classify as 3...?
  age_at_maturity = age_at_maturity,
  patcharea = patchdat$patch_area_km2
)

warmups <- 2000
total_iterations <- 4000
max_treedepth <-  12
n_chains <-  4
n_cores <- 4 
stan_model_fit <- stan(file = here::here("src","T_dep_rec_age_str_summer_flounder.stan"), # check that it's the right model!
                      data = age_data,
                      chains = 4,
                      warmup = 5000,
                      iter = 10000,
                      cores = 4,
                      refresh = 250,
                      control = list(max_treedepth = max_treedepth,
                                     adapt_delta = 0.95)
)
