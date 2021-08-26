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

dat_train_sbt <- dat %>%   # make temperature matrix 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(lat_floor, year) %>% 
  summarise(sbt = mean(btemp, na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(lat_floor %in% top_patches$lat_floor)%>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

# set fixed parameters from stock assessment
Loo = 83.6
k = 0.14
# mortality
m = 0.25
f = 0.334
z = exp(-m-f)
age_at_maturity = 2

# get time dimension
years <- sort(unique(dat_train_lengths$year)) 
ny <- length(years)
ny_proj <- 10

#get other dimensions
patches <- sort(unique(dat_train_lengths$lat_floor))
np = length(patches) 
lbins <- sort(unique(dat_train_lengths$length))
n_lbins <- length(lbins) 

# make length-at-age key 
length_at_age_key_dat <- spasm::generate_length_at_age_key(
  min_age = min(ages), 
  max_age = max(ages), 
  cv = 0.2, # "magic number" - revisit this sometime (use LIME as a guide?)
  k = k, 
  linf = Loo,
  t0 = 0, # is this right?
  time_step = 1,
  linf_buffer = 1 # doing this just to shoehorn it into the right size for the model; can make it bigger if we need to
)

# now that years are defined above, convert them into indices in the datasets
dat_train_dens$year = as.integer(as.factor(dat_train_dens$year))
dat_train_lengths$year = as.integer(as.factor(dat_train_lengths$year))

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



