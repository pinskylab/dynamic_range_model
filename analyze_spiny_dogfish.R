#############
# load packages and data
#############
# STILL NEED TO FIX:
# TIME VARYING F
# WEIGHT AT AGE
# FIT UPDATED MODEL, NOT OLD ONE

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

dat <- read_csv(here("processed-data","dogfish_catch_at_length_fall_training.csv"))
# dat %<>% filter(length >17)
dat_test <- read_csv(here("processed-data","dogfish_catch_at_length_fall_testing.csv"))

# the f-at-age data starts in 1982; fill in the previous years with the earliest year of data
# dat_f_age_prep <- read_csv(here("processed-data","spiny_dogfish_F_by_age.csv")) %>%
#   rename_with(str_to_lower)
# f_early <- expand_grid(year=seq(1972, 1981, 1), age=unique(dat_f_age_prep$age)) %>% 
#   left_join(dat_f_age_prep %>% filter(year==1982) %>% select(age, f)) 
# dat_f_age_prep <- bind_rows(dat_f_age_prep, f_early)

make_data_plots <- FALSE

if(make_data_plots==TRUE){
  # how much variation is there in length frequency over time?
  gglength <- dat %>% 
    mutate(lat_floor = floor(lat)) %>% 
    filter(between(lat_floor, 37,41)) %>% 
    group_by(length, year, lat_floor) %>% 
    summarise(total = sum(number_at_length)) %>% 
    ungroup() %>% 
    mutate(total = recode(as.numeric(total), `0` = 1e-10)) %>% # adding because there were zeros causing the weights argument below to throw an error
    mutate(year = as.character(year)) %>% 
    ggplot(aes(x=length, y=year)) +
    geom_density_ridges(aes(height=..density..,
                            weight=total),    
                        scale= 4,
                        stat="density") +
    scale_y_discrete(expand = c(0, 0))  + 
    facet_wrap(~lat_floor) + 
    coord_flip()
   gglength  
  ggsave(gglength, filename=here("results","spiny_dogfish_length_freq.png"))
  
  # by patch -- 
  ggpatchlength <- dat %>% 
    mutate(lat_floor = floor(lat)) %>% 
    group_by(length, lat_floor, year) %>% 
    summarise(total = sum(number_at_length)) %>% 
    mutate(total = recode(as.numeric(total), `0` = 1e-10)) %>% # adding because there were zeros causing the weights argument below to throw an error
    mutate(year = as.character(year)) %>% 
    ggplot(aes(x=length, y=year)) +
    geom_density_ridges(aes(height=..density..,
                            weight=total),    
                        scale= 4,
                        stat="density") +
    scale_y_discrete(expand = c(0, 0)) +
    facet_wrap(~lat_floor, scales = "free_y")
   ggpatchlength  
  ggsave(ggpatchlength, filename=here("results","spiny_dogfish_length_freq_by_patch.png"), height=8, width=5, dpi=160)
  # not sure this doesn't have any data in the other patches, maybe they are missing years? 
}

#############
# make model decisions
#############
# the 0s and 1s are for Stan
trim_to_abundant_patches=TRUE
do_dirichlet = 1
eval_l_comps = 0 # evaluate length composition data? 0=no, 1=yes
T_dep_mortality = 0 # CURRENTLY NOT REALLY WORKING
T_dep_recruitment = 0 # think carefully before making more than one of the temperature dependencies true
spawner_recruit_relationship = 0
run_forecast=0
time_varying_f=FALSE

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
    slice(1:2) # JUST FOR DEVELOPING MODEL
  
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
# library(lme4)
# sbt_lme <- lmer(btemp ~ lat + (1|year), data=bind_rows(dat, dat_test))
# sbt_test_fill <- data.frame(
#   lat_floor = c(35, 36, 37),
#   year = rep(2008, 3),
#   patch = c(1,2,3),
#   sbt = c(predict(sbt_lme, newdata=data.frame(lat=35.5, year=2008)),
#           predict(sbt_lme, newdata=data.frame(lat=36.5, year=2008)),
#           predict(sbt_lme, newdata=data.frame(lat=35.5, year=2008))
#   ))
# dat_test_sbt <- bind_rows(dat_test_sbt, sbt_test_fill)

# set fixed parameters from stock assessment (SAW 43)
# NOTE THAT THESE ALL VARY BY SEX FOR SPINY DOGFISH; VALUES BELOW ARE IN BETWEEN MALES AND FEMALES FOR loo, age_at_maturity
loo = 110 # p41
k = 0.1128 # p41
m = 0.092 #p48
age_at_maturity = 12 #p17
t0=-2.552 # p41
cv= 0.2 # guess
min_age = 0
max_age = 50 #p17

# now that we have the max_age, fill in f for years above 7 (since the f for age=7 is really for 7+)
# older_ages <- expand_grid(age=seq(max(dat_f_age_prep$age)+1, max_age, 1), year= unique(dat_f_age_prep$year)) %>% 
#   left_join(dat_f_age_prep %>% filter(age==max(age)) %>% select(year, f))
# 
# dat_f_age_prep %<>% bind_rows(older_ages)

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

# # prep f data
# dat_f_age <- dat_f_age_prep %>% 
#   filter(year %in% years) 
# 
# dat_f_age_proj <- dat_f_age_prep %>% 
#   filter(year %in% years_proj) %>%
#   bind_rows(dat_f_age %>% filter(year==max(year))) # need final year of training data to initialize projection

lbins <- unique(length_at_age_key$length_bin)
# lbins <- sort(unique(dat_train_lengths$length))
n_lbins <- length(lbins) 

n_ages <- nrow(l_at_a_mat)

# now that we have n_ages, calculate weight at age
# THIS IS VERY HARD CODED, AND ALSO NOT PROBABILISTIC
# NEED TO FIX
# wt_at_age_prep <- read_csv(here("processed-data","summer_flounder_wt_at_age.csv")) %>% 
#   filter(!Age %in% seq(7, 10, 1)) %>% 
#   mutate(Age = gsub("over7",7,Age),
#          Age = as.numeric(Age),
#          Age = Age + 1) %>% # start at 1 not 0
#   group_by(Age) %>% 
#   summarise(wt = mean(Wt)) %>% # average over all years 
#   ungroup() %>% 
#   arrange(Age)
# # fill in plus ages
# wt_at_age_add <- data.frame(Age = seq(max(wt_at_age_prep+1),n_ages,1), wt = slice_tail(wt_at_age_prep)$wt)
# wt_at_age <- rbind(wt_at_age_prep, wt_at_age_add)$wt

# now that years are defined above, convert them into indices in the datasets
# be sure all these dataframes have exactly the same year range! 

dat_train_dens$year = as.integer(as.factor(dat_train_dens$year))
dat_test_dens$year = as.integer(as.factor(dat_test_dens$year))
dat_train_lengths$year = as.integer(as.factor(dat_train_lengths$year))
#dat_test_lengths$year = as.integer(as.factor(dat_test_lengths$year))
dat_test_sbt$year= as.integer(as.factor(dat_test_sbt$year))
dat_train_sbt$year= as.integer(as.factor(dat_train_sbt$year))
# dat_f_age$year = as.integer(as.factor(dat_f_age$year))
# dat_f_age_proj$year = as.integer(as.factor(dat_f_age_proj$year))

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
# 
# 
f <- array(NA, dim=c(n_ages,ny))
for(a in min_age:max_age){
  for(y in 1:ny){
    f[a+1,y] <- 0.128 # FIX LATER, RANDOM PULL FROM STOCK ASSESSMENT
  }
}
f_proj <- array(NA, dim=c(n_ages,(ny_proj+1)))
for(a in min_age:max_age){
  for(y in 1:(ny_proj+1)){
    f_proj[a+1,y] <- 0.128# FIX LATER, RANDOM PULL FROM STOCK ASSESSMENT
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
  sel_100 = 1, # REVISIT
  age_at_maturity = age_at_maturity,
  l_at_a_key = l_at_a_mat,
  do_dirichlet = do_dirichlet,
  eval_l_comps = eval_l_comps, # evaluate length composition data? 0=no, 1=yes
  T_dep_mortality = T_dep_mortality, # CURRENTLY NOT REALLY WORKING
  T_dep_recruitment = T_dep_recruitment # think carefully before making more than one of the temperature dependencies true
)

warmups <- 1000
total_iterations <- 2000
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
