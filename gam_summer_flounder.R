library(mgcv)
library(here)
library(tidyverse)
set.seed(42)

load(here("processed-data","stan_data_prep.Rdata"))
dat <- read_csv(here("processed-data","flounder_catch_for_sdm_fall_training.csv")) %>% 
  mutate(abund_scale = scale(abundance, center=TRUE, scale=TRUE))
dat_abund_mean = mean(dat$abundance)
dat_abund_sd = sd(dat$abundance)

dat_proj <- read_csv(here("processed-data","flounder_catch_for_sdm_fall_testing.csv"))  %>% 
  mutate(abund_scale = scale(abundance, center=TRUE, scale=TRUE))
dat_proj_abund_mean = mean(dat_proj$abundance)
dat_proj_abund_sd = sd(dat_proj$abundance)
# the mean and SD are pretty different of these two which is problematic for scaling

# what fraction of rows have NA temperature values? 
nrow(dat %>% filter(is.na(btemp)))/nrow(dat) # 12.8% 
nrow(dat_proj %>% filter(is.na(btemp)))/nrow(dat_proj) # 12.0%

# is there a pattern to what years they are in? 
quantile((dat %>% filter(is.na(btemp)))$year) # slightly more in the earlier part of the time-series 

# GAM doesn't like NA predictors, need to remove
dat <- drop_na(dat)
dat_proj <- drop_na(dat_proj)

# fit GAMs
gam1 <- gam(abund_scale ~ s(btemp), data=dat, family=poisson(link = "log"))
gam.check(gam1)
gam1 <- gam(abund_scale ~ s(btemp, k=4), data=dat, family=poisson(link = "log"))
gam.check(gam1)
gam1 <- gam(abund_scale ~ s(btemp, k=10), data=dat, family=poisson(link = "log"))
gam.check(gam1) # why does the p-value go up and then down again as you increase k? is it because it's more unlikely that k would approach edf?
# also—despite having run set.seed--some of these have different results the different times I ran them! one had k=3, edf=2.96, p-value=0.41, and the other had the same k and edf but a p-value of 2e-16. why?

plot(gam1)
pred1 <- predict(gam1, newdata=dat_proj)
gam1_proj <- cbind(dat_proj, data.frame('gam_abund'=unname(pred1))) %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(lat_floor) %>% 
  summarise(med_abund_prep = median(gam_abund)) %>% 
  mutate(proj_abund = (med_abund_prep * dat_abund_sd) + dat_abund_mean)
# still getting negative fish! 

# should I be centering and scaling to avoid negatives? can then re-scale back by patch size 

gam2 <- gam(abundance ~ s(btemp) + s(lat) + s(lon) + s(year), data = dat, family=poisson(link="log"))
gam.check(gam2) # noam ross tutorial: "small p-values indicate that residuals are not randomly distributed. This often means there are not enough basis functions"
gam3 <- gam(abundance ~ s(btemp, k=10) + s(lat, k=10) + s(lon, k=10) + s(year, k=4), data = dat, family=poisson(link="log"))
gam.check(gam3) 
