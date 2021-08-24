
# load packages
library(tidyverse)
library(here)
library(stringr)
library(lubridate)
here <- here::here

# load data -- need to download all oceanadapt data from oceanadapt.rutgers.edu and store in this directory first
OApath <- "~/github/OceanAdapt_9815545/"

load(paste0(OApath,"data_clean/dat_exploded.Rdata")) # get zero-inflated survey data
load(paste0(OApath,"data_raw/neus_Survdat.RData")) # get length data
load(paste0(OApath, "data_raw/neus_SVSPP.RData")) # get taxonomy

spp_of_interest <- c("Paralichthys dentatus")
reg_of_interest <- c("Northeast US Fall")
max_yr <- 2015 # avoiding the 2016 data gap
min_yr <- 1972 # early years have some issues
# to explore more the changes in samplng over time, make annual maps of haul locations, or a tile plot of year vs. lat band 
forecast_yr_1 <- max_yr - 9

# creating a separate haul info dataframe to get the date and btemp
hauldat <- survdat %>% 
  # create a haulid for joining
  mutate(haulid = paste(formatC(CRUISE6, width=6, flag=0), formatC(STATION, width=3, flag=0), formatC(STRATUM, width=4, flag=0), sep='-')) %>% 
  select(haulid, BOTTEMP, YEAR, EST_TOWDATE, LAT, LON) %>% 
  distinct() %>%
  rename("btemp"=BOTTEMP,
         "year"=YEAR,
         "date"=EST_TOWDATE,
         "lat"=LAT,
         "lon"=LON)

# get the dat.exploded for just this region and survey
dat_exploded_reg <- dat.exploded %>% 
  filter(region == reg_of_interest) 

# tidy length data
len_flounder_prep <- survdat %>% 
  # create a haulid for joining with dat.exploded
  mutate(haulid = paste(formatC(CRUISE6, width=6, flag=0), formatC(STATION, width=3, flag=0), formatC(STRATUM, width=4, flag=0), sep='-')) %>% 
  select(haulid, SVSPP, LENGTH, NUMLEN) %>% 
  left_join(spp, by="SVSPP") %>% # get species names from species codes
  mutate(spp = str_to_sentence(SCINAME)) %>% # change to format of sppOfInt
  select(spp, haulid, LENGTH, NUMLEN) %>%
  filter(!is.na(LENGTH),
         spp == spp_of_interest,
         LENGTH>4 # get rid of a single fish at 4cm -- the next smallest is 9cm
  ) %>% 
  rename("length"=LENGTH,
         "number_at_length"=NUMLEN)

# need to create length df that still includes zeroes
len_flounder <- expand.grid(haulid=unique(dat_exploded_reg$haulid), length=seq(min(len_flounder_prep$length), max(len_flounder_prep$length), 1)) %>% # get full factorial of every haul * length bin
  mutate(spp = spp_of_interest) %>% 
  left_join(len_flounder_prep, by=c('length','haulid','spp')) %>% 
  mutate(number_at_length = replace_na(number_at_length, 0)) %>%  # fill in absences with true zeroes
  left_join(hauldat)

flounder_train <- len_flounder %>% 
  filter(year >= min_yr,
         year < forecast_yr_1)

flounder_test <- len_flounder %>% 
  filter(year >= forecast_yr_1,
         year <= max_yr)

write_csv(flounder_train, here("processed-data","flounder_catch_at_length_fall_training.csv"))
write_csv(flounder_test, here("processed-data","flounder_catch_at_length_fall_testing.csv"))