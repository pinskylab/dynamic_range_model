
# load packages
library(tidyverse)
library(here)
library(stringr)
library(lubridate)
here <- here::here

# download from: https://zenodo.org/badge/latestdoi/29789533
OApath <- "~/github/OceanAdapt-update2020/"

load(paste0(OApath,"data_clean/dat_exploded.rds")) # get zero-inflated survey data
# for some reason this line isn't working for me and I have to import it manually, not sure why

load(paste0(OApath,"data_raw/neus_Survdat.RData")) # get length data
load(paste0(OApath, "data_raw/neus_SVSPP.RData")) # get taxonomy

spp_of_interest <- c("Squalus acanthias")
reg_of_interest <- c("Northeast US Fall")


max_yr <- 2016 # 2017 is missing, for some reason 
min_yr <- 1972 # early years have some issues; the newly filtered OA data starts here anyway
# to explore more the changes in samplng over time, make annual maps of haul locations, or a tile plot of year vs. lat band 
forecast_yr_1 <- max_yr - 9

dat_exploded_neus <- dat_exploded %>% 
  filter(region %in% reg_of_interest,
         year >= min_yr,
         year <= max_yr) 

explore_data <- FALSE
if(explore_data==TRUE){
  
  # get years each stratum was sampled 
  stratum_pres <- dat_exploded_neus %>% 
    select(stratum, year) %>% 
    distinct() %>% 
    mutate(sampled=TRUE)
  
  # get all combinations of stratum*year
  stratum_all <- expand_grid(stratum=unique(dat_exploded_neus$stratum), year=unique(dat_exploded_neus$year)) %>% 
    left_join(stratum_pres) %>% 
    mutate(sampled = replace_na(sampled, FALSE)) %>% 
    left_join(dat_exploded %>% select(lat, lon, stratum) %>% distinct()) # get lat/lon 
  
  ggplot(stratum_all, aes(x=year, y=lat, fill=sampled, color=sampled)) +
    geom_tile() # looks like basically all strata are sampled every year 
  
  ggplot(stratum_pres, aes(x=stratum)) + 
    geom_histogram(stat="count") +
    theme(axis.text.x = element_text(angle=90))
  
  ggplot(stratum_all, aes(x=year, y=stratum, fill=sampled, color=sampled)) +
    geom_tile()  
  
  load("~/github/OceanAdapt_9815545/data_clean/dat_exploded.Rdata")
  
  old_dat_exploded <- dat.exploded
  rm(dat.exploded)
  
  old_dat_exploded_neus <- old_dat_exploded %>% 
    filter(region==reg_of_interest) 
  
  old_stratum_pres <- old_dat_exploded_neus %>% 
    select(stratum, year) %>% 
    distinct() %>% 
    mutate(sampled=TRUE)
  
  ggplot(old_stratum_pres, aes(x=stratum)) + 
    geom_histogram(stat="count") +
    theme(axis.text.x = element_text(angle=90))
  
  # get all combinations of stratum*year
  old_stratum_all <- expand_grid(stratum=unique(old_dat_exploded_neus$stratum), year=unique(old_dat_exploded_neus$year)) %>% 
    left_join(old_stratum_pres) %>% 
    mutate(sampled = replace_na(sampled, FALSE)) %>% 
    left_join(old_dat_exploded %>% select(lat, lon, stratum) %>% distinct()) # get lat/lon 
  
  
  ggplot(old_stratum_all, aes(x=year, y=stratum, fill=sampled, color=sampled)) +
    geom_tile()  
  
  
}



# creating a separate haul info dataframe to get the date and btemp
hauldat <- survdat %>% 
  # create a haulid for joining
  mutate(haulid = paste(formatC(CRUISE6, width=6, flag=0), formatC(STATION, width=3, flag=0), formatC(STRATUM, width=4, flag=0), sep='-')) %>% 
  filter(haulid %in% dat_exploded_neus$haulid) %>% 
  select(haulid, BOTTEMP, YEAR, EST_TOWDATE, LAT, LON) %>% 
  distinct() %>%
  rename("btemp"=BOTTEMP,
         "year"=YEAR,
         "date"=EST_TOWDATE,
         "lat"=LAT,
         "lon"=LON)

# 
# tidy length data
len_prep <- survdat %>%
  # create a haulid for joining with dat.exploded
  mutate(haulid = paste(formatC(CRUISE6, width=6, flag=0), formatC(STATION, width=3, flag=0), formatC(STRATUM, width=4, flag=0), sep='-')) %>%
  filter(haulid %in% dat_exploded_neus$haulid) %>%
  select(haulid, SVSPP, LENGTH, NUMLEN, CATCHSEX) %>% 
  left_join(spp, by="SVSPP") %>% # get species names from species codes
  mutate(spp = str_to_sentence(SCINAME)) %>% # change to format of sppOfInt
  select(spp, haulid, LENGTH, NUMLEN, CATCHSEX) %>%
  filter(!is.na(LENGTH),
         spp == spp_of_interest,
         LENGTH > 1 # get rid of a single record of 1cm dogfish
  ) %>% 
  rename("length"=LENGTH,
         "number_at_length"=NUMLEN) %>% 
  group_by_at(vars(-CATCHSEX, -number_at_length)) %>% 
  summarise(number_at_length = sum(number_at_length)) %>%  # pool counts across all sexes
  ungroup()

# need to create length df that still includes zeroes
# important to use haulids from dat_exploded_neus because it has been cleaned
len <- expand.grid(haulid=unique(dat_exploded_neus$haulid), length=seq(min(len_prep$length), max(len_prep$length), 1)) %>% # get full factorial of every haul * length bin
  mutate(spp = spp_of_interest) %>% 
  # left_joining to use only the hauls in dat_exploded_neus
  left_join(len_prep, by=c('length','haulid','spp')) %>% 
  left_join(dat_exploded %>% select(haulid, region) %>% distinct(), by="haulid") %>% # get region column back
  mutate(number_at_length = replace_na(number_at_length, 0)) %>%  # fill in absences with true zeroes
  left_join(hauldat)

explore_len <- TRUE
if(explore_len <- TRUE){
  len %>% ggplot(aes(x=length, y=number_at_length)) +
    geom_line() + 
    facet_wrap(~year)
}

len_train <- len %>% 
  filter(year >= min_yr,
         year < forecast_yr_1)

len_test <- len %>% 
  filter(year >= forecast_yr_1,
         year <= max_yr)

write_csv(len_train, here("processed-data","spiny_dogfish_catch_at_length_fall_training.csv"))
write_csv(len_test, here("processed-data","spiny_dogfish_catch_at_length_fall_testing.csv"))
