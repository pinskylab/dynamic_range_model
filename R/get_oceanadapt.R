# get data from NOAA trawl surveys
# this assumes that oceanadapt has been downloaded to the same directory as this repository, ~/github, and is stored in a directory called "OceanAdapt_9815545"

library(tidyverse)

# get NEUS and keep length frequency

load("~/github/OceanAdapt_9815545/data_raw/neus_Survdat.RData")
load("~/github/OceanAdapt_9815545/data_raw/neus_SVSPP.RData")

neus_spp <- spp %>%
  # remove some columns from spp data
  select(-ITISSPP, -COMNAME, -AUTHOR) %>% 
  mutate(SCINAME = as.character(SCINAME))

neusdat <- survdat %>% 
  left_join(neus_spp, by = "SVSPP")

saveRDS(neusdat, file=paste0(getwd(), "/processed-data/raw_neus_trawl.rds"))

# get SEUS and keep length frequency

seus_catch <- read_csv(unz("~/github/OceanAdapt_9815545/data_raw/seus_catch.csv.zip", "seus_catch.csv"), col_types = cols(.default = col_character())) %>% 
  # remove symbols
  mutate_all(list(~str_replace(., "=", ""))) %>% 
  mutate_all(list(~str_replace(., '"', ''))) %>% 
  mutate_all(list(~str_replace(., '"', '')))

seus_catch <- type_convert(seus_catch, col_types = cols(
  PROJECTNAME = col_character(),
  PROJECTAGENCY = col_character(),
  DATE = col_character(),
  EVENTNAME = col_character(),
  COLLECTIONNUMBER = col_character(),
  VESSELNAME = col_character(),
  GEARNAME = col_character(),
  GEARCODE = col_character(),
  SPECIESCODE = col_character(),
  MRRI_CODE = col_character(),
  SPECIESSCIENTIFICNAME = col_character(),
  SPECIESCOMMONNAME = col_character(),
  NUMBERTOTAL = col_integer(),
  SPECIESTOTALWEIGHT = col_double(),
  SPECIESSUBWEIGHT = col_double(),
  SPECIESWGTPROCESSED = col_character(),
  WEIGHTMETHODDESC = col_character(),
  ORGWTUNITS = col_character(),
  EFFORT = col_character(),
  CATCHSUBSAMPLED = col_logical(),
  CATCHWEIGHT = col_double(),
  CATCHSUBWEIGHT = col_double(),
  TIMESTART = col_character(),
  DURATION = col_integer(),
  TOWTYPETEXT = col_character(),
  LOCATION = col_character(),
  REGION = col_character(),
  DEPTHZONE = col_character(),
  ACCSPGRIDCODE = col_character(),
  STATIONCODE = col_character(),
  EVENTTYPEDESCRIPTION = col_character(),
  TEMPSURFACE = col_double(),
  TEMPBOTTOM = col_double(),
  SALINITYSURFACE = col_double(),
  SALINITYBOTTOM = col_double(),
  SDO = col_character(),
  BDO = col_character(),
  TEMPAIR = col_double(),
  LATITUDESTART = col_double(),
  LATITUDEEND = col_double(),
  LONGITUDESTART = col_double(),
  LONGITUDEEND = col_double(),
  SPECSTATUSDESCRIPTION = col_character(),
  LASTUPDATED = col_character()
))

saveRDS(seus_catch, file=paste0(getwd(), "/processed-data/raw_seus_trawl.rds"))
