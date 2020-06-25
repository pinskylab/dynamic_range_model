# what's the species overlap between the NOAA trawl data and GlobTherm?

library(tidyverse)
library(here)
library(stringr)
library(lubridate)
library(data.table)
here <- here::here

OApath <- "~/github/OceanAdapt_9815545/"

load(paste0(OApath,"data_clean/all-regions-full.RData"))

# https://datadryad.org/stash/dataset/doi:10.5061/dryad.1cv08
thermdat <- read_csv(here("survey-data","GlobalTherm_upload_02_11_17.csv"))

sharedspp <- thermdat %>% 
  mutate(spp = paste0(Genus, " ", Species)) %>%
  pull(spp) %>%
  intersect(unique(dat$spp))

sharedsummary <- dat %>% 
  filter(spp %in% sharedspp) %>%
  group_by(region, spp) %>%
  summarise(n=n()) 

write_csv(sharedsummary, here("processed-data","NOAA_global_therm_shared_spp.csv"))
