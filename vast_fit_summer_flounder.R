#####
## Fitting VAST model to summer flounder data for comparison with DRM 
#####


# Overview ----------------------------------------------------------------
# The objective of this script is to fit an annaul time step VAST model to summer flounder data from the Norhteast Shelf Large Marine Ecosystem, focusing on validating predictive performance within latitudinal bands and comparing with an alternative size-based dynamic range model. To accomplish this objective, the script walks through X different stages -- from initial data prep, through to making predictions from a fitted VAST model. 


# Preliminaries -----------------------------------------------------------
library(devtools)
library(VAST)
library(FishStatsUtils)
library(tidyverse)
library(lubridate)
library(sf)
library(raster)
library(here) 
library(splines)
library(patchwork)
library(akima)
library(tictoc)

# Functions -- sourced from aallyn/TargetsSDM. This includes a lot of helpful "VAST" functions.
source_url("https://raw.githubusercontent.com/aallyn/TargetsSDM/main/R/vast_functions.R")

# Initial data prep: trawl and sample data --------------------------------
# This is mostly taken from "prep_summer_flounder_data.R, with the change that for the VAST model, we are going to want to have the total biomass (e.g., numbers density) for summer flounder at each of the tow locations. This will then get multiplied by the "area" of each of the patches (e.g., latitude bins) to get the total estimated abundance for that time/patch.

# To get the data, I went to the link here and then downloaded the zipped folder. I then grabbed the dat_exploded.rds file, the neus_Survdat.RData file and the neus_SVSPP.RData file. I added a raw-data folder and then gitignored it.
dat_exploded<- readRDS(here("raw-data/dat_exploded.rds"))
load(here("raw-data/neus_Survdat.RData")) # Biomass data
load(here("raw-data/neus_SVSPP.RData")) # Taxanomy data

spp_of_interest<- c("Paralichthys dentatus")
reg_of_interest<- c("Northeast US Fall")

dat_exploded_neus<- dat_exploded %>% 
  filter(region == reg_of_interest) 

max_yr<- 2016 # 2017 is missing, for some reason 
min_yr<- 1972 # early years have some issues; the newly filtered OA data starts here anyway
# to explore more the changes in samplng over time, make annual maps of haul locations, or a tile plot of year vs. lat band 
forecast_yr_1 <- max_yr - 9

# creating a separate haul info dataframe to get the date and btemp
hauldat<- survdat %>% 
  # create a haulid for joining
  mutate(haulid = paste(formatC(CRUISE6, width=6, flag=0), formatC(STATION, width=3, flag=0), formatC(STRATUM, width=4, flag=0), sep='-')) %>% 
  dplyr::select(haulid, BOTTEMP, YEAR, EST_TOWDATE, LAT, LON) %>% 
  distinct() %>%
  rename("btemp"=BOTTEMP,
         "year"=YEAR,
         "date"=EST_TOWDATE,
         "lat"=LAT,
         "lon"=LON)

# tidy biomass data
bio_flounder_prep<- survdat %>% 
  # create a haulid for joining with dat.exploded
  mutate(haulid = paste(formatC(CRUISE6, width=6, flag=0), formatC(STATION, width=3, flag=0), formatC(STRATUM, width=4, flag=0), sep='-')) %>% 
  dplyr::select(haulid, SVSPP, LENGTH, NUMLEN, BIOMASS) %>% 
  left_join(spp, by="SVSPP") %>% # get species names from species codes
  mutate(spp = str_to_sentence(SCINAME)) %>% # change to format of sppOfInt
  dplyr::select(spp, haulid, LENGTH, NUMLEN, BIOMASS) %>%
  filter(!is.na(LENGTH),
         spp == spp_of_interest,
         LENGTH>4 # get rid of a single fish at 4cm -- the next smallest is 9cm
  ) %>% 
  rename("length"=LENGTH,
         "number_at_length"=NUMLEN) %>%
  group_by(., haulid, spp) %>% # Group by haulid, SVSPP
  summarize_at(vars(BIOMASS), sum, na.rm = TRUE) %>%
  rename("biomass" = BIOMASS) # Calculate the total biomass across all length/num length samples

# Need to create biomass df that still includes zeroes?
# important to use haulids from dat_exploded_neus because it has been cleaned
bio_flounder<- expand.grid(haulid = unique(dat_exploded_neus$haulid)) %>% # get full factorial of every haul and biomass value
  mutate(spp = spp_of_interest) %>% 
  # left_joining to use only the hauls in dat_exploded_neus
  left_join(bio_flounder_prep, by = c('haulid', 'spp')) %>% 
  mutate(biomass = replace_na(biomass, 0)) %>%
  left_join(hauldat)

# Double check that...
length(unique(dat_exploded_neus$haulid)) == length(unique(bio_flounder$haulid))

# Looks good to go. Drop NA bottom temps
bio_flounder<- bio_flounder %>%
  drop_na(btemp)

# Filter into training/testing and then write out these files
flounder_train<- bio_flounder %>% 
  filter(year >= min_yr,
         year < forecast_yr_1)

flounder_test<- bio_flounder %>% 
  filter(year >= forecast_yr_1,
         year <= max_yr)
flounder_test$biomass<- 0 # Making sure these observations won't influence things...

write_csv(flounder_train, here("processed-data","flounder_biomass_fall_training.csv"))
write_csv(flounder_test, here("processed-data","flounder_biomass_fall_testing.csv"))

# VAST data --------------------------------------
# First, a sample dataframe. With the VAST framework, we can have all the data (training and testing) and then switch Pred_TF to "1" for the "testing" data years. This means they will not be used in the likelihood, but we will predict to them.
flounder_train$Pred_TF<- 0
flounder_test$Pred_TF<- 1
flounder_all<- bind_rows(flounder_train, flounder_test)
#flounder_all<- bind_rows(flounder_train)
vast_samp_dat<- data.frame("Year" = flounder_all$year, "Lat" = flounder_all$lat, "Lon" = flounder_all$lon, "Biomass" = flounder_all$biomass, "Swept" = rep(1, nrow(flounder_all)), "Pred_TF" = flounder_all$Pred_TF)

# Next, a covariate dataframe and then rescaling bottom temperature
vast_cov_dat<- data.frame("Year" = flounder_all$year, "btemp" = flounder_all$btemp, "Lat" = flounder_all$lat, "Lon" = flounder_all$lon) 
vast_cov_dat$btemp<- as.numeric(scale(vast_cov_dat$btemp), center = TRUE, scale = TRUE)

# VAST model objectives and settings --------------------------------------
# Extrapolation grid and Strata limits. Here, we want to essentially arrive at latitude patch estimates. 
# Grabbing the grid for northwest Atlantic in FishStatsUtils
utils::data(northwest_atlantic_grid, package = "FishStatsUtils")
Data_Extrap<- northwest_atlantic_grid
str(Data_Extrap)

# Looking at this quickly...
plot(Data_Extrap$Lon, Data_Extrap$Lat)

# Compare with the data...
points(flounder_all$lon, flounder_all$lat, col = "red")

# Okay, clearly want to do some work here. First is the crazy outlier in Chesapeake Bay and then the extrapolation grid also extends into area now surveyed more commonly by NOAA SEFSC. For the Chesapeake one...
#identify(flounder_all$lon, flounder_all$lat)
plot(Data_Extrap$Lon, Data_Extrap$Lat)
points(flounder_all$lon, flounder_all$lat, col = "red")
points(flounder_all$lon[922], flounder_all$lat[922], col = "green")

flounder_all<- flounder_all[-922,]
points(flounder_all$lon, flounder_all$lat, col = "blue")

# Remove this obs from sample and covariate data
vast_samp_dat<- vast_samp_dat[-922,]
vast_cov_dat<- vast_cov_dat[-922,]

# Next, work on the extrapolation grid. Go based on min latitude of observations and a little buffer (0.25 degrees)
Data_Extrap<- Data_Extrap %>%
  filter(., Lat >= floor(min(flounder_all$lat)))
points(Data_Extrap$Lon, Data_Extrap$Lat, col = "green")

# Now, we need to have a new column that has the "STRATA" of each observation. In this case, we want those strata to be the latitude patches. From `analyze_summer_flounder` grabbing the different "lat" categories...
vast_extrap_grid<- Data_Extrap %>%
  mutate(STRATA = factor(paste0("Lat_", floor(Lat)), levels = paste0("Lat_", seq(from = min(unique(floor(Data_Extrap$Lat))), to = max(unique(floor(Data_Extrap$Lat))))))) %>%
  rename("Area_km2" = Area_in_survey_km2) %>%
  dplyr::select(., Lon, Lat, Area_km2, STRATA)

# Finally, strata.limits dataframe
#strata_use<- strata_use<- data.frame("STRATA" = as.character(unique(vast_extrap_grid$STRATA)))
strata_use<- data.frame( "STRATA" = paste0("Lat_", seq(from = min(unique(floor(Data_Extrap$Lat))), to = max(unique(floor(Data_Extrap$Lat))))),
                         "south_border" = seq(from = min(unique(floor(Data_Extrap$Lat))), to = max(unique(floor(Data_Extrap$Lat)))),
                         "north_border" = seq(from = min(unique(floor(Data_Extrap$Lat)))+0.99, to = max(unique(floor(Data_Extrap$Lat)))+0.99))

vast_extrap_info<- make_extrapolation_info(Region = "User", strata.limits = strata_use, input_grid = vast_extrap_grid, grid_dim_km = c(25, 25))
colSums(vast_extrap_info$a_el)

# Quick check...
ggplot() +
  geom_point(data = vast_extrap_grid, aes(x = Lon, y = Lat, color = STRATA))

# Field and rho configuration settings. Field config sets up the spatial/spatio-temporal components and how many factors should be estimated. Here, we are going to turn both of those "on". Rho config sets up autoregressive structure on intercepts and spatio-temporal components. Given the interest here on forecasting, going to try AR1 process for both the intercepts and spatio-temporal variability. Will make adjustments as needed given data constraints and convergence issues. 
field_config<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
rho_config<- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 2, "Epsilon2" = 2)

vast_settings<- make_settings(n_x = 400, Region = "User", purpose = "index2", FieldConfig = field_config, RhoConfig = rho_config, ObsModel = c(2, 1), OverdispersionConfig = c(0, 0), bias.correct = TRUE, knot_method = "samples", Options = c("Calculate_Range" = TRUE), strata.limits = strata_use)

# Fitting VAST model ------------------------------------------------------
nice_category_names<- "Summer flounder all"
# Model formula
hab_formula<- ~ bs(btemp, degree = 3, intercept = FALSE) 

# Build VAST model
vast_build<- fit_model("settings" = vast_settings, "Method" = vast_settings$Method, "input_grid" = vast_extrap_grid, "Lat_i" = vast_samp_dat[, 'Lat'], "Lon_i" = vast_samp_dat[, 'Lon'], "t_i" = vast_samp_dat[, 'Year'], "c_i" = rep(0, nrow(vast_samp_dat)), "b_i" = vast_samp_dat[, 'Biomass'], "a_i" = vast_samp_dat[, 'Swept'], "PredTF_i" = vast_samp_dat[, 'Pred_TF'], "covariate_data" = vast_cov_dat, "X1_formula" = hab_formula, "X2_formula" = hab_formula, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = TRUE)

# Looks good. Fit it.
vast_fit<- fit_model("settings" = vast_settings, "Method" = vast_settings$Method, "input_grid" = vast_extrap_grid, "Data_Extrap" = vast_build$extrapolation_list$Data_Extrap, "Lat_i" = vast_samp_dat[, 'Lat'], "Lon_i" = vast_samp_dat[, 'Lon'], "t_i" = vast_samp_dat[, 'Year'], "c_i" = rep(0, nrow(vast_samp_dat)), "b_i" = vast_samp_dat[, 'Biomass'], "a_i" = vast_samp_dat[, 'Swept'], "PredTF_i" = vast_samp_dat[, 'Pred_TF'], "covariate_data" = vast_cov_dat, "X1_formula" = hab_formula, "X2_formula" = hab_formula, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = TRUE, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = TRUE)

# Summarizing VAST model predictions --------------------------------------
## Covariate effects
vast_habcovs_effs<- get_vast_covariate_effects(vast_fit = vast_fit, params_plot = c("btemp"), params_plot_levels = 100, effects_pad_values = c(), nice_category_names = nice_category_names, out_dir = here("results/vast"))

# Warnings are...interesting. I think this is because of how the sampling is done to generate the plots -- basically some sampled values are beyond the boundary knots (e.g., bottom temperature values). 
vast_habcovs_plot<- plot_vast_covariate_effects(vast_covariate_effects = vast_habcovs_effs, vast_fit = vast_fit, nice_category_names = nice_category_names, out_dir = here("results/vast"))
vast_habcovs_plot # Seems to suggest preference for warmer waters to determine general distribution, and then density (biomass) within that area shows more of a "humped" distribution and preference for intermediate btemps

## Stratified biomass indices -- I think this is where the comparisons would be for the patches?
vast_bio<- get_vast_index_timeseries(vast_fit = vast_fit, all_times = unique(vast_samp_dat$Year), nice_category_names = nice_category_names, index_scale = c("raw"), out_dir = here("results/vast"))
summary(vast_bio)

# Plot...
color_pal<- rev(viridis::viridis(n = length(unique(vast_bio$Index_Region))))
plot_vast_index_timeseries(index_res_df = vast_bio, index_scale = "raw", nice_category_names = nice_category_names, nice_xlab = "Year", nice_ylab= "Biomass index (metric tons)", paneling = "none", color_pal = color_pal, out_dir = here("results/vast"))

# Tile plot similar to AF's for the DRM
patch_bio_tile_plot<- ggplot(data = vast_bio, aes(x = Year, y = Index_Region, fill = Index_Estimate)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Biomass index (metric tons)") +
  theme_bw() +
  scale_x_continuous(name = "Year", breaks = seq(from = min(vast_bio$Year), to = max(vast_bio$Year), by = 4), expand = c(0, 0)) +
  scale_y_discrete(name = "Latitiude patch", expand = c(0, 0))
ggsave(paste0(here("results/vast/"), nice_category_names, " patch biomass tile.jpg"), patch_bio_tile_plot)

# A land shapefile and setting lon/lat limits
land_use<- st_read("~/GitHub/TargetsSDM/data/supporting/land_shapefile/ne_50m_land.shp")
xlim_use<- c(-78, -65)
ylim_use<- c(34, 46)

vast_cog_plot<- vast_plot_cog(vast_fit = vast_fit, all_times = unique(vast_samp_dat$Year), summarize = TRUE, land_sf = land_use, xlim = xlim_use, ylim = ylim_use, nice_category_names = nice_category_names, land_color = "#d9d9d9", color_pal = NULL, out_dir = here("results/vast"))
vast_cog_plot

# Extracting and plotting predicted density at knot locations, smoothed over a regular grid
# Need a regional grid...
region_shape<- st_read("~/Box/RES_Data/Shapefiles/BottomTrawlStrata/BTS_Strata.shp") %>%
  st_transform(., crs = 4326) %>%
  st_union()
vast_density_plot<- vast_fit_plot_spatial(vast_fit = vast_fit, spatial_var = "D_gct", nice_category_names = nice_category_names, mask = region_shape, all_times = as.character(unique(vast_samp_dat$Year)), plot_times = NULL, land_sf = land_use, xlim = xlim_use, ylim = ylim_use, panel_or_gif = "panel", out_dir = here("results/vast"), land_color = "#d9d9d9", lab_lat = 37, lab_lon = -68, panel_cols = 7, panel_rows = 8)
vast_density_plot
