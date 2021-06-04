# library(tidyverse)
library(VAST)
# library(here)
# unlink(here("Kmeans_extrapolation-2000.RData")) # delete old spatial files if changing regions / spatial extent or resolution
# devtools::install_version("Matrix", version = "1.2.8") # downgrade version of Matrix to match TMB


#########
# VAST with size structure (https://github.com/James-Thorson-NOAA/VAST/wiki/Expand-age-and-length-composition)
#########

# this is a quick and dirty way to pull in length comps for summer flounder, drawn from process_survey_data.Rmd, but be sure that those methods/data aren't updated later and left out of here! uncomment this part to set up data
# library(tidyverse)
# library(here)
# OApath <- "~/github/OceanAdapt_9815545/"
# load(paste0(OApath,"data_raw/neus_Survdat.RData"))
# load(paste0(OApath, "data_raw/neus_SVSPP.RData"))
# load(paste0(OApath,"data_clean/dat_exploded.RData"))
# 
# dat_exploded_neus <- dat.exploded %>% 
#   filter(region %in% c("Northeast US Spring", "Northeast US Fall"))
# 
# neus_len <- survdat %>% 
#   # create a haulid for joining with dat.exploded
#   mutate(haulid = paste(formatC(CRUISE6, width=6, flag=0), formatC(STATION, width=3, flag=0), formatC(STRATUM, width=4, flag=0), sep='-')) %>% 
#   select(haulid, SVSPP, LENGTH, NUMLEN) %>% 
#   filter(haulid %in% dat_exploded_neus$haulid) %>% # pare down to hauls of interest 
#   left_join(spp, by="SVSPP") %>% # get species names from species codes
#   mutate(spp = str_to_sentence(SCINAME)) %>% # change to format of sppOfInt
#   select(spp, haulid, LENGTH, NUMLEN) %>%
#   filter(!is.na(LENGTH))
# 
# fluke_neus_len_agg <- neus_len %>% filter(spp=="Paralichthys dentatus")  %>% 
#   mutate(Length_bin = ifelse(LENGTH<18, "0-18cm", ifelse(LENGTH<27,"18-27cm","27cm-and-up"))) %>% 
#   group_by(haulid, spp, Length_bin) %>% 
#   summarise(Counts = sum(NUMLEN))  
# 
# DF <- expand_grid(haulid=unique(dat_exploded_neus$haulid),  Length_bin=unique(fluke_neus_len_agg$Length_bin),
#                                     spp="Paralichthys dentatus") %>%  
#   left_join(fluke_neus_len_agg, by=c("spp","haulid","Length_bin")) %>% 
#   left_join(dat_exploded_neus, by=c("spp","haulid")) %>% 
#   mutate(Counts = replace_na(Counts, 0),
#          AreaSwept_km2 = 0.01,
#          year = as.numeric(year)) 
# write_csv(DF, here("processed-data","fluke_for_VAST_length_comps.csv"))

DF <- read.csv(paste0(getwd(),"/processed-data/fluke_for_VAST_length_comps.csv"))
DF$Length_bin_ID <- as.numeric(as.factor(DF$Length_bin))-1 # note that this works without specifying the order of length bins because it just happens to correspond to how R orders the character strings of the length bins I used. for other cases, double-check that the IDs are in the right order!

# Look into overdispersion options in make_settings later 

Region="Northwest_Atlantic"

Version=get_latest_version("VAST_v13_0_0")

n_x = 100 # number of knots in the VAST mesh 
fine_scale = FALSE # spatial variables interpolated between knots, as opposed to assigned the nearest knot 

Options =  c("Calculate_Range"=1) # turn on range calculations

strata.limits <- list('All_areas' = 1:1e5) # see Issue https://github.com/James-Thorson-NOAA/VAST/issues/302
FieldConfig = c("Omega1"=0, # spatial variation in first linear predictor (encounter probability or zero-inflation)
                "Epsilon1"=1, # spatiotemporal variation in same
                "Omega2"=0, # spatial variation in second linear predictor (catch rates or count-data intensity)
                "Epsilon2"=1 # spatiotemporal variation in same
)

# specify model structure; note that when models do not converge, we sequentially reduce the number of parameters being estimated by setting some of these parameters to 1
RhoConfig = c("Beta1"=4, # "autoregressive"; intercept estimated as a fixed effect
              "Beta2"=4, # "autoregressive"; intercept estimated as a fixed effect
              "Epsilon1"=4, # "autoregressive";  autocorrelation in temporal variation  estimated as a fixed effect
              "Epsilon2"=4 # "autoregressive";  autocorrelation in temporal variation  estimated as a fixed effect
)

# bind together all the model settings we specified earlier
settings = make_settings( "n_x"=n_x, "Region"=Region, purpose="index2", bias.correct=FALSE, use_anisotropy=TRUE, "strata.limits"=strata.limits,
                          
                          "Options"=Options,"Version"=Version,
                          "FieldConfig"=FieldConfig,
                          "RhoConfig"=RhoConfig,
                          ObsModel=c(4,2) # from ?make_data: "Poisson-link delta-model, but fixing encounter probability=1 for any year where all samples encounter the species and encounter probability=0 for any year where no samples encounter the species" -- NOT SURE I ACTUALLY INPUT THIS CORRECTLY THOUGH! 
                          )


# derived objects that are constant among all models in a region
Year_Set <- sort(unique(DF$Year))

# Run model
fit = fit_model( settings = settings,
                 Lat_i = DF[,'lat'],
                 Lon_i = DF[,'lon'],
                 t_i = DF[,'year'],
      #           c_i = as.numeric(DF[,'Length_bin'])-1, # this isn't in fit_model's documentation anymore, so maybe deprecated?
                 b_i = DF[,'Counts'],
                 c_iz = matrix(DF[,'Length_bin_ID'], ncol = 1), # https://github.com/James-Thorson-NOAA/VAST/issues/255
                 a_i = DF[,'AreaSwept_km2'],
      #           Npool = 40, # see https://github.com/James-Thorson-NOAA/VAST/issues/220
                 newtonsteps = 1,
                 test_fit = FALSE
      )

#########
# "regular VAST" (no stage structure)
#########
# DF <- as.data.frame(readRDS("/Users/afh/github/range-edge-niches/processed-data/neus_catch_rates.rds")) %>% 
#   filter(Sci=="Paralichthys dentatus") %>% 
#   mutate(Year=as.numeric(Year))
# DF <- read.csv(paste0(getwd(),"/processed-data/fluke_for_VAST.csv"))

# Data_Geostat = data.frame( "spp"=DF[,"Sci"], "Year"=DF[,"Year"], "Catch_KG"=DF[,"Wt"], "AreaSwept_km2"=0.01, "Vessel"=0, "Lat"=DF[,"Lat"], "Lon"=DF[,"Long"] )
# 
# Region="Northwest_Atlantic"
# 
# Version=get_latest_version("VAST_v13_0_0")
# 
# n_x = 100 # number of knots in the VAST mesh 
# fine_scale = FALSE # spatial variables interpolated between knots, as opposed to assigned the nearest knot 
# 
# # predictors used in the model: setting any to 0 removes them for the model structure. we did this for spatial predictors to force the model to attribute variation to spatiotemporal predictors; otherwise it sometimes predicts perfectly static ranges over time (i.e., entirely spatial variation)
# FieldConfig = c("Omega1"=0, # spatial variation in first linear predictor (encounter probability or zero-inflation) 
#                 "Epsilon1"=1, # spatiotemporal variation in same
#                 "Omega2"=0, # spatial variation in second linear predictor (catch rates or count-data intensity) 
#                 "Epsilon2"=1 # spatiotemporal variation in same
# )
# 
# # specify model structure; note that when models do not converge, we sequentially reduce the number of parameters being estimated by setting some of these parameters to 1
# RhoConfig = c("Beta1"=4, # "autoregressive"; intercept estimated as a fixed effect
#               "Beta2"=4, # "autoregressive"; intercept estimated as a fixed effect
#               "Epsilon1"=4, # "autoregressive";  autocorrelation in temporal variation  estimated as a fixed effect
#               "Epsilon2"=4 # "autoregressive";  autocorrelation in temporal variation  estimated as a fixed effect
# ) 
# 
# OverdispersionConfig = c("Eta1"=0, "Eta2"=0) # turning off vessel-level catchability factors 
# ObsModel = c(2,1) # specifies functional form of encounter probabilities--here, lognormal
# 
# Options =  c("Calculate_Range"=1) # turn on range calculations
# 
# strata.limits <- list('All_areas' = 1:1e5) # see Issue https://github.com/James-Thorson-NOAA/VAST/issues/302
# 
# # bind together all the model settings we specified earlier
# settings = make_settings( "n_x"=n_x, "Region"=Region, purpose="index2", bias.correct=FALSE, use_anisotropy=TRUE, "strata.limits"=strata.limits,
#                           "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig,
#                           "Options"=Options,"ObsModel"=ObsModel, "Version"=Version  )
# 
# 
# # derived objects that are constant among all models in a region
# Year_Set <- sort(unique(DF$Year))
# 
# fit = try(
#   fit_model("settings"=settings, "Lat_i"=Data_Geostat[,'Lat'],
#             "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'],
#             "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'],
#             "a_i"=Data_Geostat[,'AreaSwept_km2'], 
#             "getReportCovariance"=FALSE, "getJointPrecision"=TRUE,
#             lower=-Inf, upper=Inf,
#             test_fit = FALSE,
#             fine_scale=fine_scale,
#             anisotropy=TRUE,
#             Use_REML=TRUE,
#             getsd=TRUE,
#             newtonsteps=1
#   ))
# 
# plot(fit) # this isn't working right now: "Error in (function (plot_set = 3, Obj = NULL, PlotDF, Sdreport = NULL,  : object 'Nyears' not found"
# 
# # get density out 
# get_density <- function(fit.model, Years2Include){
#   
#   # Get object
#   obj<- fit.model$tmb_list[["Obj"]]
#   
#   # Extract
#   report<- obj$report()
#   
#   # Slots worth describing -- density at each knot? Can project this after to get east/northing, or could pull east/northing out of the fit.model$extrapolation_list$Data_Extrap object...
#   dens.df<- data.frame("lon" = fit.model$spatial_list$latlon_g[,"Lon"],
#                        "lat" = fit.model$spatial_list$latlon_g[,"Lat"],
#                        # "density" = report$D_gcy # THIS IS FROM AN OLDER VERSION OF VAST
#                        "density" = report$Index_gctl # I'm not 1000% sure this is biomass at each knot. double-check! 
#                        )
#   
#   dens.df <- reshape2::melt(dens.df, id.vars=c('lon','lat'))
#   colnames(dens.df) <- c('lon','lat','year',
#                          'density')
#   
#   dens.df$year <- as.numeric(gsub("density.","", dens.df$year))
#   dens.df$year <- fit$year_labels[dens.df$year] # convert to real years 
#   
#   dens.df <- dens.df[dens.df$year %in% Years2Include,] # keep only years with real data 
#   
#   return(dens.df)
# }
# 
# flukedat <- get_density(fit, Year_Set)
# 
# flukeround <- flukedat %>% 
#   mutate(latround = floor(lat)) %>% 
#   group_by(latround, year) %>% 
#   summarise(biomass = sum(density))
