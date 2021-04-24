library(tidyverse)
library(VAST)
library(here)
# unlink(here("Kmeans_extrapolation-2000.RData")) # delete old spatial files if changing regions / spatial extent or resolution
# devtools::install_version("Matrix", version = "1.2.8") # downgrade version of Matrix to match TMB
DF <- as.data.frame(readRDS("/Users/afh/github/range-edge-niches/processed-data/neus_catch_rates.rds")) %>% 
  filter(Sci=="Paralichthys dentatus") %>% 
  mutate(Year=as.numeric(Year))

# starting with the same settings as for edge thermal niche paper, except purpose=index2. Look into overdispersion options in make_settings

Data_Geostat = data.frame( "spp"=DF[,"Sci"], "Year"=DF[,"Year"], "Catch_KG"=DF[,"Wt"], "AreaSwept_km2"=0.01, "Vessel"=0, "Lat"=DF[,"Lat"], "Lon"=DF[,"Long"] )

Region="Northwest_Atlantic"

Version=get_latest_version("VAST_v13_0_0")

n_x = 100 # number of knots in the VAST mesh 
fine_scale = FALSE # spatial variables interpolated between knots, as opposed to assigned the nearest knot 

# predictors used in the model: setting any to 0 removes them for the model structure. we did this for spatial predictors to force the model to attribute variation to spatiotemporal predictors; otherwise it sometimes predicts perfectly static ranges over time (i.e., entirely spatial variation)
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

OverdispersionConfig = c("Eta1"=0, "Eta2"=0) # turning off vessel-level catchability factors 
ObsModel = c(2,1) # specifies functional form of encounter probabilities--here, lognormal

Options =  c("Calculate_Range"=1) # turn on range calculations

strata.limits <- data.frame('STRATA'="All_areas")


# bind together all the model settings we specified earlier
settings = make_settings( "n_x"=n_x, "Region"=Region, purpose="index2", bias.correct=FALSE, use_anisotropy=TRUE,
                          "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig,
                          "Options"=Options,"ObsModel"=ObsModel, "Version"=Version  )


# derived objects that are constant among all models in a region
Year_Set <- sort(unique(DF$Year))

fit = try(
  fit_model("settings"=settings, "Lat_i"=Data_Geostat[,'Lat'],
            "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'],
            "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'],
            "a_i"=Data_Geostat[,'AreaSwept_km2'], 
            "getReportCovariance"=FALSE, "getJointPrecision"=TRUE,
            lower=-Inf, upper=Inf,
            test_fit = FALSE,
            fine_scale=fine_scale,
            anisotropy=TRUE,
            Use_REML=TRUE,
            getsd=TRUE,
            newtonsteps=1
  ))

plot(fit)