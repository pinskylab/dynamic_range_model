#############
# load packages and data
#############

# may not need all of these anymore since splitting out prep_summer_flounder
set.seed(42)
library(tidyverse)
library(tidybayes)
library(Cairo)
library(here)
library(magrittr)
library(rstan)
library(Matrix)
library(rstanarm)

rstan_options(javascript=FALSE, auto_write =TRUE)
load(here("processed-data","stan_data_prep.Rdata"))

run_name <- "process_error_stock_recruit_Tmov"

results_path <- file.path("results",run_name)

if (!dir.exists(results_path)){
  dir.create(results_path, recursive = TRUE)
}


#############
# make model decisions and prep for model 
#############
do_dirichlet = 1
eval_l_comps = 0 # evaluate length composition data? 0=no, 1=yes
T_dep_mortality = 0 # CURRENTLY NOT REALLY WORKING
T_dep_recruitment = 0 # think carefully before making more than one of the temperature dependencies true
T_dep_movement = 1
spawner_recruit_relationship = 1
run_forecast=1
process_error_toggle = 1
exp_yn = 0

# note that many more model decisions are made in the data reshaping in prep_summer_flounder.R!

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
  length_50_sel_guess=length_50_sel_guess,
  n_lbins = n_lbins, 
  age_sel = age_sel,
  bin_mids=bin_mids,
  sel_100=sel_100,
  age_at_maturity = age_at_maturity,
  l_at_a_key = l_at_a_mat,
  wt_at_age = wt_at_age,
  do_dirichlet = do_dirichlet,
  eval_l_comps = eval_l_comps, # evaluate length composition data? 0=no, 1=yes
  T_dep_mortality = T_dep_mortality, # CURRENTLY NOT REALLY WORKING
  T_dep_recruitment = T_dep_recruitment, # think carefully before making more than one of the temperature dependencies true,
  T_dep_movement = T_dep_movement,
  spawner_recruit_relationship = spawner_recruit_relationship, 
  run_forecast=run_forecast,
  exp_yn = exp_yn,
  process_error_toggle = process_error_toggle
)
nums <- 100 * exp(-.2 * (0:(n_ages - 1)))
check <- t(l_at_a_mat) %*% matrix(nums,nrow = n_ages, ncol = 1)
round(sum(nums)) == round(sum(check))
plot(check)

######
# fit model
######

warmups <- 1000
total_iterations <- 2000
max_treedepth <-  10
n_chains <-  4
n_cores <- 4

######################### RUN THE MODEL
stan_model_fit <- stan(file = here::here("src","process_sdm.stan"), # check that it's the right model!
                       data = stan_data,
                       chains = n_chains,
                       warmup = warmups,
                       #     init = list(list(log_mean_recruits = log(1000),
                       #                      theta_d = 1,
                       #                     ssb0=1000000)),
                       iter = total_iterations,
                       cores = n_cores,
                       refresh = 10,
                       control = list(max_treedepth = max_treedepth,
                                      adapt_delta = 0.85),
                       init = lapply(1:n_cores, function(x) list(Topt = jitter(12,4),
                                                                 log_r0 = jitter(10,5),
                                                                 beta_obs = jitter(1e-6,4),
                                                                 beta_obs_int = jitter(-10,2)))
)
write_rds(stan_data, file.path(results_path, "stan_model_data.rds"))
readr::write_rds(stan_model_fit, file.path(results_path,
                                                  "stan_model_fit.rds"))

