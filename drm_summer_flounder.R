# ADD FLAG TO TELL ANNOTATE/AMAREL TO RUN WITH R

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
library(parallel)

load(here("processed-data","stan_data_prep.Rdata"))

# SET UP FUNCTION SO IT TAKES ONE ARGUMENT THAT IS A CHARACTER STRING 
# OR ONE ARGUMENT FOR EACH DATA FILE THAT NEEDS TO LOAD? PLUS THE MODEL NAME
# THEN SPECIFY MODELS TO RUN IN BASH 
# LOOK FOR GUIDES TO ARGUMENTS TO R SCRIPTS FROM BASH

# set up different models 
process_error = list(
  do_dirichlet = 1,
  eval_l_comps = 0, 
  T_dep_mortality = 0,
  T_dep_recruitment = 0, 
  T_dep_movement = 0,
  spawner_recruit_relationship = 0,
  run_forecast=1,
  process_error_toggle = 1,
  exp_yn = 0
)

process_error_l_comps = list(
  do_dirichlet = 1,
  eval_l_comps = 1, 
  T_dep_mortality = 0,
  T_dep_recruitment = 0, 
  T_dep_movement = 0,
  spawner_recruit_relationship = 0,
  run_forecast=1,
  process_error_toggle = 1,
  exp_yn = 0
)

run_names <- c("process_error",
               "process_error_l_comps")

run_specs <- list(process_error = process_error, 
                    process_error_l_comps = process_error_l_comps)

# run models in parallel 
run_drm <- function(run_spec){
  
  rstan_specs(javascript=FALSE, auto_write =TRUE)
  
  run_name <- names(run_spec)
  results_path <- file.path("results",run_name)
  
  
  if (!dir.exists(results_path)){
    dir.create(results_path, recursive = TRUE)
  }
  
  #############
  # make model decisions and prep for model 
  #############
  
  # unpack listed model specs into environment 
  list2env(run_spec, .GlobalEnv)
  
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
    eval_l_comps = eval_l_comps,
    T_dep_mortality = T_dep_mortality, 
    T_dep_recruitment = T_dep_recruitment, 
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
  
  warmups <- 100
  total_iterations <- 200
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
  
}


# implement function in parallel 
cl <- parallel::makeCluster(4)
parallel::parLapply(cl,
                   run_specs,
                    run_drm)
parallel::stopCluster(cl)
