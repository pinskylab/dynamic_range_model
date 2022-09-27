set.seed(42)
library(tidyverse)
library(tidybayes)
library(Cairo)
library(here)
library(magrittr)
library(rstan)
library(Matrix)
library(rstanarm)

rstan_options(javascript = FALSE, auto_write = TRUE)

funs <- list.files("functions")
sapply(funs, function(x)
  source(file.path("functions", x)))

ctrl_file <- read_csv("control_file.csv")


ctrl_file <-  ctrl_file %>%
  slice(2) %>% 
  mutate(fits = pmap(
    list(
      run_name = id,
      do_dirichlet = do_dirichlet,
      eval_l_comps = eval_l_comps,
      T_dep_movement = T_dep_movement,
      spawner_recruit_relationship = spawner_recruit_relationship,
      process_error_toggle = process_error_toggle,
      exp_yn = exp_yn,
      known_f = known_f
    ),
    fit_drm,
    warmup = 1000,
    iter = 2000,
    chains = 1,
    cores = 1,
    run_forecast = 1
  ))


# process results ---------------------------------------------------------

tmp <- ctrl_file$fits[[1]]

rstan::check_hmc_diagnostics(tmp)


load(here("processed-data","stan_data_prep.Rdata"))


# stan_data <- list(
#   np=np,
#   n_ages=n_ages,
#   ny_train=ny,
#   ny_proj=ny_proj,
#   n_lbins=n_lbins,
#   n_p_l_y = len,
#   abund_p_y = dens,
#   sbt = sbt,
#   sbt_proj=sbt_proj,
#   m=m,
#   f=f,
#   f_proj=f_proj,
#   k=k,
#   loo=loo,
#   t0=t0,
#   cv=cv,
#   length_50_sel_guess=length_50_sel_guess,
#   n_lbins = n_lbins, 
#   age_sel = age_sel,
#   bin_mids=bin_mids,
#   sel_100=sel_100,
#   age_at_maturity = age_at_maturity,
#   l_at_a_key = l_at_a_mat,
#   wt_at_age = wt_at_age,
#   do_dirichlet = do_dirichlet,
#   eval_l_comps = eval_l_comps, # evaluate length composition data? 0=no, 1=yes
#   T_dep_mortality = T_dep_mortality, # CURRENTLY NOT REALLY WORKING
#   T_dep_recruitment = T_dep_recruitment, # think carefully before making more than one of the temperature dependencies true,
#   T_dep_movement = T_dep_movement,
#   spawner_recruit_relationship = spawner_recruit_relationship, 
#   run_forecast=run_forecast,
#   exp_yn = exp_yn,
#   process_error_toggle = process_error_toggle,
#   number_quantiles = number_quantiles,
#   quantiles_calc = quantiles_calc
# )

abund_p_y <-  dat_train_dens %>%
  mutate(abundance = mean_dens * meanpatcharea) 

abund_p_y_hat <- tidybayes::spread_draws(tmp, dens_p_y_hat[patch,year])


abundance_v_time <- abund_p_y_hat %>% 
  ggplot(aes(year, dens_p_y_hat)) + 
  stat_lineribbon() +
  geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch, scales = "free_y") +
  labs(x="Year",y="Abundance") + 
  scale_fill_brewer()

abund_p_y <-  dat_test_dens %>%
  mutate(abundance = mean_dens * meanpatcharea) 

observed_abund_posterior_predictive <- tidybayes::spread_draws(tmp, dens_p_y_obs_proj[patch,year])

abund_posterior_predictive <- tidybayes::spread_draws(tmp, dens_p_y_proj[patch,year])


observed_abundance_forecast <- observed_abund_posterior_predictive %>% 
  ggplot(aes(year, dens_p_y_obs_proj)) + 
  stat_lineribbon() +
  geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch) +
  labs(x="Year",y="Abundance") + 
  scale_fill_brewer()

abundance_forecast <- abund_posterior_predictive %>% 
  ggplot(aes(year, dens_p_y_proj)) + 
  stat_lineribbon() +
  geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch) +
  labs(x="Year",y="Abundance") + 
  scale_fill_brewer()






