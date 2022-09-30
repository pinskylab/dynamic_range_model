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

fit_drms <- TRUE

if (fit_drms){
drm_fits <-  ctrl_file %>%
  filter(id == "v0.1") %>% 
  mutate(fits = pmap(
    list(
      run_name = id,
      do_dirichlet = do_dirichlet,
      eval_l_comps = eval_l_comps,
      T_dep_movement = T_dep_movement,
      spawner_recruit_relationship = spawner_recruit_relationship,
      process_error_toggle = process_error_toggle,
      exp_yn = exp_yn,
      known_f = known_f,
      known_historic_f = known_historic_f
    ),
    fit_drm,
    warmup = 1000,
    iter = 2000,
    chains = 1,
    cores = 1,
    run_forecast = 1
  ))
} else {
  drm_fits <- ctrl_file %>%
    mutate(fits = pmap(list(run_name = id), ~ purrr::safely(readr::read_rds)(
        here("results", .x, "stan_model_fit.rds")
    )))
  
  drm_worked <- map_lgl(map(drm_fits$fits,"error"), is.null)
  
  drm_fits <- drm_fits %>% 
    filter(drm_worked) %>% 
    mutate(fits = map(fits, "result"))
  
}

# process results ---------------------------------------------------------

diagnostic_fit <- drm_fits$fits[[drm_fits$id == "v0.6"]]

rstan::check_hmc_diagnostics(diagnostic_fit)

load(here("processed-data","stan_data_prep.Rdata"))


abund_p_y <-  dat_train_dens %>%
  mutate(abundance = mean_dens * meanpatcharea) 

abund_p_y_hat <- tidybayes::spread_draws(diagnostic_fit, dens_p_y_hat[patch,year])


abundance_v_time <- abund_p_y_hat %>% 
  ggplot(aes(year, dens_p_y_hat)) + 
  stat_lineribbon() +
  geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch, scales = "free_y") +
  labs(x="Year",y="Abundance") + 
  scale_fill_brewer()

abund_p_y <-  dat_test_dens %>%
  mutate(abundance = mean_dens * meanpatcharea) 

observed_abund_posterior_predictive <- tidybayes::spread_draws(diagnostic_fit, dens_p_y_obs_proj[patch,year])

abund_posterior_predictive <- tidybayes::spread_draws(diagnostic_fit, dens_p_y_proj[patch,year])

observed_abundance_forecast <- observed_abund_posterior_predictive %>% 
  ggplot(aes(year, dens_p_y_obs_proj)) + 
  stat_lineribbon() +
  geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch, scales = "free_y") +
  labs(x="Year",y="Abundance") + 
  scale_fill_brewer()

abundance_forecast <- abund_posterior_predictive %>% 
  ggplot(aes(year, dens_p_y_proj)) + 
  stat_lineribbon() +
  geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch, scales = "free_y") +
  labs(x="Year",y="Abundance") + 
  scale_fill_brewer()

n_p_l_y <- data.table::as.data.table(len) %>% 
  rename(year = V3, patch = V1, length = V2) %>% 
  group_by(year, patch) %>% 
  mutate(plength = value / sum(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(year == min(year) | year == max(year))


n_p_l_y_hat <- rstan::extract(diagnostic_fit,"n_p_l_y_hat")[[1]]

sigma_obs_hat <-  rstan::extract(diagnostic_fit,"sigma_obs")[[1]]

theta_d <-  rstan::extract(diagnostic_fit,"sigma_obs")[[1]]


n_p_l_y_hat <- n_p_l_y_hat[sample(1:1000, 100, replace = FALSE),c(1,35),,]

n_p_l_y_hat <- data.table::as.data.table(n_p_l_y_hat) %>% 
  rename(year = V3, patch = V2, length = V1) %>% 
  group_by(year, patch,iterations) %>% 
  mutate(plength = value / sum(value)) %>% 
  ungroup()

n_p_l_y_hat$year[n_p_l_y_hat$year ==  max(n_p_l_y_hat$year)] <-  max(n_p_l_y$year)

n_p_l_y_hat %>% 
  ggplot(aes(length, plength)) +
  stat_lineribbon(size = .1) +
  geom_point(data = n_p_l_y, aes(length, plength), color = "red", alpha = 0.25) +
  facet_grid(patch ~ year, scales = "free_y") + 
  scale_x_continuous(limits = c(0, 50))




