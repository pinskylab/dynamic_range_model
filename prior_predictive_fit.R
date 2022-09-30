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

load(here("processed-data","stan_data_prep.Rdata"))

ctrl_file <- read_csv("control_file.csv")

sample_run <- ctrl_file$id[1]

sample_fit <- read_rds(here("results",sample_run,"stan_model_fit.rds"))


sigma_obs <- mean(extract(sample_fit, "sigma_obs")[[1]])

sample_abund_p_y <- (extract(sample_fit, "dens_p_y_hat")[[1]][1,,]) # get a random draw from posterior of 'abundance'

sample_n_p_l_y <- round(extract(sample_fit, "n_p_l_y_hat")[[1]][1,,,]) # get a random draw from posterior of length comps
year = 12
patch = 5
before <- sample_n_p_l_y[year,patch,]

# ugh. n_p_l_y_hat is different dimensions than n_p_l_y. reshaping to match inputs
being_lazy <- len

no_error <- sample_abund_p_y

for (y in 1:ny){
  
  for (p in 1:np){
    
    sample_abund_p_y[p,y] <- exp(rnorm(1, log(sample_abund_p_y[p,y] + 1e-6), sigma_obs)) # manually add observation error since no posterior predictive yet, change this once new fits are in
    
    being_lazy[p,,y] <- sample_n_p_l_y[y,p,]
    
  }
}

sample_n_p_l_y <- being_lazy

plot(before)

after <- sample_n_p_l_y[patch,,year]

lines(after)

drm_name <-  "process_sdm"

do_dirichlet <- 1

eval_l_comps <-  1

T_dep_mortality = 0

T_dep_recruitment = 0

T_dep_movement = 1

spawner_recruit_relationship = 1

run_forecast = 1

process_error_toggle = 1

exp_yn = 0

warmup = 1000

iter = 2000

max_treedepth =  10

chains =  1

refresh = 10

cores = 1

drm_name = "process_sdm"

number_quantiles = 3

quantiles_calc = c(0.05, 0.5, 0.95)



stan_data <- list(
  np=np,
  n_ages=n_ages,
  ny_train=ny,
  ny_proj=ny_proj,
  n_lbins=n_lbins,
  n_p_l_y = sample_n_p_l_y,
  abund_p_y = sample_abund_p_y,
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
  process_error_toggle = process_error_toggle,
  number_quantiles = number_quantiles,
  quantiles_calc = quantiles_calc,
  sigma_obs_cv = 0.1,
  h = 0.8
)

diagnostic_fit <- stan(file = here::here("src",paste0(drm_name,".stan")), # check that it's the right model!
                       data = stan_data,
                       chains = chains,
                       warmup = warmup,
                       iter = iter,
                       cores = cores,
                       refresh = refresh,
                       control = list(max_treedepth = max_treedepth,
                                      adapt_delta = 0.85),
                       init = lapply(1:chains, function(x) list(Topt = jitter(12,4),
                                                                log_r0 = jitter(10,5),
                                                                beta_obs = jitter(1e-6,4),
                                                                beta_obs_int = jitter(-10,2)))
)




dens_p_y <- sample_abund_p_y %>% 
  as.data.frame() %>% 
  mutate(patch = 1:nrow(.)) %>% 
  pivot_longer(-patch, names_to = "year", values_to = "density",
               names_prefix = "V", names_transform = list(year = as.integer))

dens_p_y_hat <- tidybayes::spread_draws(diagnostic_fit, dens_p_y_hat[patch,year], ndraws = 2)


abundance_v_time <- dens_p_y_hat %>% 
  ggplot(aes(year, dens_p_y_hat)) + 
  stat_lineribbon() +
  geom_point(data = dens_p_y, aes(year, density), color = "red") +
  facet_wrap(~patch, scales = "free_y") +
  labs(x="Year",y="Density") + 
  scale_fill_brewer()


n_p_l_y <- data.table::as.data.table(sample_n_p_l_y) %>% 
  rename(year = V3, patch = V1, length = V2) %>% 
  group_by(year, patch) %>% 
  mutate(plength = value / sum(value)) %>% 
  ungroup() %>% 
  filter(year == min(year) | year == max(year))
  

n_p_l_y_hat <- rstan::extract(diagnostic_fit,"n_p_l_y_hat")[[1]]

sigma_obs_hat <-  rstan::extract(diagnostic_fit,"sigma_obs")[[1]]

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

# abundance_v_time