library(spasm)
library(tidyverse)
library(FishLife)
library(spasm)
library(ggridges)
library(sampling)
library(gganimate)
library(rstan)
library(ggdist)

fish <-
  create_fish(
    scientific_name = "Atractoscion nobilis",
    query_fishlife = T,
    mat_mode = "length",
    time_step = 1,
    sigma_r = 0.4,
    price = 10,
    price_cv = 0,
    price_ac = 0,
    price_slope = 0,
    steepness = 0.9,
    r0 = 10000,
    rec_ac = .4,
    density_movement_modifier = 0.1,
    adult_movement = 20,
    larval_movement = 3,
    density_dependence_form = 3
  )



np = 2

fleet <- create_fleet(
  fish = fish,
  cost_cv =  0,
  cost_ac = 0,
  cost_slope = 0,
  q = .1,
  q_cv = 0,
  q_ac = .7,
  q_slope = 0,
  fleet_model = "constant-effort",
  target_catch = 200,
  sigma_effort = 0,
  length_50_sel = .5 * fish$length_50_mature,
  initial_effort = (fish$m / .1) * np,
  profit_lags =  1,
  beta = 1,
  max_cr_ratio = 0.8,
  max_perc_change_f = 2,
  effort_allocation = 'profit-gravity',
  b_ref_oa = .9
)



a <- Sys.time()
sim <- spasm::sim_fishery(
  fish = fish,
  fleet = fleet,
  manager = create_manager(mpa_size = 0, year_mpa = 5),
  num_patches = np,
  sim_years = 200,
  burn_years = 1,
  time_step = fish$time_step,
  est_msy = FALSE,
  random_mpas = TRUE,
  min_size = 0.05,
  mpa_habfactor = 1,
  sprinkler = FALSE,
  keep_burn = FALSE,
  tune_costs = FALSE
)
Sys.time() - a
b = plot_spasm(sim)
a = plot_spasm(sim)


sim <- sim %>% 
  filter(year > 175) 


sim$year <- sim$year - min(sim$year) + 1

plot_spasm(sim)


n_at_age <- sim %>%
  select(year,patch, age, numbers, numbers_caught) %>%
  group_by(year,patch) %>%
  nest() %>%
  mutate(
    tmp  = map(data,~ sample_lengths(
      n_at_age = .x,
      cv = 0.1,
      k = fish$vbk,
      linf = fish$linf,
      t0 = fish$t0,
      sample_type = 'catch',
      percent_sampled = .5,
      time_step = 1
    ))
  )

length_samples <- n_at_age %>% 
  select(year,patch,tmp) %>% 
  unnest(cols = tmp) %>% 
  filter(length_bin < 1.25 * fish$linf)
  

length_samples %>%
  # filter(length_bin > 10) %>% 
  group_by(year) %>% 
  mutate(numbers = numbers / sum(numbers)) %>% 
  ggplot() +
  geom_density_ridges(aes(length_bin, y = factor(year),height = numbers),stat = "identity")


length_comps <- length_samples %>% 
  group_by(year, length_bin) %>% 
  summarise(numbers = sum(numbers)) %>% 
  ungroup() %>% 
  filter(length_bin < 120)
  

length_comps %>%
  filter(length_bin > 10) %>% 
  group_by(year) %>% 
  mutate(numbers = numbers / sum(numbers)) %>% 
  ggplot() +
  geom_density_ridges(aes(length_bin, y = factor(year),height = numbers),stat = "identity")


# length_samples %>%
#   ggplot(aes(length_bin, numbers)) +
#   geom_col() + 
#   facet_wrap(~year)

abundance_index <- n_at_age %>% 
  select(year,patch, data) %>% 
  unnest(cols = data) %>% 
  group_by(year, patch) %>% 
  summarise(numbers = sum(numbers_caught)) %>% 
  ungroup()

abundance_index %>% 
  ggplot(aes(year, numbers, color = factor(patch))) + 
  geom_line()
  
total_abundance <- abundance_index %>% 
  group_by(year) %>% 
  summarise(total_abundance = sum(numbers))

save(list = c("fish","length_comps", "total_abundance"),file =  "for_lime.Rdata")
# prepare data for model --------------------------------------------------

min_age <- min(sim$age)

max_age <- max(sim$age)

length_at_age_key <-
  generate_length_at_age_key(
    min_age = min_age,
    max_age = max_age,
    cv = 0.1,
    linf = fish$linf,
    k = fish$vbk,
    t0 = fish$t0,
    time_step = 1,
    linf_buffer = 1.5
  )

length_at_age_key %>% 
  filter(age > 5) %>% 
  ggplot(aes(age, length_bin, fill = p_bin)) + 
  geom_hline(aes(yintercept = fish$linf)) +
  geom_tile() + 
  scale_fill_viridis_c()

l_at_a_mat <- length_at_age_key %>% 
  select(age, length_bin, p_bin) %>% 
  pivot_wider(names_from = length_bin, values_from = p_bin) %>% 
  ungroup() %>% 
  select(-age) %>% 
  as.matrix()
dim(l_at_a_mat)

# nn <- 100 * exp(-.2 * (min_:20))
# 
# t(l_at_a_mat) %*% nn -> b

# get time dimension
years <- sort(unique(abundance_index$year)) 
ny <- length(years)
ny_proj <- 10

#get other dimensions
patches <- sort(unique(length_samples$patch))
np = length(patches) 

lbins <- unique(length_at_age_key$length_bin)
# lbins <- sort(unique(dat_train_lengths$length))
n_lbins <- length(lbins) 

n_ages <- nrow(l_at_a_mat)

# now that years are defined above, convert them into indices in the datasets
# dat_train_dens$year = as.integer(as.factor(dat_train_dens$year))
# dat_train_lengths$year = as.integer(as.factor(dat_train_lengths$year))
# dat_train_sbt$year= as.integer(as.factor(dat_train_sbt$year))

# make matrices/arrays from dfs

# length_samples <- length_samples %>%
#   filter(patch ==1, year %in% c(1,5,6,10,12,15))

len <- array(0, dim = c(np, n_lbins, ny)) 
for(p in 1:np){
  for(l in 1:n_lbins){
    for(y in 1:ny){
      tmp <- length_samples %>% filter(patch==p, round(length_bin)==lbins[l], year==y) 
      if (nrow(tmp) > 0){
        len[p,l,y] <- tmp$numbers
      }
    }
  }
}

sim2 <- sim %>% 
  filter(patch ==1, year == 10)


ages <- array(0, dim = c(np, n_ages, ny)) 
for(p in 1:np){
  for(l in unique(sim$age)){
    for(y in 1:ny){
      tmp <- sim2 %>% filter(patch==p, age == l, year==y) 
      if (nrow(tmp) > 0){
        ages[p,l,y] <- round(tmp$numbers)
      }
    }
  }
}

plot(len[1,,2])

# plot(ages[1,,3])


dens <- array(NA, dim=c(np, ny))
for(p in 1:np){
  for(y in 1:ny){
    tmp2 <- abundance_index %>% filter(patch==p, year==y) 
    dens[p,y] <- rlnorm(1,log(tmp2$numbers), .05)
  }
}

sbt <- array(NA, dim=c(np,ny))
for(p in 1:np){
  for(y in 1:ny){
    sbt[p,y] <- rnorm(1)
  }
}

a <- seq(min_age, max_age)

catches <- sim %>% 
  group_by(patch, year) %>% 
  summarise(catch= sum(numbers_caught)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "year", values_from = "catch") %>% 
  select(-patch) %>% 
  as.matrix()


# fit simulated data ------------------------------------------------------

stan_data <- list(
  np=np,
  n_ages=n_ages,
  ny_train=ny,
  n_lbins=n_lbins,
  n_p_l_y = len,
  n_p_a_y = ages,
  c_p_y = catches,
  abund_p_y = dens,
  sbt = sbt,
  m = fish$m,
  k= fish$vbk,
  loo= fish$linf,
  t0= fish$t0,
  cv= .1,
  length_50_sel_guess=20, # THIS IS A RANDOM GUESS, I can't find help in the stock assessment
  n_lbins = n_lbins, 
  age_sel = 0,
  bin_mids=lbins+0.5, # also not sure if this is the right way to calculate the midpoints
  sel_100 = 3, # not sure if this should be 2 or 3. it's age 2, but it's the third age category because we start at 0, which I think Stan will classify as 3...?
  age_at_maturity = round(fish$age_mature),
  patcharea = rep(1, np),
  l_at_a_key = l_at_a_mat,
  do_dirichlet = 1
)

warmups <- 1000
total_iterations <- 2000
max_treedepth <-  10
n_chains <-  1
n_cores <- 1
stan_model_fit <- stan(file = here::here("src","process_sdm.stan"), # check that it's the right model!
                       data = stan_data,
                       init = list(list(log_mean_recruits = rep(log(100000), np),
                                        theta_d = 1)),
                       chains = n_chains,
                       warmup = warmups,
                       iter = total_iterations,
                       cores = n_cores,
                       refresh = 250,
                       control = list(max_treedepth = max_treedepth,
                                      adapt_delta = 0.85)
)


a = extract(stan_model_fit, "theta_d")

# b = extract(stan_model_fit, "scalar")

hist(a$theta_d)

# rstanarm::launch_shinystan(stan_model_fit)

# assess abundance fits


abund_p_y_hat <- tidybayes::spread_draws(stan_model_fit, dens_p_y_hat[patch,year])

abund_p_y_hat %>% 
  ggplot(aes(year, (dens_p_y_hat ))) + 
  stat_lineribbon() + 
  geom_point(data = abundance_index, aes(year, (numbers)), color = "red") +
  facet_wrap(~patch, scales = "free_y") + 
  scale_x_continuous(name = "Year") + 
  scale_y_continuous(name = "Abundance")

n_p_l_y_hat <- tidybayes::gather_draws(stan_model_fit, n_p_l_y_hat[year,patch,length], n = 500)

p = 1

dat_train_lengths <- length_samples %>% 
  group_by(patch, year) %>% 
  mutate(p_length = numbers / sum(numbers))

n_p_l_y_hat %>% 
  ungroup() %>% 
  filter(patch == p) %>% 
  group_by(patch, year, .iteration) %>% 
  mutate(pvalue = .value / sum(.value)) %>% 
  ggplot(aes(length, pvalue)) + 
  stat_lineribbon() + 
  geom_point(data = dat_train_lengths %>% filter(patch == p), aes(length_bin,p_length), color = "red", alpha = 0.2) +
  facet_wrap(~year, scales = "free_y")


