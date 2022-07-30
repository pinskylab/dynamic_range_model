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

run_name <- "without_lcomps"

results_path <- file.path("results",run_name)


if (!dir.exists(results_path)){
  dir.create(results_path, recursive = TRUE)
}

rstan_options(javascript=FALSE, auto_write =TRUE)
load(here("processed-data","stan_data_prep.Rdata"))

#############
# make model decisions and prep for model 
#############
do_dirichlet = 1
eval_l_comps = 0 # evaluate length composition data? 0=no, 1=yes
T_dep_mortality = 0 # CURRENTLY NOT REALLY WORKING
T_dep_recruitment = 0 # think carefully before making more than one of the temperature dependencies true
T_dep_movement = 1
spawner_recruit_relationship = 1
run_forecast=0

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
  exp_yn = 0
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

####### IGNORE THIS PART IT DOESN'T REALLY WORK YET
# # exploring T-dep movement model from Jim
# n_g = np
# 
# # Domain characteristics
# lat_g = patches
# Temp_g = sbt
# 
# # Parameters
# diffusion = 0.8^2
# preference_g = -0.1 * (Temp_g - mean(Temp_g))^2
# 
# # Movement operator
# A_gg = ifelse( round(abs(outer(lat_g,lat_g,"-")),2) == round(mean(diff(lat_g)),2), 1, 0 ) # make a matrix that is 1 when the distance between patches equals the minimum possible nonzero distance (i.e., they are adjacent) and 0 otherwise 
# # Diffusion
# diffusion_gg = A_gg * diffusion 
# diag(diffusion_gg) = -1 * colSums(diffusion_gg) # fill in the diagonal with within-patch "diffusion" to balance out diffusion between adjacent patches
# # Taxis
# taxis_gg = array(NA, dim=np, np, ny)
# mrate_gg = array(NA, dim=c(ny, np, np))
# mfraction_gg = array(NA, dim=c(ny, np, np))
# for(y in 1:ny){
#   taxis_gg[,y] = A_gg * outer(preference_g[,y], preference_g[,y], "-") # multiply the adjacency matrix by the preference difference between adjacent patches
# diag(taxis_gg[y]) = -1 * colSums(taxis_gg[y]) # fill in the diagonal with within-patch "taxis" (don't get this either)
# taxis_gg[y] = expm(taxis_gg[y] )# edit as per Dan/Devin/Jim convo
# # Total
# mrate_gg[y] = diffusion_gg + taxis_gg[y] # movement as a sum of diffusion and taxis (can cancel each other out)
# # Annualized
# mfraction_gg[y] = Matrix::expm(mrate_gg[y])
# }
# # Stationary distribution
# stationary_g = eigen(mfraction_gg)$vectors[,1]
# stationary_g = stationary_g / sum(stationary_g)
# 
# #
# matplot( x=lat_g, y=cbind(preference_g-min(preference_g),stationary_g), type="l", col=c("black","blue") )

# n(t+1) = Mrate_gg * n(t)

######################### RUN THE MODEL
stan_model_fit <- stan(file = here::here("src","process_sdm_T_dep_movement.stan"), # check that it's the right model!
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
                       init = lapply(1:n_cores, function(x) list(Topt = jitter(12,4)))
)

readr::write_rds(stan_model_fit, file = file.path(results_path,
                                                  "stan_model_fit.rds"))


# a = rstan::extract(stan_model_fit, "theta_d")
# write_rds(stan_model_fit,"sigh.rds")
# hist(a$sigma_obs)
# rstanarm::launch_shinystan(stan_model_fit)


# plot important parameters 
plot(stan_model_fit, pars=c('sigma_r','sigma_obs','d','width','Topt','alpha','beta_obs','theta_d'))
plot(stan_model_fit, pars=c('sigma_r','sigma_obs','d','alpha','beta_obs','theta_d'))


# assess abundance fits

hist(extract(stan_model_fit, "mean_recruits")$mean_recruits)

abund_p_y <- dat_train_dens %>%
  mutate(abundance = mean_dens * meanpatcharea) 
# left_join(patchdat, by = c("lat_floor", "patch")) %>% 
# group_by(patch, year) %>% 
# summarise(abundance = sum(mean_dens *meanpatcharea)) %>% 
# ungroup()


abund_p_y_hat <- tidybayes::spread_draws(stan_model_fit, dens_p_y_hat[patch,year])

abundance_v_time <- abund_p_y_hat %>% 
  ggplot(aes(year, dens_p_y_hat)) + 
  stat_lineribbon() + 
  geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch, scales = "free_y") +
  labs(x="Year",y="Abundance") + 
  scale_fill_brewer()
abundance_v_time
ggsave(abundance_v_time, filename=here("results","density_v_time_no_length_comps.png"), width=7, height=4)

# assess length comp fits

n_p_l_y_hat <- tidybayes::gather_draws(stan_model_fit, n_p_l_y_hat[year,patch,length], n = 200)

# neff <- tidybayes::gather_draws(stan_model_fit, n_eff[patch,year], n = 500)

#neff <- tidybayes::gather_draws(stan_model_fit, n_eff[patch,year], n = 500)

p = 2

dat_train_lengths <- dat_train_lengths %>% 
  group_by(patch, year) %>% 
  mutate(p_length = sum_num_at_length / sum(sum_num_at_length))

dat_train_lengths %>% 
  group_by(year,patch) %>% 
  summarise(n = sum(sum_num_at_length)) %>% 
  ggplot(aes(year, n, color =factor(patch))) + 
  geom_point()

n_p_l_y_hat %>% 
  ungroup() %>% 
  filter(patch == p, year > 30) %>% 
  group_by(patch, year, .iteration) %>% 
  mutate(pvalue = .value / sum(.value)) %>% 
  ggplot(aes(length, pvalue)) + 
  stat_lineribbon() +
  geom_point(data = dat_train_lengths %>% filter(patch == p, year > 30), aes(length,p_length), color = "red", alpha = 0.2) +
  facet_wrap(~year, scales = "free_y")

# length frequency over time
l_freq_time <- n_p_l_y_hat %>% 
  filter(year > 30) %>% 
  ungroup() %>% 
  ggplot(aes(x=length, y=..density.., weight=.value)) + 
  geom_histogram(bins=50) +
  facet_grid(patch~year)

ggsave(l_freq_time, filename=here("results","length_freq_time_flounder.png"), scale=1.5, width=15, height=4)

# is there a temperature - recruitment relationship? 

## first need to figure out which lengths correspond to age 0 
length_at_age_key %>% 
  ggplot(aes(x=length_bin, y=p_bin)) +
  geom_line() + 
  facet_wrap(~age)

## let's call recruits <5cm

selectivity_at_bin <- gather_draws(stan_model_fit, selectivity_at_bin[n_lbins], n=500)

mean_sel <- selectivity_at_bin %>% 
  group_by(n_lbins) %>% 
  summarise(mean_sel = mean(.value)) %>% 
  rename(length = n_lbins)

# get proportion of recruits adjusted by selectivity
recruits <- n_p_l_y_hat %>% 
  left_join(mean_sel, by="length") %>% 
  mutate(recruit = ifelse(length <=5, "yes", "no"),
         val_post_sel = .value / mean_sel) %>% 
  group_by(recruit, patch, year, .draw) %>% 
  summarise(sumcount = sum(val_post_sel)) %>% 
  ungroup() %>% 
  group_by(patch, year, .draw) %>% 
  mutate(prop_recruit = sumcount / sum(sumcount)) %>% 
  filter(recruit == "yes")

recruits %>% group_by(patch, year) %>% 
  mutate(mean_prop_rec = mean(prop_recruit)) %>% 
  left_join(dat_train_sbt, by=c('patch','year')) %>% 
  ggplot(aes(x=sbt, y=mean_prop_rec, color=year, fill=year)) +
  geom_point()

# plot recruitment deviates
# note that the length of rec_dev is actually 34 not 35
gather_draws(stan_model_fit, rec_dev[ny]) %>% 
  ggplot(aes(x=ny, y=.value)) + 
  stat_lineribbon() + 
  scale_fill_brewer()

# plot raw
raw_v_time <- gather_draws(stan_model_fit, raw[ny]) %>%
  ggplot(aes(x=ny, y=.value)) + 
  stat_lineribbon() + 
  scale_fill_brewer()
ggsave(raw_v_time, filename=here("results","raw_v_time.png"))

# plot actual recruitment
n_p_a_y_hat <- gather_draws(stan_model_fit, n_p_a_y_hat[np, n_ages, ny])

n_p_a_y_hat %>%
  filter(n_ages==1) %>%
  ggplot(aes(x=ny, y=.value)) +
  stat_lineribbon() +
  scale_fill_brewer() +
  facet_wrap(~np)

# plot all ages over time
n_p_a_y_hat %>% 
  ggplot(aes(x=ny, y=.value)) +
  stat_lineribbon() +
  scale_fill_brewer() +
  facet_grid(np~n_ages)

proj_n_p_a_y_hat %>% 
  ggplot(aes(x=`(ny_proj + 1)`, y=.value)) +
  stat_lineribbon() +
  scale_fill_brewer() +
  facet_grid(np~n_ages)

# detection stats 
spread_draws(stan_model_fit, theta[patch,year]) %>% 
  ggplot(aes(x=year, y=theta)) + 
  stat_lineribbon() + 
  facet_wrap(~patch) +
  scale_fill_brewer()

beta_obs <- extract(stan_model_fit, "beta_obs")$beta_obs
hist(beta_obs)

# save model run if desired
saveRDS(stan_model_fit, here("results","summer_flounder_stan_fit.rds"))

###########
# calculate summary statistics and evaluate range shifts
###########

# centroid position by year 
dat_centroid <- abund_p_y %>% 
  group_by(year) %>% 
  summarise(centroid_lat = weighted.mean(x=patch, w=abundance))

# model fit centroid -- should eventually estimate in model for proper SE -- just exploring here
est_centroid <- abund_p_y_hat %>% 
  group_by(year, .iteration) %>%  # IS THIS SUPPOSED TO BE .ITERATION? CHECK WHEN MODEL IS RUN FOR LONGER 
  summarise(centroid_lat = weighted.mean(x=patch, w=dens_p_y_hat)) %>% 
  ungroup()

gg_centroid <- est_centroid %>% 
  ggplot(aes(year, centroid_lat)) + 
  stat_lineribbon() + 
  scale_fill_brewer() +
  geom_point(data = dat_centroid, aes(year, centroid_lat), color = "red") +
  theme(legend.position = "none") +
  labs(x="Year", y="Patch", title="Centroid Position", fill="Credible Interval")
gg_centroid
ggsave(gg_centroid, filename=here("results","centroid_v_time_no_length_comps.png"))

# centroid didn't shift at all!

# patch abundance fraction every year 
est_patch_abund <- abund_p_y_hat %>% 
  group_by(year, patch) %>% 
  summarise(abundance = mean(dens_p_y_hat))

observed_abundance_tile <- abund_p_y %>% 
  ggplot(aes(x=year, y=patch, fill=abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0, 36, 4)) +
  scale_y_continuous(breaks=seq(1, np, 1)) +
  labs(title="Observed", x="Year", y="Patch", fill="Abundance")


estimated_abundance_tile <- est_patch_abund %>% 
  ggplot(aes(x=year, y=patch, fill=abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0, 36, 4)) +
  scale_y_continuous(breaks=seq(1, np, 1)) +
  labs(title="Estimated")

ggsave(observed_abundance_tile, filename=here("results","observed_abundance_v_time_tileplot_no_length_comps.png"))
ggsave(estimated_abundance_tile, filename=here("results","estimated_abundance_v_time_tileplot_no_length_comps.png"))

# who's doing the colonizing?
dat_train_lengths %>% 
  group_by(patch, length) %>% 
  arrange(year) %>% 
  mutate(logratio = log(sum_num_at_length / lag(sum_num_at_length))) %>% 
  filter(logratio < Inf, logratio > -Inf) %>% 
  ungroup() %>% 
  ggplot(aes(x=year, y=logratio, group=length, color=length)) +
  geom_point() + 
  geom_line() +
  scale_color_viridis_c() +
  facet_wrap(~patch)

dat_train_lengths %>% 
  ggplot(aes(x=year, y=sum_num_at_length, fill=length)) + 
  geom_bar(stat="identity") +
  facet_wrap(~patch) +
  scale_fill_viridis_c() +
  theme_bw()


dat_train_lengths %>% 
  ggplot(aes(x=year, y=sum_num_at_length, fill=length)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_viridis_c() +
  theme_bw()  

# plot temperature difference from optimum over space and time
Topt <- extract(stan_model_fit, "Topt")$Topt
hist(Topt) # note, not normally distributed

dat_train_sbt %>%
  mutate(Tdiff = sbt - median(Topt)) %>%
  ggplot(aes(x=year, y=patch, fill=Tdiff)) +
  geom_tile() + 
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)

########
# evaluate forecast
########
proj_dens_p_y_hat <- gather_draws(stan_model_fit, proj_dens_p_y_hat[np, ny_proj])

proj_abund_p_y <- dat_test_dens %>%
  mutate(abundance = mean_dens * meanpatcharea)
# left_join(patchdat, by = c("lat_floor", "patch")) %>% 
# group_by(patch, year) %>% 
# summarise(abundance = sum(mean_dens *patch_area_km2)) %>% 
#  ungroup()

proj_observed_abundance_tile <- proj_abund_p_y %>% 
  mutate(Year = (year + min(years_proj) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>% 
  ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="Observed")


proj_est_patch_abund <- proj_dens_p_y_hat %>% 
  group_by(ny_proj, np) %>% 
  summarise(abundance = mean(.value))

proj_estimated_abundance_tile <- proj_est_patch_abund %>% 
  mutate(Year = (ny_proj + min(years_proj) - 1), Latitude = (np + min(patches) - 1), Abundance=abundance) %>% 
  ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="Estimated")
ggsave(proj_estimated_abundance_tile, filename=here("results","proj_estimated_abundance_v_time_tileplot.png"), scale=0.9)
ggsave(proj_observed_abundance_tile, filename=here("results","proj_observed_abundance_v_time_tileplot.png"), scale=0.9)

proj_abundance_v_time <- proj_dens_p_y_hat %>% 
  rename(patch=np) %>% 
  ggplot(aes(ny_proj, .value)) + 
  stat_lineribbon() + 
  geom_point(data = proj_abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch, scales = "free_y") +
  labs(x="Year",y="Abundance") + 
  scale_x_continuous(breaks=seq(0, 10, 2), limits=c(0, 10)) +
  theme(legend.position="none") +
  scale_fill_brewer()
ggsave(proj_abundance_v_time, filename=here("results","proj_density_v_time.png"), width=7, height=4)

# centroid position by year 
dat_centroid_proj <- proj_abund_p_y %>% 
  group_by(year) %>% 
  summarise(centroid_lat = weighted.mean(x=patch, w=abundance)) %>%
  mutate(Year = (year + min(years_proj) - 1), Latitude = (centroid_lat + min(patches) - 1))

# model fit centroid -- should eventually estimate in model for proper SE -- just exploring here
est_centroid_proj <- proj_dens_p_y_hat %>% 
  group_by(ny_proj, .iteration) %>%  # IS THIS SUPPOSED TO BE .ITERATION? CHECK WHEN MODEL IS RUN FOR LONGER 
  summarise(centroid_lat = weighted.mean(x=np, w=.value)) %>% 
  ungroup()

gg_centroid_proj_prep <- dat_centroid_proj %>% 
  ggplot(aes(Year, Latitude)) + 
  geom_point(color = "red") +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(37.6, 39.2, 0.2), limits=c(37.6, 39.2)) +
  theme(legend.position = "none") +
  labs(title="Centroid Position") 

gg_centroid_proj <- est_centroid_proj %>% 
  mutate(Year = (ny_proj + min(years_proj) - 1), Latitude = (centroid_lat + min(patches) - 1)) %>% 
  ggplot(aes(Year, Latitude)) + 
  stat_lineribbon() + 
  scale_fill_brewer() +
  geom_point(data = dat_centroid_proj, aes(Year, Latitude), color = "red") +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(37.6, 39.2, 0.2)) +
  theme(legend.position = "none") +
  labs(title="Centroid Position") 
ggsave(gg_centroid_proj_prep, filename=here("results","proj_centroid_v_time_prep.png"), scale=0.9)
ggsave(gg_centroid_proj, filename=here("results","proj_centroid_v_time.png"), scale=0.9)
