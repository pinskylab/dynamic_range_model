fit_drm <- function(run_name = "test",
                    do_dirichlet = 1,
                    eval_l_comps = 0,
                    T_dep_mortality = 0,
                    T_dep_recruitment = 0,
                    T_dep_movement = 0,
                    spawner_recruit_relationship = 1,
                    run_forecast = 0,
                    process_error_toggle = 1,
                    exp_yn = 0,
                    warmup = 1000,
                    iter = 2000,
                    max_treedepth =  10,
                    chains =  1,
                    refresh = 10,
                    cores = 1,
                    drm_name = "process_sdm",
                    number_quantiles = 3,
                    quantiles_calc = c(0.05, 0.5, 0.95)) {
  
  
  results_path <- file.path("results",run_name)
  
  
  if (!dir.exists(results_path)){
    dir.create(results_path, recursive = TRUE)
  }
  
  
  load(here("processed-data","stan_data_prep.Rdata"))
  
  
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
    process_error_toggle = process_error_toggle,
    number_quantiles = number_quantiles,
    quantiles_calc = quantiles_calc
  )
  nums <- 100 * exp(-.2 * (0:(n_ages - 1)))
  
  stan_model_fit <- stan(file = here::here("src",paste0(drm_name,".stan")), # check that it's the right model!
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
  
  readr::write_rds(stan_model_fit, file = file.path(results_path,
                                                    "stan_model_fit.rds"))
  
  abund_p_y <- dat_train_dens %>%
    mutate(abundance = mean_dens * meanpatcharea) 
  
  abund_p_y_hat <- tidybayes::spread_draws(stan_model_fit, dens_p_y_hat[patch,year])
  
  
  abundance_v_time <- abund_p_y_hat %>% 
    ggplot(aes(year, dens_p_y_hat)) + 
    stat_lineribbon() +
    geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
    facet_wrap(~patch, scales = "free_y") +
    labs(x="Year",y="Abundance") + 
    scale_fill_brewer()
  # abundance_v_time
  ggsave(abundance_v_time, filename=file.path(results_path,"abundance_fits.pdf"), width=7, height=4)
  
  
  return(stan_model_fit)
  
} # close fit_drm