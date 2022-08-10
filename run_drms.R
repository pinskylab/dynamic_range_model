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

# ctrl_file <- ctrl_file %>%
#   pivot_longer(-c(id, description)) %>%
#   group_by(id,description) %>%
#   nest()

ctrl_file <-  ctrl_file %>%
  slice(1) %>% 
  mutate(fits = pmap(
    list(
      run_name = id,
      do_dirichlet = do_dirichlet,
      eval_l_comps = 1,
      T_dep_movement = T_dep_movement,
      spawner_recruit_relationship = spawner_recruit_relationship,
      process_error_toggle = process_error_toggle,
      exp_yn = exp_yn
    ),
    fit_drm,
    warmup = 2,
    iter = 4
  ))