library(tidyverse)

library(ggdist)

library(rstanarm)

library(sdmTMB)

library(here)

library(reshape2)

library(ggridges)

normal <- rstanarm::normal

run_name <- "test"

results_path <- file.path("results", run_name)


load(here("processed-data", "stan_data_prep.Rdata"))

stan_model_fit <- read_rds(file.path(results_path,
                                     "stan_model_fit.rds"))


# prepare data ------------------------------------------------------------

density<- dens %>%
  as_tibble() %>%
  mutate(patch = factor(1:nrow(.))) %>%
  pivot_longer(
    -patch,
    names_to = "year",
    values_to = "density",
    names_transform = list(year = as.integer),
    names_prefix = "V"
  )

temperature <- sbt %>%
  as_tibble() %>%
  mutate(patch = factor(1:nrow(.))) %>%
  pivot_longer(
    -patch,
    names_to = "year",
    values_to = "sbt",
    names_transform = list(year = as.integer),
    names_prefix = "V"
  )



train <- temperature %>% 
  filter(year < 30)

moments <- train %>% 
  summarise(mean_sbt = mean(sbt),
            sd_sbt = sd(sbt))

temperature <- temperature %>% 
  mutate(cs_sbt = (sbt - moments$mean_sbt) / moments$sd_sbt) %>% 
  group_by(patch) %>% 
  mutate(l1_cs_sbt = lag(cs_sbt,1),
         l2_cs_sbt = lag(cs_sbt,2),
         l3_cs_sbt = lag(cs_sbt,3),
         l4_cs_sbt = lag(cs_sbt,4))



temperature_gradients <- expand_grid(p1 = 1:np, p2 = 1:np) %>% 
  mutate(distance = abs(p2 - p1)) %>% 
  filter(distance == 1) %>% 
  right_join(temperature %>% select(year, patch, cs_sbt) %>% mutate(patch = as.integer(patch)), by = c("p1" = "patch")) %>% 
  right_join(temperature %>% select(year, patch, cs_sbt) %>% mutate(patch = as.integer(patch)), by = c("p2" = "patch", "year")) %>% 
  mutate(gradient = cs_sbt.y - cs_sbt.x) %>% 
  select(p1, p2, year, gradient) %>% 
  group_by(p1, year) %>% 
  mutate(p2 = 1:n_distinct(p2)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "p2", values_from = "gradient",values_fill = 0,
              names_prefix = "neighbor") %>% 
  mutate(patch = factor(p1),
         year = year + 2)

temperature %>% 
  ggplot(aes(year, cs_sbt, color = patch)) + 
  geom_line()

reg_data <- density %>% 
  left_join(temperature, by = c("year", "patch")) %>% 
  left_join(temperature_gradients, by = c("year", "patch")) %>% 
  na.omit() %>% 
  group_by(patch) %>% 
  mutate(has_any = any(density > 0)) %>% 
  ungroup() %>% 
  filter(has_any) %>% 
  mutate(density = density + 1e-3) %>% 
  mutate(training = year < 30)


abund_p_y <- dat_train_dens %>%
  mutate(abundance = mean_dens * meanpatcharea) 

# fit GAM -----------------------------------------------------------------


density_gam <-
  stan_gamm4(density ~ patch + s(year) + s(cs_sbt) + s(l1_cs_sbt) + s(l2_cs_sbt) + s(neighbor1) + s(neighbor2),
             data = reg_data %>% filter(training), family = Gamma(link = "log"), 
             chains = 4, cores = 4, 
             iter = 5000)



density_ranger <-
  ranger(density ~ year + patch + cs_sbt + l1_cs_sbt + l2_cs_sbt + neighbor1 + neighbor2,
             data = reg_data %>% filter(training))

reg_data$ranger_pred = predict(density_ranger, data = reg_data)$predictions


predictions <- tidybayes::add_epred_draws(reg_data,density_gam, ndraws = 1000)


abund_p_y_hat <- tidybayes::spread_draws(stan_model_fit, dens_p_y_hat[patch,year])


gam_fit_plot <- predictions %>% 
  mutate(patch = as.numeric(as.character(patch))) %>% 
  ggplot() + 
  stat_lineribbon(aes(year, .epred, color = "GAM")) +
  geom_vline(aes(xintercept = 30), linetype = 2) +
  geom_line(aes(year, ranger_pred, color = "random forest"), size = 1.5) +
  # stat_lineribbon(data = abund_p_y_hat, aes(year, dens_p_y_hat, color = "fancy model")) +
  geom_point(data = abund_p_y %>% filter(patch %in% unique(predictions$patch)), aes(year, abundance), color = "orange") +
  facet_wrap(~patch, scales = "free_y") +
  labs(x="Year",y="Abundance") + 
  scale_fill_brewer() + 
  scale_color_manual(values = c("red","black"))

gam_fit_plot



# explore recruits --------------------------------------------------------

lcomps <- reshape2::melt(len) %>% 
  rename(patch = Var1, length = Var2, year = Var3)

lcomps %>% 
  ggplot(aes(length, value)) + 
  geom_col() + 
  facet_wrap(~patch)

a = lcomps %>% 
  ggplot(aes(x = length, y = forcats::fct_rev(as.factor(year)), height = value)) + 
  geom_density_ridges(stat = "identity", scale = 4) + 
  facet_wrap(~patch)
  
lcomps <- lcomps %>% 
  mutate(recruits = between(length, 25,40))

recruits <- lcomps %>% 
  filter(recruits) %>% 
  group_by(year, patch) %>% 
  summarise(recruits = sum(value)) %>% 
  ungroup() %>% 
  group_by(patch) %>% 
  mutate(recruits = recruits / max(recruits)) %>% 
  mutate(patch = as.factor(patch)) %>% 
  left_join(temperature, by = c("patch", "year"))

recruits %>%
  pivot_longer(contains("cs_sbt"), names_to = "temp_thing") %>%
  ggplot(aes(abs(value), recruits, color = patch)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ temp_thing) 

# fit sdmTMB --------------------------------------------------------------


