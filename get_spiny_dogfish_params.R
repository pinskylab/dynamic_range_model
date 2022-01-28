# this script extracts F values from PDFs from the summer flounder stock assessment 
# while tabulizer is great, please actually check the outputs against the PDFs!
library(tidyverse)
library(tabulizer)
library(here)


# from 2018 Status Report p23: Table 9a. Summary of stochastic fishing mortality rates expressed as the mean of full F on the exploitable biomass of female and male spiny dogfish, 1990-2017. Estimates for 2013 are not available. Year represents the year of the catch (landings plus dead discards). Estimates for 2016 are based on survey biomass from 2015 - 2017. Estimates for 2017 are based on biomass estimates from 2016- 2018. Sampling distribution of F estimates for females are given in Figure 11a,b. Fthreshold for females is 0.2439.
F_path <- here("spiny_dogfish_F_table.pdf")
F_out <- extract_areas(F_path, output="data.frame") # extract table -- have to do this one manually for some reason
F_df <- F_out[[1]] %>% # make into data frame
  slice(4:31) %>%  # get rid of the first row which was null 
  mutate(F1 = as.numeric(F1..Female),
         F2 = as.numeric(F2..Male)) %>% 
  rename(., Year = X)
F_df[F_df$Year==2013,]$F1 <- mean(c(F_df[F_df$Year==2012,]$F1,F_df[F_df$Year==2014,]$F1)) # fill in the missing year by sex
F_df[F_df$Year==2013,]$F2 <- mean(c(F_df[F_df$Year==2012,]$F2,F_df[F_df$Year==2014,]$F2))

F_df <- F_df %>% 
  rowwise() %>% 
  mutate(`F` = mean(c(F1, F2))) %>% 
  select(Year, `F`)

write_csv(F_df, here("processed-data","spiny_dogfish_F.csv"))
