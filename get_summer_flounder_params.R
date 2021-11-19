# this script extracts F values from PDFs from the summer flounder stock assessment 
# while tabulizer is great, please actually check the outputs against the PDFs!
library(tidyverse)
library(tabulizer)
library(here)

# from SAW66 p239: "Table A88. 2018 SAW­66 assessment fishing mortality (F) estimates at age; F2018_BASE_V2 model run."
F_age_path <- here("summer_flounder_F_table_by_age.pdf")
F_age_out <- extract_tables(F_age_path, output="data.frame") # extract table
F_age_df <- F_age_out[[1]] %>% # make into data frame
  set_names(c('Year',seq(0, 7, 1))) %>% # add column names back in -- get from the pdf
  pivot_longer(cols=2:9, names_to="Age", values_to="F") # tidy df 

# from SAW66 p238: Table A87. 2018 SAW­66 assessment summary results for Spawning Stock Biomass (SSB) in metric tons (mt); Recruitment (R) at age 0 (000s); Fishing Mortality (F) for fully recruited (peak) age 4; F2018_BASE_V2 model run.
F_path <- here("summer_flounder_F_table.pdf")
F_out <- extract_areas(F_path, output="data.frame") # extract table -- have to do this one manually for some reason
F_df <- F_out[[1]] %>% # make into data frame
  slice(2:37) %>%  # get rid of the first row which was null 
  select(Year, `F`) # also has SSB and R columns that we don't need at present

# from SAW66 p173-174 Table A33. Mean weight (kg) at age of summer flounder catch, Maine­North Carolina. Includes ‘New’ MRFSS/MRIP.
wt_path <- here("summer_flounder_weight_at_age_tables.pdf")
wt_out <- extract_tables(wt_path, output="data.frame")
wt_df1 <- wt_out[[1]]
wt_df2 <- wt_out[[2]]
wt_df <- rbind(wt_df1, wt_df2) %>% 
  select(1:12, 14) %>% # drop total column
  set_names(c('Year',seq(0, 10, 1), "over7")) %>% # add column names back in -- get from the pdf
  pivot_longer(cols=2:13, names_to="Age", values_to="Wt") # tidy df 

write_csv(F_age_df, here("processed-data","summer_flounder_F_by_age.csv"))
write_csv(F_df, here("processed-data","summer_flounder_F.csv"))
write_csv(wt_df, here("processed-data","summer_flounder_wt_at_age.csv"))
