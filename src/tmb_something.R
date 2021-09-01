library(LIME)
library(tidyverse)


fish <- fish

fleet <- pfo$prepped_fishery[[1]]$fleet

scrooge_lh <- create_lh_list(vbk= fish$vbk,
                             linf= fish$linf,
                             t0= fish$t0,
                             lwa= fish$weight_a,
                             lwb= fish$weight_b,
                             S50= fleet$length_50_sel,
                             S95=fleet$length_95_sel,
                             selex_input="length",
                             selex_type=c("logistic"),
                             M50= fish$length_50_mature,
                             M95= fish$length_95_mature,
                             maturity_input="length",
                             M= fish$m,
                             binwidth=1,
                             CVlen= fish$cv_len,
                             SigmaR= fish$sigma_r + .001,
                             SigmaF= fleet$sigma_effort + .001,
                             SigmaC=0.2,
                             SigmaI=0.2,
                             R0= fish$r0,
                             qcoef=1e-5,
                             start_ages=0,
                             rho=0.43,
                             nseasons=1)

temp_LF_matrix <- pfo$prepped_fishery[[1]]$length_comps


temp_LF_matrix %>%
  gather(length, numbers, -year) %>%
  mutate(length = as.numeric(length)) %>%
  ggplot(aes(length, numbers, fill = factor(year))) +
  geom_density(stat = "identity", alpha = 0.5)

LF_matrix <- temp_LF_matrix %>%
  as.matrix()

rownames(LF_matrix) <- LF_matrix[,"year"]

LF_matrix <- LF_matrix[, -c(1)]

scrooge_data_LF <-
  list("years" = 1:nrow(LF_matrix), "LF" = LF_matrix)

start <- Sys.time()
res <- run_LIME(modpath=NULL,
                lh=scrooge_lh,
                input_data=scrooge_data_LF,
                est_sigma="log_sigma_R",
                data_avail="LC",
                newtonsteps=3)
end <- Sys.time() - start


check <- res$df

## check for other issues
issues <- res$issues

## check TMB inputs
Inputs <- res$Inputs

## Report file
Report <- res$Report

## Standard error report
Sdreport <- res$Sdreport

##----------------------------------------------------
## Step 4: Plot fits
## ---------------------------------------------------
## plot length composition data
plot_LCfits(LFlist=list("LF"=LF_matrix),
            Inputs=Inputs,
            Report=Report)

true_f <- pfo$prepped_fishery[[1]]$simed_fishery %>%
  group_by(year) %>%
  summarise(f = unique(f))

true_f$predicted_f <- Report$F_y

predicted_f <- data_frame(f = Report$F_y)

true_f %>%
  ggplot() +
  geom_point(aes(year,f,color = "True")) +
  geom_line(aes(year, predicted_f,color = "Predicted"))

