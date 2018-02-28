# ------------------------------------------------------------------------------
# Title: Time-series 2010 to 2015 analysis
# Author: Ryan Gan
# Date Created: 2018-02-28
# ------------------------------------------------------------------------------

# Purpose: Running time-series analysis between morbidity/mortality outcomes
# and PM and smoke estimates

# libraries ----
library(tidyverse)
library(lme4) # loading lme4 mixed model package
library(parallel)

# read in time series ----
ts <- read_csv("./data/health/1015-morbidity_pm_ts.csv") %>% 
  filter(!is.na(pm_krig)) %>% 
  mutate(day = as.factor(weekdays(date)),
         weekend = ifelse(day %in% c("Saturday", "Sunday"), 1, 0),
         month = as.factor(lubridate::month(date)),
         year = as.factor(lubridate::year(date)),
         season = as.factor(case_when(month %in% c(12, 1, 2) ~ "winter",
                            month %in% c(3:5) ~ "spring",
                            month %in% c(6:8) ~ "summer",
                            month %in% c(9:11)~ "fall")))

# same-day association of binary smoke variables, accounting for county, poulation
# month, year, state, weekend
# vector of names of outcomes to regress

outcomes <- c("resp", "asthma", "copd", "acute_bronch", "pneum", "cvd", 
              "arrhythmia", "cereb_vas", "hf", "ihd", "mi", "cardiopulm_n")

exposure <- c("smoke0", "smoke5", "smoke10", "smoke15", "smoke_wave")

exp_out_combo <- expand.grid(exposure, outcomes) %>% arrange(Var1)

print(head(exp_out_combo))
# set up cluster of 8 cores to parallelize models
cores <- parallel::detectCores()
# check to see if cores detected
print(cores)
cl <- makeCluster(cores)
# check cluster
print(cl)

# load packages on each processor of the node/cluster
clusterCall(cl, function() c(library(tidyverse), library(lme4), 
                             library(splines)))

# export global objects to cluster
clusterExport(cl, c("ts", "exp_out_combo"), 
              envir = .GlobalEnv)

# binary result test results ----
# start time
start <- Sys.time()
# parallel apply
binary_smoke_results <- parApply(cl, exp_out_combo,1, function(x){
  # define outcome and exposure
  outcome <- x[2]
  exposure <- x[1]
  # run model
  mod <- glmer(as.formula(paste0(outcome, "~", 
    exposure, "+weekend+month+year+state+(1|fips)+offset(log(population))")),
      family = "poisson"(link="log"), data = ts,
      control = glmerControl(optimizer = "bobyqa"))

  # broom:tidy the model
  result <- broom::tidy(mod) %>% 
    filter(term == exposure) %>% 
    mutate(rr = round(exp(estimate),3),
           lower95 = round(exp(estimate-(1.96*std.error)),3),
           upper95 = round(exp(estimate+(1.96*std.error)),3),
           p = round(p.value,3)) %>% 
    select(rr:p)
  # result dataframe
  exp_out <- as_data_frame(cbind(exposure, outcome))
  output_result <- bind_cols(exp_out, result)
  return(output_result)
  }) %>% 
  # bind list together as rows
  map_dfr(., bind_rows)

# stop time
stop <- Sys.time()
time <- stop - start
# print time
print(time)
# check
print(head(binary_smoke_results))

# write file 
write_csv(binary_smoke_results, "./data/health/1015-ts_binary_smoke_results.csv")

# close cluster
stopCluster(cl)

