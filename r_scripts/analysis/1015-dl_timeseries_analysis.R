# ------------------------------------------------------------------------------
# Title: Distributed Lag time-series 2010 to 2015 analysis
# Author: Ryan Gan
# Date Created: 2018-02-28
# ------------------------------------------------------------------------------

# Purpose: Running time-series analysis between morbidity/mortality outcomes
# and PM and smoke estimates

# libraries ----
library(tidyverse)
library(lme4) # loading lme4 mixed model package
library(splines)
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


exposure <- c("smoke10", "smoke15", "smoke_wave")
# expand grid so each outcome/exposure combo is modeled 
exp_out_combo <- expand.grid(exposure, outcomes) %>% arrange(Var1)

# read county pm timeseries
county_pm_ts <- read_csv("./data/smoke/1015-county_popwt_pm.csv")

# binary smoke lagged dataframe
smoke10_lag <- county_pm_ts %>% 
  filter(str_sub(fips,start=1,end=2) %in% c("08", "53")) %>% 
  select(fips, date, smoke10) %>% 
  arrange(fips, date) %>% 
  group_by(fips) %>% 
  mutate(smk_lag1 = lag(smoke10, 1, order_by = fips),
         smk_lag2 = lag(smoke10, 2, order_by = fips),
         smk_lag3 = lag(smoke10, 3, order_by = fips),
         smk_lag4 = lag(smoke10, 4, order_by = fips),         
         smk_lag5 = lag(smoke10, 5, order_by = fips),
         smk_lag6 = lag(smoke10, 6, order_by = fips))


# join with morbidity data
ts_lag <- ts %>%
  select(-smoke10) %>% 
  left_join(smoke10_lag, by = c("fips", "date")) %>% 
  filter(complete.cases(.)) 

# parallel distributed lag computing ----
# output matrix of smoke and lag variables
smk_matrix <- as.matrix(select(ts_lag, smoke10:smk_lag6))

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

# export new set of objects to global objects to cluster
clusterExport(cl, c("ts_lag", "smk_matrix", "outcomes", "exposure"), 
              envir = .GlobalEnv)
# Distributed lag association for smoke10 --------------------------------------
# start time
start <- Sys.time()

# distributed lag function
smoke_dl_results <- parApply(cl, exp_out_combo,1, function(x){
  # define outcome and exposure
  outcome <- x[2]
  exposure <- x[1]
  # define basis b using natural spline function
  smk_b <- ns(0:(ncol(smk_matrix)-1), df = 3, intercept = T)
  # multiply lagged pm matrix by basis
  smk_basis <- smk_matrix %*% smk_b
  # fit mixed model
  mod <- glmer(as.formula(paste0(outcome,"~", exposure, 
    "smk_basis + as.factor(weekend) + state + as.factor(month) +",
     "as.factor(year) + (1|fips) + offset(log(population))")),
               family = "poisson"(link="log"), data = ts_lag,
               control = glmerControl(optimizer = "bobyqa"))

  # calculate estimates ----
  # output pm basis parameters
  dl_parms <- broom::tidy(mod) %>% 
    filter(stringr::str_detect(term, "smk_basis")) %>% 
    select(estimate) %>% 
    as_vector()
  
  # estimate distributed lag values for each day
  estimate <-  smk_b %*% dl_parms
  # time variable
  time <- ((rep(1:length(estimate))-1))
  # covariance matrix for knots 
  # fix the matrix
  cov_matrix <- as.matrix(vcov(mod))[2:(3+1), 2:(3+1)]
  # estimate variance of spline
  variance <- smk_b %*% cov_matrix %*% t(smk_b)
  # estimate lag ----
  # estimate standard error for each lag day
  l_se <- sqrt(diag(variance))
  # calculate lower and upper bound
  l_lower95 <- estimate+(l_se*qnorm(1-0.975))
  l_upper95 <- estimate+(l_se*qnorm(0.975))
  l_type <- as.character(rep("lag", times = length(estimate)))
  
  # l_estimate
  l_estimate <- data.frame(outcome, l_type, time, exp(estimate), 
                           exp(l_lower95), exp(l_upper95), row.names = NULL) 
  # assign column names
  colnames(l_estimate) <- c("outcome", "type","time", "estimate",
                            "lower95","upper95") 
  
  # estimate cumulative ----
  c_type <- as.character(rep("cumulative", times = length(estimate)))
  # sequential cumulative estimate
  cumulative_estimate <- sapply(seq_along(estimate), function(x){
    sum(estimate[1:x])
  })
  # stderr cumulative effect
  cumulative_se <- sapply(seq_along(estimate), function(y){
    sqrt(sum(variance[1:y,1:y]))
  })
  # cumulative 95CI
  c_lower95 <- cumulative_estimate+(cumulative_se*qnorm(1-0.975))
  c_upper95 <- cumulative_estimate+(cumulative_se*qnorm(0.975))
  # return dataframe
  c_estimate <- data.frame(outcome, c_type, time, exp(cumulative_estimate), 
                           exp(c_lower95), exp(c_upper95), row.names = NULL)
  # assign column names
  colnames(c_estimate) <- c("outcome", "type","time","estimate",
                            "lower95","upper95") 
  # bind lag and cumulative estimates together
  return_estimate <- rbind(c_estimate, l_estimate)
  # return estimate  
  return(return_estimate)
  }) %>%   # end parLapply
  # bind list together as rows
  map_dfr(., bind_rows)

# check
print(head(smoke_dl_results))

# write file ----
write_csv(smoke_dl_results, "./data/health/1015-ts_smoke_dl_results.csv")

# stop time
stop <- Sys.time()
time <- stop - start
# print time
print(time)

