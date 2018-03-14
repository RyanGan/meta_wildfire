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

# read in time series and limit to colorado and washington ----
ts <- read_csv("./data/health/1015-morbidity_pm_ts.csv") %>% 
  filter(!is.na(pm_krig)) %>% 
  rename(pm_diff = pm_smk) %>% 
  mutate(day = as.factor(weekdays(date)),
         weekend = ifelse(day %in% c("Saturday", "Sunday"), 1, 0),
         month = as.factor(lubridate::month(date)),
         year = as.factor(lubridate::year(date)),
         season = as.factor(case_when(month %in% c(12, 1, 2) ~ "winter",
                                      month %in% c(3:5) ~ "spring",
                                      month %in% c(6:8) ~ "summer",
                                      month %in% c(9:11)~ "fall")),
         smoke10_hms = if_else(pm_diff > 10 & month %in% c(3:11) & hms > 0.1,
                               1, 0),
         smoke15_hms = if_else(pm_diff > 15 & month %in% c(3:11) & hms > 0.1,
                               1, 0))

# read county pm timeseries
county_pm_ts <- read_csv("./data/smoke/1015-county_popwt_pm.csv") %>% 
  filter(str_sub(fips,start=1,end=2) %in% c("08", "53")) %>% 
  rename(pm_diff = pm_smk) %>% 
  mutate(smoke10_hms = if_else(pm_diff > 10 & month %in% c(3:11) & hms > 0.1,
                        1, 0),
         smoke15_hms = if_else(pm_diff > 15 & month %in% c(3:11) & hms > 0.1,
                      1, 0))

# defining a lag function
lags <- function(var, n=10){
  var <- enquo(var)
  indices <- seq_len(n)
  map( indices, ~quo(lag(!!var, !!.x)) ) %>% 
    set_names(sprintf("%s_lag%d", rlang::quo_text(var), indices))
  
}

# output krig pm lag matrix
pm_lag <- county_pm_ts %>% 
  filter(str_sub(fips,start=1,end=2) %in% c("08", "53")) %>% 
  select(fips, date, pm_krig, smoke10_hms, smoke15_hms) %>% 
  arrange(fips, date) %>% 
  group_by(fips) %>% 
  mutate(., !!!lags(pm_krig,6), !!!lags(smoke10_hms,6), !!!lags(smoke15_hms,6)) 

# join with morbidity data
ts_lag <- ts %>%
  select(-c(pm_krig, smoke10_hms, smoke15_hms)) %>% 
  left_join(pm_lag, by = c("fips", "date")) %>% 
  filter(!is.na(pm_krig_lag6)) %>% 
  arrange(fips, date) %>% 
  mutate(hms_percent = hms *100)

# output pm lag matrix
pm_mat <- as.matrix(select(ts_lag, contains("pm_krig")))/10

# creating a matrix of outcomes and binary smoke combos to run
# vector of names of outcomes to regress
outcomes <- c("resp", "asthma", "copd", "acute_bronch", "pneum", "cvd", 
              "arrhythmia", "cereb_vas", "hf", "ihd", "mi")

exposure <- c("smoke10_hms", "smoke15_hms")
# expand grid so each outcome/exposure combo is modeled 
exp_out_combo <- expand.grid(exposure, outcomes) %>% arrange(Var1)

# define time spline to adjust for
time_spl <- ns(ts_lag$date, df=12)

# parallel distributed lag computing ----
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
clusterExport(cl, c("ts_lag", "outcomes", "exposure", "exp_out_combo", "pm_mat",
                    "time_spl"), 
              envir = .GlobalEnv)

# Distributed lag association for smoke ----------------------------------------
# start time
start <- Sys.time()

# distributed lag function
smoke_dl_results <- parApply(cl, exp_out_combo, 1, function(x){
  # define outcome and exposure
  outcome <- x[2]
  exposure <- x[1]

  # output matrix of binary smoke
  smk_mat <- as.matrix(select(ts_lag, contains(exposure)))

  # define basis b using natural spline function
  exp_b <- ns(0:(ncol(smk_mat)-1), df = 3, intercept = T)
  # multiply lagged pm matrix by basis
  pm_basis <- pm_mat %*% exp_b
  smoke_basis <- smk_mat %*% exp_b
  
  # set up pm_smk and pm_nosmk basis
  pm_smk_basis <- pm_basis * smoke_basis
  pm_nosmk_basis <- pm_basis * ifelse(smoke_basis == 1, 0, 1)
  
  # fit mixed model
  mod <- glmer(as.formula(paste0(outcome,"~pm_smk_basis + pm_nosmk_basis +",
      "smoke_basis + as.factor(weekend) + temp_f + time_spl +",
      "offset(log(population)) + (1|fips)")),
               family = "poisson"(link="log"), data = ts_lag,
               control = glmerControl(optimizer = "bobyqa"))
  # test mod
  # mod <- glm(as.formula(paste0(outcome,"~pm_smk_basis + pm_nosmk_basis +",
  #    "smoke_basis + offset(log(population))")),
  #     family = "poisson"(link="log"), data = ts_lag)
  # summary(mod)
  
  # calculate estimates ----
  # output pm basis parameters
  # calculate estimates ----
  # output pm basis parameters
  pm_smk_parms <- broom::tidy(mod) %>% 
    filter(stringr::str_detect(term, "pm_smk")) %>% 
    select(estimate) %>% 
    as_vector()
  # extract names
  pm_smk_names <- broom::tidy(mod) %>% 
    filter(stringr::str_detect(term, "pm_smk")) %>% 
    select(term) %>% 
    as_vector()
  # no smoke basis params
  pm_nosmk_parms <- broom::tidy(mod) %>% 
    filter(stringr::str_detect(term, "pm_nosmk")) %>% 
    select(estimate) %>% 
    as_vector()
  # no smoke names
  pm_nosmk_names <- broom::tidy(mod) %>% 
    filter(stringr::str_detect(term, "pm_nosmk")) %>% 
    select(term) %>% 
    as_vector()
  #binary smoke basis
  smk_parms <- broom::tidy(mod) %>% 
    filter(stringr::str_detect(term, "smoke_basis")) %>% 
    select(estimate) %>% 
    as_vector()
  # binary smoke basis
  smk_names <- broom::tidy(mod) %>% 
    filter(stringr::str_detect(term, "smoke_basis")) %>% 
    select(term) %>% 
    as_vector()

  # estimate distributed lag values for each day
  pm_smk_estimate <- exp_b %*% pm_smk_parms
  pm_nosmk_estimate <- exp_b %*% pm_nosmk_parms
  smk_estimate <- exp_b %*% smk_parms
  
  # time variable
  time <- ((rep(1:length(smk_estimate))-1))
  
  # covariance matrix for knots 
  # output smoke basis matrix and nosmoke basis matrix
  pm_smk_cov_mat <- as.matrix(vcov(mod))[pm_smk_names, pm_smk_names]
  pm_nosmk_cov_mat <- as.matrix(vcov(mod))[pm_nosmk_names, pm_nosmk_names]
  smk_cov_mat <- as.matrix(vcov(mod))[smk_names, smk_names]

  # estimate variance of spline
  pm_smk_var <- exp_b %*% pm_smk_cov_mat %*% t(exp_b)
  pm_nosmk_var <- exp_b %*% pm_nosmk_cov_mat %*% t(exp_b)
  smk_var <- exp_b %*% smk_cov_mat %*% t(exp_b)
  
  # estimate lag ----
  # estimate standard error for each lag day for smoke
  l_pm_smk_se <- sqrt(diag(pm_smk_var))
  # calculate lower and upper bound for smoke
  l_pm_smk_l95 <- pm_smk_estimate+(l_pm_smk_se*qnorm(1-0.975))
  l_pm_smk_u95 <- pm_smk_estimate+(l_pm_smk_se*qnorm(0.975))
  
  # estimate standard error for each lag dat for nosmoke
  l_pm_nosmk_se <- sqrt(diag(pm_nosmk_var))
  # calculate lower and upper bound for smoke
  l_pm_nosmk_l95 <- pm_nosmk_estimate+(l_pm_nosmk_se*qnorm(1-0.975))
  l_pm_nosmk_u95 <- pm_nosmk_estimate+(l_pm_nosmk_se*qnorm(0.975))
  # smoke
  l_smk_se <- sqrt(diag(smk_var))
  # calculate lower and upper bound for smoke
  l_smk_l95 <- smk_estimate+(l_smk_se*qnorm(1-0.975))
  l_smk_u95 <- smk_estimate+(l_smk_se*qnorm(0.975))
  
  # list lag type of estimate
  l_type <- as.character(rep("lag", times = length(smk_estimate)*3))

  # smoke strata
  smk_strata <- as.character(c(rep("pm_smk", times = length(smk_estimate)),
                              rep("pm_nosmk", times = length(smk_estimate)),
                              rep("smk", times = length(smk_estimate))))

  # l_estimate
  l_estimate <- data.frame(exposure, outcome, l_type, smk_strata, time, 
    exp(c(pm_smk_estimate, pm_nosmk_estimate, smk_estimate)), 
    exp(c(l_pm_smk_l95, l_pm_nosmk_l95, l_smk_l95)), 
    exp(c(l_pm_smk_u95, l_pm_nosmk_u95, l_smk_u95)), 
                               row.names = NULL) 
  # assign column names
  colnames(l_estimate) <- c("exposure", "outcome", "type", "smoke", "time", 
                            "estimate", "lower95","upper95") 
  
  # estimate cumulative ----
  c_type <- as.character(rep("cumulative", times = length(smk_estimate)*2))
  # sequential cumulative estimate for smoke
  c_pm_smk_est <- sapply(seq_along(pm_smk_estimate), function(x){
    sum(pm_smk_estimate[1:x])
  })
  # stderr cumulative effect for smoke
  c_pm_smk_se <- sapply(seq_along(pm_smk_estimate), function(y){
    sqrt(sum(pm_smk_var[1:y,1:y]))
  })
  
  # sequential cumulative estimate for nosmoke
  c_pm_nosmk_est <- sapply(seq_along(pm_nosmk_estimate), function(x){
    sum(pm_nosmk_estimate[1:x])
  })
  # stderr cumulative effect for smoke
  c_pm_nosmk_se <- sapply(seq_along(pm_nosmk_estimate), function(y){
    sqrt(sum(pm_nosmk_var[1:y,1:y]))
  })
  # sequential cumulative estimate for binary smoke
  c_smk_est <- sapply(seq_along(smk_estimate), function(x){
    sum(smk_estimate[1:x])
  })
  # stderr cumulative effect for smoke
  c_smk_se <- sapply(seq_along(smk_estimate), function(y){
    sqrt(sum(smk_var[1:y,1:y]))
  })
  
  # cumulative 96CI for pm smoke
  c_pm_smk_l95 <- c_pm_smk_est+(c_pm_smk_se*qnorm(1-0.975))
  c_pm_smk_u95 <- c_pm_smk_est+(c_pm_smk_se*qnorm(0.975))
  # cumulative 96CI for pm no smoke
  c_pm_nosmk_l95 <- c_pm_nosmk_est+(c_pm_nosmk_se*qnorm(1-0.975))
  c_pm_nosmk_u95 <- c_pm_nosmk_est+(c_pm_nosmk_se*qnorm(0.975))
  # cumulative 95CI for binary smoke
  c_smk_l95 <- c_smk_est+(c_smk_se*qnorm(1-0.975))
  c_smk_u95 <- c_smk_est+(c_smk_se*qnorm(0.975))
  
  # list lag type of estimate
  c_type <- as.character(rep("cumulative", times = length(smk_estimate)*3))

  # c_estimate
  c_estimate <- data.frame(exposure, outcome, c_type, smk_strata, time, 
    exp(c(c_pm_smk_est, c_pm_nosmk_est, c_smk_est)), 
    exp(c(c_pm_smk_l95, c_pm_nosmk_l95, c_smk_l95)), 
    exp(c(c_pm_smk_u95, c_pm_nosmk_u95, c_smk_u95)), 
    row.names = NULL) 
  # assign column names
  colnames(c_estimate) <- c("exposure", "outcome", "type", "smoke", "time", 
                            "estimate", "lower95","upper95") 
  # bind lag and cumulative estimates together
  return_estimate <- rbind(c_estimate, l_estimate)
  # print check
  print(return_estimate)
  # return estimate  
  return(return_estimate)
  }) %>%   # end parLapply
  # bind list together as rows
  map_dfr(., bind_rows)

# check
print(head(smoke_dl_results))

# print warnings
warnings()

# write file ----
write_csv(smoke_dl_results, paste0("./data/health/",
  "1015-ts_dl_int_tspl_adjsmk_results.csv"))

# stop time
stop <- Sys.time()
time <- stop - start
# print time
print(time)

