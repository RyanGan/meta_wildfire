# ------------------------------------------------------------------------------
# Title: Distributed Lag time-stratified case-crossover 2010 to 2015 analysis
# Author: Ryan Gan
# Date Created: 2018-03-16
# ------------------------------------------------------------------------------

# Purpose: Running time-series analysis between morbidity/mortality outcomes
# and PM and smoke estimates

# libraries ----
library(tidyverse)
library(survival) # loading lme4 mixed model package
library(splines)
library(stringr)
library(lubridate)

# load case-crossover dataframe list
load("./data/health/1015-morbidity_casecross_list.RData")

# list of outcome names
out_names <- names(casecross_list)

# filter list to just washington and colorado and months april through october
casecross_analysis_list <- casecross_list %>% 
  map(~filter(., state != "oregon") %>% 
        mutate(month = as.factor(lubridate::month(date))) %>% 
        filter(month %in% 4:10)) %>% 
  # add outcome name to each dataframe
  map2(.x = ., .y = out_names, ~mutate(.x, out_name = .y))

# remove casecross_list to save memory
rm(casecross_list)

# load environmental PM data and create lags each year -------------------------
# defining a lag function
funlag <- function(var, n=7){
  var <- enquo(var)
  indices <- seq_len(n)
  map( indices, ~quo(lag(!!var, !!.x)) ) %>% 
    set_names(sprintf("%s_lag%d", rlang::quo_text(var), indices))
}

# read pm dataframe and apply lag function 
pm <- read_csv("./data/smoke/1015-county_popwt_pm.csv") %>% 
  filter(str_sub(fips,start=1,end=2) %in% c("53", "08")) %>% 
  rename(pm_diff = pm_smk) %>% 
  mutate(day = as.factor(weekdays(date)),
         weekend = ifelse(day %in% c("Saturday", "Sunday"), 1, 0),
         month = as.factor(lubridate::month(date)),
         year = as.factor(lubridate::year(date)),
         season = as.factor(case_when(month %in% c(12, 1, 2) ~ "winter",
                                      month %in% c(3:5) ~ "spring",
                                      month %in% c(6:8) ~ "summer",
                                      month %in% c(9:11)~ "fall")),
         smoke0_hms = ifelse(pm_diff > 0 & month %in% c(3:10) & hms > 0.5, 1, 0),
         smoke5_hms = ifelse(pm_diff > 5 & month %in% c(3:10) & hms > 0.5, 1, 0),
         smoke10_hms = ifelse(pm_diff > 10 & month %in% c(3:10) & hms > 0.5, 1, 0),
         smoke15_hms = ifelse(pm_diff > 15 & month %in% c(3:10) & hms > 0.5, 1, 0),
         pm = pm_krig/10) %>% 
  arrange(fips, date) %>% 
  group_by(fips) %>% 
  mutate(., !!!funlag(pm,6), !!!funlag(smoke0_hms,6), !!!funlag(smoke5_hms,6),
         !!!funlag(smoke10_hms,6), !!!funlag(smoke15_hms,6)) %>% 
  select(-state, -month)

# function to extract dl estimates ---------------------------------------------
distribute_that_lag <- function(lag_mod, strata) {
  # output pm basis estimates
  parms <- broom::tidy(lag_mod) %>% 
    filter(stringr::str_detect(term, strata)) %>% 
    select(estimate) %>% 
    as_vector()
  # output estimate names for cov matrix
  names <- stringr::str_subset(names(lag_mod$coefficients), strata)
  # estimate associations
  est <- exp_b %*% parms
  # estimate standard error for each interval
  # time variable
  time <- ((rep(1:length(est))-1))
  # covariance matrix for knots 
  cov_mat <- as.matrix(vcov(lag_mod))[names, names]
  # estimate variance of spline
  var <- exp_b %*% cov_mat %*% t(exp_b)
  # estimate lag ----
  # estimate standard error for each lag day for smoke
  l_se <- sqrt(diag(var))
  # calculate lower and upper bound for smoke
  l_est_l95 <- est + (l_se*qnorm(1-0.975))
  l_est_u95 <- est + (l_se*qnorm(0.975))
  l_type <- "lag"
  # lag dataframe
  l_df <- data.frame(strata, l_type, time, 
                     exp(est), exp(l_est_l95), exp(l_est_u95), 
                     row.names = NULL) 
  # assign column names
  colnames(l_df) <- c("strata", "type", "time", 
                      "odds_ratio", "lower_95", "upper_95")
  # cumulative estimates
  c_est <- sapply(seq_along(est), function(x){
    sum(est[1:x])
  })
  # stderr cumulative effect smk
  c_se <- sapply(seq_along(c_est), function(y){
    sqrt(sum(var[1:y,1:y]))
  })
  # estimate 95% CI
  c_l95 <- c_est+(c_se*qnorm(1-0.975))
  c_u95 <- c_est+(c_se*qnorm(0.975))
  # type
  c_type <- "cumulative"
  # return dataframe
  c_df <- data.frame(strata, c_type, time, exp(c_est), 
                     exp(c_l95), exp(c_u95), row.names = NULL) 
  # assign column names
  colnames(c_df) <- c("strata", "type", "time", 
                      "odds_ratio", "lower_95", "upper_95")
  # bind lagged and cumulative 
  lag_est <- rbind(l_df, c_df) %>% 
    mutate(strata = as.character(strata),
           type = as.character(type))
  # return lagged estimate
  return(lag_est)
} # end lag estimate function

# continous PM2.5 distributed lag results --------------------------------------
# start time
start <- Sys.time()
# distributed lag function
casecross_dl_pm_results  <- lapply(casecross_analysis_list, function(x){
  # output dataframe from list
  data <- x %>% 
    mutate(date = as.Date(date),
           outcome = as.numeric(as.character(outcome))) %>%
    left_join(pm, by = c("date", "fips")) %>% 
    # remove missing lagged data
    filter(!is.na(pm_lag6))

  # output outcome name
  out_name <- as.character(unique(data$out_name))
  print(out_name)
  # create lagged matrix
  pm_mat <- as.matrix(select(data, pm, contains("pm_lag")))
  # define lagged basis spline
  exp_b <- ns(0:(ncol(pm_mat)-1), df = 3, intercept = T)
  # pm basis
  pm_basis <- pm_mat %*% exp_b
  # run lagged model
  lag_mod <- clogit(outcome ~ pm_basis + temp_f + strata(id), data = data)

  # estimate lag estimate
  lag_est <- distribute_that_lag(lag_mod, strata = "pm") %>% 
    mutate(outcome = out_name) %>% select(outcome, strata:upper_95)
  return(lag_est)
  }) %>%  #end lappply
  # bind rows
  map_dfr(.,rbind)
# stop time
stop <- Sys.time()
time <- stop - start
# print time
print(time)

# saving pm results
write_csv(casecross_dl_pm_results, "./data/health/1015-cc_dl_pm_results.csv")

# binary smoke 0 ---------------------------------------------------------------
start <- Sys.time()
# distributed lag function
dl_smk10_results  <- lapply(casecross_analysis_list, function(x){
  # output dataframe from list
  data <- x %>% 
    mutate(date = as.Date(date),
           outcome = as.numeric(as.character(outcome))) %>%
    left_join(pm, by = c("date", "fips")) %>% 
    # remove missing lagged data
    filter(!is.na(pm_lag6))
  
  # output outcome name
  out_name <- as.character(unique(data$out_name))
  print(out_name)
  # create lagged matrix
  smk_mat <- as.matrix(select(data, contains("smoke10_hms")))
  # define lagged basis spline
  exp_b <- ns(0:(ncol(smk_mat)-1), df = 3, intercept = T)
  # pm basis
  smk_basis <- smk_mat %*% exp_b
  # run lagged model
  lag_mod <- clogit(outcome ~ smk_basis + temp_f + strata(id), data = data)
  # estimate lag estimate
  lag_est <- distribute_that_lag(lag_mod, strata = "smk") %>% 
    mutate(outcome = out_name) %>% select(outcome, strata:upper_95)
  return(lag_est)
}) %>%  #end lappply
  # bind rows
  map_dfr(.,rbind)
# stop time
stop <- Sys.time()
time <- stop - start
# print time
print(time)

# saveing smoke10 results
write_csv(dl_smk10_results, "./data/health/1015-cc_dl_smk10_results.csv")

# Interaction model ------------------------------------------------------------

# start time
start <- Sys.time()
# distributed lag function
casecross_dl_int_results  <- lapply(casecross_analysis_list, function(x){
  # output dataframe from list
  data <- x %>% 
    mutate(date = as.Date(date),
      outcome = as.numeric(as.character(outcome))) %>%
      left_join(pm, by = c("date", "fips")) %>% 
    # remove missing lagged data
    filter(!is.na(pm_lag6))

  # output outcome name
  out_name <- unique(data$out_name)
  # create lagged matrix
  pm_mat <- as.matrix(select(data, pm, contains("pm_lag")))
  # define lagged basis spline
  exp_b <- ns(0:(ncol(pm_mat)-1), df = 3, intercept = T)
  # create vector of smoke var to estimate over
  smoke_var <- c("smoke0_hms", "smoke5_hms", "smoke10_hms", "smoke15_hms")
  
  smoke_results <- smoke_var %>% 
    map(function(s){ 
      smk_mat <- as.matrix(select(data, contains(s)))
      # lagged pm x basis
      pm_basis <- pm_mat %*% exp_b
      smk_basis <- smk_mat %*% exp_b
      # smoke basis for interaction
      pm_smk_b <- pm_basis * smk_basis 
      pm_nosmk_b <- pm_basis * ifelse(smk_basis == 1, 0, 1)
      # run lagged model
      lag_mod <- clogit(outcome ~ pm_smk_b + pm_nosmk_b + smk_basis + temp_f +
                        strata(id), data = data)
      # define strata terms
      strata_terms <- c("pm_smk_b", "pm_nosmk_b", "smk_basis")
      # estimate cumulative and lagged effect for each basis
      
      lagged_estimates <- strata_terms %>% 
        map_dfr(~distribute_that_lag(lag_mod = lag_mod, strata = .)) %>% 
        mutate(outcome = out_name, smoke = s) %>% 
        select(outcome, smoke, strata:upper_95)
                
      return(lagged_estimates)      
    }) %>%  # end smoke map
      # bind rows
      map_dfr(.,rbind)
  }) %>%  #end lappply
  # bind rows
  map_dfr(.,rbind)

# stop time
stop <- Sys.time()
time <- stop - start
# print time
print(time)

# write file
write_csv(casecross_dl_int_results, "./data/health/1015-cc_dl_int_results.csv")

        