# ------------------------------------------------------------------------------
# Title: Distributed lag random-effects models for Meta-Wildfire project
# Author: Ryan Gan
# Date Created: 2018-02-07
# ------------------------------------------------------------------------------

# Purpose: This script performs the distributed lag analyses on county-specific
# time series data. The random-effects models take a while to run.
# I've split this code out to run on the server; consider parallel for efficiency.
# Results are saved as .csv files.

# load libraries ----
library(tidyverse)
library(lme4)
library(splines)


# import population denominators ----
# import pm2.5 data -----
# read pm data and subset to colorado
pm_path <- paste0("./data/smoke/2015-smoke_wus_county_popwt.csv")
# read pm
pm <- read_csv(pm_path) %>% 
  mutate(smoke0 = ifelse(gwr_smk > 0, 1, 0),
         smoke5 = ifelse(gwr_smk > 5, 1, 0),
         smoke10 = ifelse(gwr_smk > 10, 1, 0),
         smoke15 = ifelse(gwr_smk > 15, 1, 0))

# colorado census
co_denom <- read_csv("./data/census/2016-colorado_population.csv") %>% 
  select(GEO.id2, 'GEO.display-label', respop72010:respop72016) %>% 
  rename(fips = GEO.id2, county = 'GEO.display-label') %>% 
  rename_at(vars(respop72010:respop72016), 
            funs(paste0("pop",str_sub(., start=8)))) %>% 
  # take out state name and "county" from county name
  mutate(county = str_split_fixed(county, " County", n=2)[,1])

# washington cesus
wa_denom <- read_csv("./data/census/2016-washington_population.csv") %>% 
  select(GEO.id2, 'GEO.display-label', respop72010:respop72016) %>% 
  rename(fips = GEO.id2, county = 'GEO.display-label') %>% 
  rename_at(vars(respop72010:respop72016), 
            funs(paste0("pop",str_sub(., start=8)))) %>% 
  # take out state name and "county" from county name
  mutate(county = str_split_fixed(county, " County", n=2)[,1])

county_denom <- wa_denom %>% 
  mutate(fips = as.character(fips)) %>% 
  bind_rows(co_denom) %>% 
  # add state variable in
  mutate(state = if_else(str_sub(fips, 1, 2)=="08",
                         "Colorado", "Washington")) %>% 
  select(fips, state, county, pop2015) %>% 
  rename(population = pop2015)

# import outcomes data ----
# define relative path to data folder
co_path <- paste0("./data/health/2010-2015_morbidity_co_ts.csv")
# read in dataframe using read_csv function from tidyverse package
co_ts <- read_csv(co_path, progress = F) 
# define relative path
wa_path <- paste0("./data/health/2011-2015_morbidity_wa_ts.csv")
# read in dataframe
wa_ts <- read_csv(wa_path, progress = F) %>% 
  mutate(fips = as.character(fips))

# bind time series together
ed_ts <- co_ts %>% 
  # bind colorado and washington time series together
  bind_rows(wa_ts) %>% 
  # filter to dates I have PM2.5 estimates for
  filter(date >= "2015-06-01" & date <= "2015-09-30") %>% 
  mutate(year = lubridate::year(date)) %>% 
  # join with county denoms
  left_join(county_denom, by = "fips") %>% 
  # join with pm values
  left_join(pm, by = c("date", "fips")) %>% 
  # create weekend variable
  mutate(day = as.factor(weekdays(date)),
         weekend = ifelse(day %in% c("Saturday", "Sunday"), 1, 0),
         month = as.factor(lubridate::month(date)),
         # create pm 10 unit change v ariables
         gwr_smk10 = gwr_smk/10,
         krig_smk10 = krig_smk/10,
         gwr10 = gwr/10,
         krig10 = krig/10)

# subset to total ed strata
total_ed <- ed_ts %>% 
  filter(strata == "total") %>% 
  filter(!is.na(gwr_smk10))

# background PM spline ---------
# define degrees of freedom
spl_vec <- 5:20

# fit spline by date
pm_spline_fits <- spl_vec %>% 
  map_dfr(function(df_val){
    date <- total_ed$date
    # run spline based on date sequence
    pm_spl <- ns(date, df = df_val)
    # model gwr_smoke units regressed on the spline
    pm_mod <- lm(gwr_smk10 ~ pm_spl, data = total_ed)
    # estimated predicted values based on the spline and bind with dates
    pm_pred <- predict(pm_mod, type="response")
    # rep df val by length of prediction
    degree_freedom <- rep(df_val, length(pm_pred))
    # rep aic val by lenth of prediction 
    aic <- round(rep(AIC(pm_mod), length(pm_pred)),2)
    # join with total_ed 
    plot_df <- total_ed %>% 
      select(date, fips, gwr_smk10) %>% 
      cbind(pm_pred, degree_freedom, aic)
  } # end function
) # end map 

# save spline fits
write_csv(pm_spline_fits, "./data/health/2015-ts_pm_spline_fit_results.csv")

# distributed lag ----
# define outcomes to loop through
outcomes <- total_ed %>% 
  select(resp:mi) %>% 
  colnames()
         
# general code to create lagged variables
# may be worth creating in a function one day
gwr_lag <- pm %>% 
  select(fips, date, gwr_smk) %>% 
  arrange(fips, date) %>% 
  group_by(fips) %>% 
  mutate(gwr_smk10 = gwr_smk/10,
         gwr_smk_lag1 = lag(gwr_smk10, 1, order_by = fips),
         gwr_smk_lag2 = lag(gwr_smk10, 2, order_by = fips),
         gwr_smk_lag3 = lag(gwr_smk10, 3, order_by = fips),
         gwr_smk_lag4 = lag(gwr_smk10, 4, order_by = fips),         
         gwr_smk_lag5 = lag(gwr_smk10, 5, order_by = fips),
         gwr_smk_lag6 = lag(gwr_smk10, 6, order_by = fips)) %>% 
  select(-gwr_smk)

# join with morbidity data
ed_ts_lag <- total_ed %>%
  select(-gwr_smk10) %>% 
  left_join(gwr_lag, by = c("fips", "date")) %>% 
  filter(complete.cases(.)) 

# output matrix of smoke and lag variables
smk_matrix <- as.matrix(select(ed_ts_lag, gwr_smk10:gwr_smk_lag6))

# define custom functions for dl fit and dl results
# FIT function
dl_fit <- function(data, exp_mat, outcome, 
                   lag_df, pm_spline) {
  # create list of lengh pm_df
  mod_aic <- sapply(lag_df, function(x){
    pm_b <- splines::ns(0:(ncol(exp_mat)-1), df=x, intercept=T)
    # create pm basis
    pm_basis <- exp_mat %*% pm_b
    # run model
    mod <- glmer(as.formula(paste0(outcome, 
                                   "~pm_basis + pm_spline + (1|fips) + offset(log(population))")),
                 data, family="poisson"(link="log"), 
                 control = glmerControl(optimizer = "bobyqa"))
    aic <- AIC(mod)
    return(aic)
  }) # end sapply
  # create aic df
  pm_df_aic <- data.frame(outcome, lag_df, mod_aic)
  return(pm_df_aic)
}

# defining exposure spline with 13 df
pm_spl <- ns(ed_ts_lag$date, df = 13)
# find aic fit for each outcome
dl_fit_df <- outcomes %>% 
  map_dfr(~dl_fit(data=ed_ts_lag, exp_mat = smk_matrix, outcome = .,
                  lag_df = 2:5, pm_spline = pm_spl)) %>% 
  mutate(outcome = forcats::fct_relevel(outcome, outcomes))

# distributed lag spline fit ----

# using a pm spline with 13 degrees of freedom
# defining exposure spline with 13 df
pm_spl <- ns(ed_ts_lag$date, df = 13)
# find aic fit for each outcome
dl_fit_df <- outcomes %>% 
  map_dfr(~dl_fit(data=ed_ts_lag, exp_mat = smk_matrix, outcome = .,
                  lag_df = 2:5, pm_spline = pm_spl)) %>% 
  mutate(outcome = forcats::fct_relevel(outcome, outcomes))

# filter to minimum aic degrees of freedom minimum aic 
min_aic <- dl_fit_df %>% 
  group_by(outcome) %>% 
  slice(which.min(mod_aic))

# save dl_fit_df
write_csv(dl_fit_df, "./data/health/2015-ts_morbidity_spline_dl_fit.csv")

# distributed lag results -----
# distributed lag function
distributed_lag <- function(data, exp_mat, outcome, 
                            lag_df, pm_spline){
  # define basis b using natural spline function
  pm_b <- ns(0:(ncol(exp_mat)-1), df = lag_df, intercept = T)
  # multiply lagged pm matrix by basis
  pm_basis <- exp_mat %*% pm_b
  # fit mixed model
  mod <- glmer(as.formula(paste0(outcome, 
      "~pm_basis + pm_spline + weekend + (1|fips) + offset(log(population))")),
               data, family="poisson"(link="log"), 
               control = glmerControl(optimizer = "bobyqa"))
  
  # calculate estimates ----
  # output pm basis parameters
  dl_parms <- broom::tidy(mod) %>% 
    filter(stringr::str_detect(term, "pm_basis")) %>% 
    select(estimate) %>% 
    as_vector()
  
  # estimate distributed lag values for each day
  estimate <-  pm_b %*% dl_parms
  # time variable
  time <- ((rep(1:length(estimate))-1))
  # covariance matrix for knots 
  # fix the matrix
  cov_matrix <- as.matrix(vcov(mod))[2:(lag_df+1), 2:(lag_df+1)]
  # estimate variance of spline
  variance <- pm_b %*% cov_matrix %*% t(pm_b)
  # estimate lag ----
  # estimate standard error for each lag day
  l_se <- sqrt(diag(variance))
  # calculate lower and upper bound
  l_lower95 <- estimate+(l_se*qnorm(1-0.975))
  l_upper95 <- estimate+(l_se*qnorm(0.975))
  l_type <- as.character(rep("lag", times = length(estimate)))
  # l_estimate
  l_estimate <- data.frame(outcome, l_type, time, exp(estimate), 
                           exp(l_lower95), exp(l_upper95)) 
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
  cumulative_se <- sapply(seq_along(estimate), function(x){
    sqrt(sum(variance[1:x,1:x]))
  })
  # cumulative 95CI
  c_lower95 <- cumulative_estimate+(cumulative_se*qnorm(1-0.975))
  c_upper95 <- cumulative_estimate+(cumulative_se*qnorm(0.975))
  # return dataframe
  c_estimate <- data.frame(outcome, c_type, time,
                           exp(cumulative_estimate), exp(c_lower95), exp(c_upper95))
  # assign column names
  colnames(c_estimate) <- c("outcome", "type","time","estimate",
                            "lower95","upper95") 
  
  # bind lag and cumulative estimates together
  return_estimate <- rbind(c_estimate, l_estimate)
  # return estimate  
  return(return_estimate)
} # end of function


# distributed lag results for presentation ------
# defining exposure spline 
pm_spl <- ns(ed_ts_lag$date, df = 13)
# create outcome vector and degrees of freedom vector to use map2 function
outcome_vector <- min_aic$outcome
df_vector <- min_aic$lag_df
# estimate cumulative results
dl_results <- map2(.x=outcome_vector, .y=df_vector,
  ~distributed_lag(data=ed_ts_lag, exp_mat = smk_matrix, outcome = .x, 
                    lag_df = .y, pm_spline = pm_spl)) %>% 
  plyr::rbind.fill()

# save results so I don't have to run them everytime i want to knit doc
write_csv(dl_results, "./data/health/2015-morbidity_dl_results.csv")

