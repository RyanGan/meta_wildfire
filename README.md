# Meta-Wildfire

## Repository Purpose

This repository contains code and some data elements for our project addressing the relationship between wildfire smoke exposure and cardiorespiratory morbidity and mortality in the Western United States. 

This repository does not include county-level mortality and morbidity time-series data as these data contain protected health information and are covered under the Health Information and Portability Acountability Act (HIPAA). 

## General Overview
Uses the state health outcomes and smoke data from Colorado, Oregon, and Washington for multiple fire seasons.

## Outcomes Variable Key
Each state time series of outcome data should have the following variables (make sure they are all lower case):

**fips**: 5 digit character variable where the first two digits are the state fips code and the last 3 digits are the county fips codes
**statefp**: 2 digit state fips code
**countyfp**: 3 digit county fips code
**date**: Date of counts of outcomes. %Y-%m-%d format. 
**n_obs**: Number of observations for that date and county
**

## To Do List:
- Population-weight county-level estimates of PM2.5
- Add 2015 Washington and 2015 Colorado, time series data to Washington 2012 and Oregon 2013 data.
- Pilot 2015 Colorado Mortality analysis for Rish.
- Rename files and organize repo.
- Change paths on files moved.

