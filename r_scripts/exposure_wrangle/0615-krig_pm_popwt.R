# ------------------------------------------------------------------------------
# Title: County-level population-weighting of kriged PM2.5 values 2006-2015
# Author: Ryan Gan
# Date Created: 2018-02-23
# ------------------------------------------------------------------------------

# This script population-weights gridded estimates of PM2.5 to county for the
# continental united states

# libraries ----
library(tidyverse) # all purpose tidy data
library(ncdf4) # working with netcdf files
library(parallel) # for parallel; computing not sure i need this yet

# proportion intersect setup ------
# reading in files and creating necissary vectors outside of apply functions

# read in proportion intersect grid 
prop_int <- read_csv("./data/smoke/1016-countygrid_propint_us.csv", 
                     col_types = cols(.default = "c")) %>% 
  mutate_all(funs(as.numeric))

# proportion intersect matrix
pi_mat <- as.matrix(prop_int[,2:415])
# assign grid name as row name
rownames(pi_mat) <- prop_int$GRID_ID

# read in grid population density ---
popden2015 <- read_csv("./data/smoke/2015-popdensity_us_grid.csv") %>% 
  arrange(GRID_ID) %>% 
  # set missing to 0
  mutate(popden2015 = ifelse(is.na(popden2015), 0, popden2015))

# population density vector
pop_vec <- popden2015$popden2015

# estimating county population density and inverse for use in formula
# estimate of population density per county
popden_county <- t(pi_mat) %*% pop_vec
# calcuate the inverse of county population density
popden_county_inverse <- 1/popden_county

# extract fips names
fips <- substr(colnames(prop_int)[-1], 6, 10)

# open nc connection to pm grid -----
# open connection pm grid
pm_nc_path <- paste0("./data/smoke/krigedPM25_06-15/krigedPM25_2015.nc")
pm_nc <- nc_open(pm_nc_path)

# extract latitude and longitude values
pm_lat <- ncvar_get(nc = pm_nc, varid = "lat")
pm_lon <- ncvar_get(nc = pm_nc, varid = "lon")
# create cell id values
id <- seq(1:58401)
# bind columns together for spatial coordinates of pm grid and cell id
pm_coords <- data.frame(cbind(id, as.vector(pm_lon), as.vector(pm_lat))) %>% 
  rename(cell_id = id, lon = V2, lat = V3)

# create Cell_ID key that matches the shapefile ----
cell_id_mat <- t(apply(matrix(58401:1, nrow=189, ncol=309, byrow=T),1,rev))
# check lower right section of the matrix (should go lower left to right)
cell_id_mat[180:189,1:10]
# vector of cell id
cell_id <- as.vector(cell_id_mat)

# GRID_ID matrix based on cell location
grid_id_mat <- apply(matrix(1:58401, nrow=189, ncol=309), 2, rev)
# check lower right seciton of matrix (should go lower left up)
grid_id_mat[180:189,1:10]
# vector of GRID_ID
GRID_ID <- as.vector(grid_id_mat)

# create grid id key
grid_id_key <- data.frame(cbind(cell_id, GRID_ID)) %>% 
  left_join(pm_coords, by = "cell_id") %>% 
  arrange(cell_id)

# prep to make sure grids line up correctly
# extract vector of arranged WRFGRID IDs that match cell IDs for rownames
GRID_ID <- arrange(grid_id_key, id)$GRID_ID

# i think the loop should start here
# output 3d arrays from .nc files ----
krig_pm <- ncvar_get(pm_nc, varid = "PM25")
bg_pm <- ncvar_get(pm_nc, varid = "Background PM25")
hms_smk <- ncvar_get(pm_nc, varid = "HMS Smoke")

# kates date in the .nc file is just one value; i'm just going to create
# my own date vector that will account for leap year too
year <- "2015"
date <- paste0("d", gsub("-", "", 
  seq.Date(from = as.Date(paste0(year, "-01-01")), 
    to = as.Date(paste0(year, "-12-31")), by = "day"))) # end date seq

# convert pm2.5 to 2d matrix
pm_mat <- apply(krig_pm, 3, as.vector)
bg_pm_mat <- apply(bg_pm, 3, as.vector)
hms_mat <- apply(hms_smk, 3, as.vector)

# bind in to a list to itterate through population weighting function
pm_list <- list(pm_mat, bg_pm_mat, hms_mat)
# assign names to list elements
names(pm_list) <- c("krig_pm", "background_pm", "hms_smk")

# lapply row and column names and sort by GRID_ID Cell
pm_list <- lapply(pm_list, function(x){
  colnames(x) <- date
  rownames(x) <- GRID_ID
  x <- x[order(as.numeric(row.names(x))),]
})

pm_list[[1]][1:10,1:10]

popwt_pm_list <- lapply(pm_list, function(x){
  # multiply population vector by pm matrices
  pm_pop_mat <- pop_vec * x
  # multiply pm population matrix by county/grid intersection matrix 
  # produces a summed estimate of pm values in a county for each day
  county_pm_sum_mat <- t(pi_mat) %*% pm_pop_mat
  # multiple county summed PM values by inverse of county population density
  # produces county population weighted PM estimates for each day
  county_popwt_pm <- county_pm_sum_mat * as.vector(popden_county_inverse)
  
  # create dataframe and join in FIPS id
  pm_county_wt <- as_data_frame(cbind(fips, county_popwt_pm)) %>% 
    # convert characters to numeric
    mutate_at(vars(-fips), funs(as.numeric))
  # return dataframe
  return(pm_county_wt)
})


# population weighting for each 
# multiply population vector by regridded environmental vector
e_pop_vec <- pop_vec * pm_mat
# multiply population matrix by county/grid intersection matrix 
# produces a summed estimate of temp values 
county_e_sum_mat <- t(pi_mat) %*% e_pop_vec
# multiple county summed temp values by inverse of county population density
# produces county population weighted temp estimates f
county_popwt_e <- county_e_sum_mat * as.vector(popden_county_inverse)

# hms smoke
s_pop_vec <- pop_vec * hms_mat
c_s_sum_mat <- t(pi_mat) %*% s_pop_vec
county_popwt_s <- c_s_sum_mat * as.vector(popden_county_inverse)

county_popwt_s[1:10,1:10]

# side note: one day I should put all these estimates in to a database
