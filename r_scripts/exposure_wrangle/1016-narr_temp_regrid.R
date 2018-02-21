# ------------------------------------------------------------------------------
# Title: Regridding of NCEP North American Regional Reanalysis Temperature
# Author: Ryan Gan
# Date Created: 2018-02-21
# ------------------------------------------------------------------------------

# This script regrids temperature estimates from NCEP NARR to Kate's PM grid

# libraries ----
library(tidyverse) # all purpose tidy data
library(ncdf4) # working with netcdf files
library(raster) # creating rasters

# open nc connections -----
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
cell_id_mat <- t(apply(matrix(58401:1, nrow=309, ncol=189, byrow=T),1,rev))
# check lower right section of the matrix (should go lower left to right)
cell_id_mat[300:309,1:10]
# vector of cell id
cell_id <- as.vector(cell_id_mat)

# GRID_ID matrix based on cell location
grid_id_mat <- apply(matrix(1:58401, nrow=309, ncol=189), 2, rev)
# check lower right seciton of matrix (should go lower left up)
grid_id_mat[300:309,1:10]
# vector of GRID_ID
GRID_ID <- as.vector(grid_id_mat)

# create grid id key
grid_id_key <- data.frame(cbind(cell_id, GRID_ID)) %>% 
  left_join(pm_coords, by = "cell_id") %>% 
  arrange(cell_id)

# checks commented out after confirmed it worked
# head(grid_id_key)  
# head(pm_coords)  
# # check to make sure grid id key and pm_coords columns are equal (except GRID_ID)
# all_equal(pm_coords, grid_id_key[,c(1,3:4)])
# grid_id_key can be used to join with shapefile

# regrid irregular temp raster data -----
# extract lat/lon coords from temperature nc

