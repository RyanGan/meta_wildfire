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
# open connection to temp netcdf (may eventually itterate this for different years)
temp_nc <- nc_open("./data/smoke/air.2m.2015.nc")
temp_nc
# extract lat/lon coords from temperature nc
# extract temp lat lon
t_lat <- ncvar_get(nc = temp_nc, var = "lat")
t_lon <- ncvar_get(nc = temp_nc, var = "lon")
# extract air temp for the year
t_temp <- ncvar_get(nc  =temp_nc, var = "air")

# set temp coordinates as matrix
temp_coords <- matrix(cbind(as.vector(t_lon), as.vector(t_lat)), ncol=2)
# set new extent
# i need to rebuild temp raster with appropriate lat lon coords
# set up empty extent
e <- extent(temp_coords[,1:2])
# create empty raster of point centroids I want to assign temp values to
r <- raster(e, nrow=277, ncol=349)
# new temp raster with latlon coords 
r_temp <- rasterize(temp_coords[,1:2], r, as.vector(t_temp[,,1]), fun=mean)
# extract values to pm points
temp_k_regrid <- extract(r_temp, pm_coords[,2:3])

# apply to each date
test_mat <- t_temp[,,1:2]

# apply function and time it
system.time(
  regrid_temp_mat <- apply(t_temp, 3, function(x){
    rast_temp <- rasterize(temp_coords[,1:2], r, as.vector(x), fun=mean)
    temp_k_regrid <- extract(rast_temp, pm_coords[,2:3])
    return(temp_k_regrid)
})
)
# i think it might be a good idea to population-weight all in the same script
# and save the county-specific daily temperature estimates in a .nc file to save
# space. Took about 26 minutes to regrid 1 year of temp data on my local machine
# so for 6 years it would take over 2 hours. 

summary(regrid_temp_mat)
