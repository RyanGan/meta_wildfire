# ------------------------------------------------------------------------------
# Title: Regridding of NCEP North American Regional Reanalysis Temperature to 
#         Colorado Grid Corridnates
# Author: Ryan Gan
# Date Created: 2018-02-21
# ------------------------------------------------------------------------------

# This script regrids temperature estimates from NCEP NARR to Kate's PM grid
# for colorado

# libraries ----
library(tidyverse) # all purpose tidy data
library(ncdf4) # working with netcdf files
library(raster) # creating rasters
library(parallel) # for parallel computing

# read in colorado pm coordinates ----
grid_coords <- read_csv("./data/smoke/colo_krig_grid_coords.csv")

# regrided temperature raster ----
# open temp nc
temp_nc <- nc_open(paste0("./data/smoke/air_temp_2m/air.2m.2015.nc"))
# extract temp lat lon
t_lat <- as.vector(ncvar_get(nc = temp_nc, var = "lat"))
t_lon <- as.vector(ncvar_get(nc = temp_nc, var = "lon"))
# temp_coordinates as matrix
temp_coords <- as.matrix(cbind(t_lon, t_lat), ncol =2)
# set names to lat and lon
colnames(temp_coords) <- c("lon", "lat")
# view first obs of temp coordinates
head(temp_coords)
# create extent object
e <- extent(temp_coords[,1:2])
# create new raster
r <- raster(e, nrow=277, ncol=349)

# parallel computation -----
#define air_temp file path
temp_path <- paste0("./data/smoke/air_temp_2m/")
# create list of temp files
temp_nc_files<- list.files(temp_path)
# print temp files
print(temp_nc_files)

# set up cluster of 6 cores to parallelize across each .nc file 
cl <- makeCluster(6)

# load packages on each processor of the node/cluster
clusterCall(cl, function() c(library(tidyverse), library(ncdf4), 
                             library(raster)))

# export global sf objects and empty tibble to each core
clusterExport(cl, c("temp_nc_files", "r", "temp_coords", "grid_coords", 
                    "temp_path"), 
              envir = .GlobalEnv)

# parallel sapply nc read and write function over list of .nc files
parSapply(cl, temp_nc_files, function(meow){
  
  # open nc connection
  temp_nc <- nc_open(paste0(temp_path, meow))
  # get temp_matrix
  temp_mat <- ncvar_get(temp_nc, varid = "air")
  # extract date values to assign to column header
  date <- paste0("d",
    gsub("-", "", as.Date(ncvar_get(temp_nc, varid = "time")/24, 
      origin = "1800-01-01")))
  
  # extract year from name of list
  year <- substring(meow, 8, 11)
  
    # regrid apply for each day
    co_temp_mat <- apply(temp_mat, 3, function(x){
    # regrid temp raster 
    rast_temp <- rasterize(temp_coords[,1:2], r, as.vector(x), fun=mean)
    # extract values from raster to grid_id coordinates
    temp_k_regrid <- extract(rast_temp, grid_coords[,2:3])
    }) # end innter apply function

  # assign rownnames and colunn names and create dataframe
  temp_df <- data.frame(co_temp_mat)
  # assign column names of dates
  colnames(temp_df) <- date
  # output grid ids
  GRID_ID <- grid_coords$GRID_ID
  # bind in grid_id column
  temp_df <- cbind(GRID_ID, temp_df)
  
  
  # print first obs of temp_df
  print(temp_df[1:5, 1:5])
  
  # write name 
  write_name <- paste0("./data/smoke/co_air_temp/",
                       year,"-co_grid_temp.csv")
  
  # write final dataset to air temp file
  write_csv(temp_df, path = write_name)
  
  # close nc connection
  nc_close(temp_nc)
}) # end parsapply

# close cluster
stopCluster(cl)

# check to see if files were created
list.files("./data/smoke/co_air_temp/", pattern = ".csv")

