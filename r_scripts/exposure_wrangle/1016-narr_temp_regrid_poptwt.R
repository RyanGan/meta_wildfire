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
library(parallel) # for parallel computing

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

# create grid id key; I will bind regridded points to this raster
grid_id_key <- data.frame(cbind(cell_id, GRID_ID)) %>% 
  left_join(pm_coords, by = "cell_id") %>% 
  arrange(GRID_ID)

# regrid irregular temp raster data -----
# open connection to temp netcdf; opening to 2015 since grids are the same 
# across years
temp_nc <- nc_open("./data/smoke/air_temp_2m/air.2m.2015.nc")

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

# commented out since I figured this out
# # new temp raster with latlon coords 
# r_temp <- rasterize(temp_coords[,1:2], r, as.vector(t_temp[,,1]), fun=mean)
# # extract values to pm points
# temp_k_regrid <- extract(r_temp, pm_coords[,2:3])

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

# population weighting of temperature ------------------------------------------
# itterate through nc files, opening one at a time
temp_nc_files <- list.files("./data/smoke/air_temp_2m/")

# set up cluster of 6 cores to parallelize across each .nc file
cl <- makeCluster(6)

# load packages on each processor of the node/cluster
clusterCall(cl, function() c(library(tidyverse), library(ncdf4), 
                             library(raster)))

# export global sf objects and empty tibble to each core
clusterExport(cl, c("temp_nc_files", "r", "temp_coords", "pm_coords",
                    "pop_vec", "pi_mat", "popden_county_inverse", "fips"), 
              envir = .GlobalEnv)

system.time( # start system time
# parallel sapply nc read and write function over list of .nc files
parSapply(cl, temp_nc_files, function(meow){
  # open nc connection
  temp_nc <- nc_open(paste0("./data/smoke/air_temp_2m/", meow))
  # get temp_matrix
  temp_mat <- ncvar_get(temp_nc, varid = "air")
  # extract date values to assign to column header
  date <- paste0("d",
    gsub("-", "", as.Date(ncvar_get(temp_nc, varid = "time")/24, 
                 origin = "1800-01-01")))
  # extract year from name of list
  year <- substring(meow, 8, 11)

# apply function to extract and regrid temp and population weight to county
  # innner function to regrid and population weight    
    pop_wt_temp_mat <- apply(temp_mat, 3, function(x){
      # regrid temp raster 
      rast_temp <- rasterize(temp_coords[,1:2], r, as.vector(x), fun=mean)
      # extract values from raster to grid_id coordinates
      # I sorted coordinates by GRID_ID to match the population density and
      # proportion intersect ids
      temp_k_regrid <- extract(rast_temp, grid_id_key[,3:4])
      
      # population weigting
      # multiply population vector by regridded temp vector
      temp_pop_vec <- pop_vec * temp_k_regrid
      # multiply population matrix by county/grid intersection matrix 
      # produces a summed estimate of temp values 
      county_temp_sum_mat <- t(pi_mat) %*% temp_pop_vec
      # multiple county summed temp values by inverse of county population density
      # produces county population weighted temp estimates f
      county_popwt_temp <- county_temp_sum_mat * as.vector(popden_county_inverse)
      # return county population-weight temp
      return(county_popwt_temp)
    }
  ) # end inner function

# assign rownames and column names and create dataframe
temp_df <- data.frame(pop_wt_temp_mat)
# assign column names
colnames(temp_df) <- date
# create final dataframe to write
final_df <- cbind(fips, temp_df)

# write name 
write_name <- paste0("./data/smoke/air_temp_2m/county_air_temp/",
                     year,"-county_popwt_temp.csv")
# write final dataset to air temp file
write_csv(final_df, path = write_name)

# close nc connection
nc_close(temp_nc)

    } # end nc read function
  ) # end sapply
) # end system time

# close cluster
stopCluster(cl)

# check to see if files were created
list.files("./data/smoke/air_temp_2m/county_air_temp/", pattern = ".csv")
