# ------------------------------------------------------------------------------
# Title: Proportion Intersection of Western US Grid 1015 and County Shapefiles
# Author: Ryan Gan
# Date Created: 2018-02-16
# ------------------------------------------------------------------------------

# Purpose: This script was designed to calculate the proportion interest on the
# 2010 to 2015 grid provided by Kate O'Dell and Bonne Ford. I plan to run this 
# on the ozone server, so have written it in a script rather than markdown doc.

# libraries ----
library(tidyverse)
library(sf)

# read sf objects ----
# grid path
grid_path <- paste0("./data/shapefile/1016-wus_grid")
# read grid sf
grid_sf <- st_read(dsn = grid_path, layer = "1016-wus_grid") 

# county path
county_path <- paste0("./data/shapefile/us_county")
# read county sf
county_sf <- st_read(dsn = county_path, layer = "us_county")

# assign crs of county_sf to grid_sf
wgs84 <- st_crs(county_sf)
# assign crs
st_crs(grid_sf) <- wgs84

# subset western county sf
# STATE FIPS I want to subset
western_state_fp <- c("04","06","08","16","53","41","32","49","56",
                      "35", "30")

# subset simple features
western_county_sf <- county_sf %>% 
  filter(STATEFP %in% western_state_fp) %>% 
  mutate(FIPS = paste0(STATEFP, COUNTYFP))


# propotion intersetc ----
# create empty dataframe
prop_int_tibble <- grid_sf$GRID_ID %>% 
  tibble() %>% 
  rename(grid_id = ".")

# start compute time
start_time <- Sys.time()

for(i in 1:length(western_county_sf$FIPS)){
  # subset county to find intersect
  county <- slice(western_county_sf, i)
  # extract fips number for variable name
  fips_id <- paste0("fips_", county$FIPS)
  # subset grid cells that touch any part of the county
  grid <- grid_sf[county,]
  # subset the intersected polygon
  inter_poly <- st_intersection(grid, county) %>% 
    # filter only to polygon or multipolygon type 
    # to avoid errors with point or line types
    filter(st_is(., c("POLYGON", "MULTIPOLYGON"))==T)
  # filter grid ids to only grids in the inter_poly object
  grid2 <- grid %>% filter(GRID_ID %in% inter_poly$GRID_ID) 
  # find proportion intersect with original grid
  prop_int <- as.numeric(st_area(inter_poly)/st_area(grid2))
  # subset grid id
  grid_id <- grid2$GRID_ID
  # make a tibble
  county_grid_int <- tibble(grid_id, prop_int) %>% 
    set_names(c("grid_id", fips_id))
  # join with full tibble
  prop_int_tibble <- prop_int_tibble %>% 
    left_join(county_grid_int, by = "grid_id")
} # end loop

# compute time
stop_time <- Sys.time()
compute_time <- stop_time - start_time
# run time
compute_time
# 5.5 minutes

# Check of LA county ----
# subset la county (FIPS 06037)
la_county <- western_county_sf %>% 
  filter(FIPS == "06037")

# subset proportion intersect values for LA county
la_prop_int <- prop_int_tibble %>% 
  select(grid_id, fips_06037) %>% 
  filter(!is.na(fips_06037))

# subset la grid
la_grid <- grid_sf %>% 
  filter(GRID_ID %in% la_prop_int$grid_id) %>% 
  left_join(la_prop_int, by = c("GRID_ID" = "grid_id")) %>% 
  rename(prop_int = fips_06037)

# make a dataframe with rounded proportion values
grid_names <- la_grid %>% 
  mutate(proportion = round(prop_int,2)) %>% 
  group_by(GRID_ID) %>% 
  mutate(lon = unlist(st_centroid(geometry))[1],
         lat = unlist(st_centroid(geometry))[2]) %>% 
  ungroup() %>% 
  select(GRID_ID, lon, lat, proportion)
# convert simple featrues to dataframe
st_geometry(grid_names) <- NULL

# plot grid overlayed on county
plot <- ggplot(la_county) + 
  geom_sf() +
  geom_sf(data=la_grid, aes(fill=prop_int), alpha=0.7) +
  geom_text(data=grid_names, aes(x=lon, y=lat, label=proportion), 
            size = 2.5, angle=45) +
  theme_bw()

plot
# plot looks good

# writing final proportion intersect file ----
# replace NAs with 0s
wus_prop_int <- prop_int_tibble %>% 
  rename(GRID_ID = grid_id) %>% 
  mutate_all(funs(replace(., is.na(.), 0)))

# save file
write_csv(wus_prop_int, paste0("./data/smoke/1016-countygrid_propint_us.csv"))

