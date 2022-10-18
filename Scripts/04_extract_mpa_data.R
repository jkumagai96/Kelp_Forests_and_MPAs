# Date: October 16th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Create a data table with MPA protection per point with kelp data
# BIO 202: Ecological Statistics

##### Set up: packages #########################################################
# load packages 
library(sf)
library(tidyverse)

##### Set up: load data ########################################################
mpas_original <- read_sf("Data/Filtered_MPAs/MPAs_CA_Giant_Kelp.shp")
station_data <- read.csv("Processed_data/data_tables/PixelID_reference.csv")
kelp_data <- read.csv("Processed_data/data_tables/kelp_data.csv")

##### Formatting ###############################################################
# Select needed attribuets from MPAs
mpas <- mpas_original %>% 
  dplyr::select(Site_ID_12, Estab_Yr_1, AreaMar_12, Status_12)  %>% # I selected Area_Mar12, and not Area? QUESTION
  rename("mpa_status" = "Status_12") %>% 
  # Match projections with kelp data 
  st_transform(crs = 4326) 

# convert station points into spatial data to intersect with mpas 
station_points <- st_as_sf(station_data, 
                           coords = c("x", "y"), 
                           crs = 4326) 

# Make sure the projections match 
crs(station_points) == crs(mpas)


##### Processing ###############################################################
# Spatial intersect between mpas and points to get all points within mpas 
points_in_mpas <- st_intersection(station_points, mpas) %>% 
  dplyr::select(PixelID, mpa_status, Site_ID_12) %>% 
  st_drop_geometry()

# Join data with kelp data 
kelp_w_mpas <- left_join(kelp_data, points_in_mpas, by = "PixelID") 

# Format values of mpa_status so NA is none, no take is full, and uniform multiple use is partial
kelp_w_mpas$mpa_status[is.na(kelp_w_mpas$mpa_status)] <- "None"

# Formatting findal data 
final_data <- kelp_w_mpas %>% 
  rename("Mpa_ID" = "Site_ID_12") %>% 
  relocate(Mpa_ID, .after = depth) %>% 
  relocate(mpa_status, .after = Mpa_ID)

##### Export ###################################################################

write.csv(final_data, "Processed_data/data_tables/kelp_data_w_mpas.csv", row.names = F)
# End of script 


