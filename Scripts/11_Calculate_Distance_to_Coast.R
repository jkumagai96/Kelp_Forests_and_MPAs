# Date: December 13th 2022
# Author: Joy A. Kumagai (kumagaij@stanford.edu)
# Purpose: Calculate distance to coast 
# Kelp Forests, MPAs, and Heat waves Project

##### Set up: packages #########################################################
# load packages 
library(sf)
library(tidyverse)

# load data
station_data <- read.csv("Processed_data/data_tables/PixelID_reference.csv")
CA <- read_sf("Data/gadm41_USA_shp/gadm41_USA_1.shp") %>% 
  st_transform(crs = 3310) %>%           # NAD83 / California Albers (equal area)
  filter(NAME_1 == "California")

##### Station Data to Points ###################################################
station_points <- st_as_sf(station_data, 
                           coords = c("x", "y"), 
                           crs = 4326) 
station_points <- st_transform(station_points, crs = 3310)

##### Calculate Distance #######################################################
distance_to_coast <- st_distance(station_points, CA)

##### Export distance per station point ########################################
data_final <- cbind(station_data, distance_to_coast)

write.csv(data_final, "Processed_data/distances_to_coast.csv", row.names = F)

