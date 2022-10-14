# Date: October 14th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Extract raster data into table and create pixel ID dataset
# BIO 202: Ecological Statistics

##### Set up: packages #########################################################

# Load packages
library(raster)
library(tidyverse)
library(here)
library(sf)

##### Create Pixel ID spatial dataset ##########################################
# load depth data 
depth_r <- raster(here(
  "Processed_data", 
  "kelp_variables",
  "LandsatKelp_Quarterly_depth.tif"
))

# Create points out of depth raster
depth_p <- rasterToPoints(depth_r) %>% 
  as.data.frame %>% 
  mutate(PixelID = 1:nrow(.))

# export
write.csv(depth_p, 
          "Processed_data/data_tables/PixelID_reference.csv",
          row.names = F
)

##### Create Datables for each variable ########################################

area_files <- list.files("Processed_data/kelp_variables", pattern = "area_20", full.names = T)
years <- 2009:2021
PixelID <- 1:nrow(depth_p)

# load data 
r <- raster(area_files[1])

# convert to df
area_data <- rasterToPoints(r) %>% 
  as.data.frame() %>% 
  rename("area" = "layer") %>% 
  mutate(year = years[1],
         PixelID = PixelID) %>% 
  select(x, y, PixelID, year, area)

# Need to make a for loop that ends up with all the years in this data frame

# All tested rasters have the same dimensions, let's try to join them the proper way. Cbind also works 
new_dat2 <- cbind(as.data.frame(test), as.data.frame(test3)$depth)


