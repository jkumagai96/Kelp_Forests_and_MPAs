# Date: November 23rd 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Rasterize MHWs and downsample
# BIO 202: Ecological Statistics

##### Set up ###################################################################
# load packages
library(tidyverse)
library(raster)

# load data
MHW_df <- readRDS("Processed_data/SST/MHW_1984_2021.rds")
base_grid <- raster("Data/standard_grid.tif")

##### Format data ##############################################################
MHW_intensity <- MHW_df %>% 
  arrange(year) %>% 
  dplyr::select(-c(mhw_events, mhw_days)) %>% # Filter data
  pivot_wider(names_from = year,
              values_from = mhw_int_cumulative) # rasterfrom XYZ requires wide data
  

##### Process data #############################################################
# Step 1: Create MHW rasters at 0.25 degree grid
MHW_025_grid <- rasterFromXYZ(MHW_intensity, 
                              crs = 4326) 

# Step 2: Downsample to base_grid 
MHW_001_grid <- disaggregate(MHW_025_grid, fact = 25)

# Step 3: Create table that can be joined by lat long
MHW_001_table <- rasterToPoints(MHW_001_grid) %>% 
  as.data.frame() %>% 
  mutate_all(~replace(., is.na(.), 0)) # Set NAs to zero where there was no MHWs detected

# Format final data 
final_MHW_data <- MHW_001_table %>% 
  rename("long" = x, "lat" = y) %>% 
  pivot_longer(cols = X1984:X2021, names_to = "year", values_to = "MHW_cummulative") %>% 
  mutate(year = substring(year, 2))

##### Export ###################################################################
saveRDS(object = final_MHW_data,
        file = "Processed_data/SST/MHW_cummulative_intensity_1km.rds")
