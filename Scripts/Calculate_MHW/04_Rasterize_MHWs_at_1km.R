# Date: November 23rd 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Rasterize MHWs, Cold spells and downsample
# BIO 202: Ecological Statistics

##### Set up ###################################################################
# load packages
library(tidyverse)
library(raster)

print('packages loaded')

# load data
MHW_df <- readRDS("Processed_data/SST/MHW_1983_2021.rds")
CS_df <- readRDS("Processed_data/SST/CS_1983_2021.rds")
base_grid <- raster("Data/standard_grid.tif")

print('data loaded')

# First we will go through formatting, processing, and exporting for marine heat
# waves, then the last section will do the same process for cold spells.
##### Format data ##############################################################
MHW_intensity <- MHW_df %>% 
  arrange(year) %>% 
  relocate(lon, lat) %>% 
  dplyr::select(-c(total_days)) %>% # Filter data
  pivot_wider(names_from = year,
              values_from = total_icum) # rasterfrom XYZ requires wide data
  

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
  pivot_longer(cols = X1983:X2021, names_to = "year", values_to = "MHW_cummulative") %>% 
  mutate(year = substring(year, 2))

##### Export ###################################################################
# Export tif to visualize in QGIS 
writeRaster(MHW_001_grid, "Processed_data/SST/MHW_cummulative_intensity_1km.tif", overwrite = T) 

saveRDS(object = final_MHW_data,
        file = "Processed_data/SST/MHW_cummulative_intensity_1km.rds")

print('marine heat wave data exported')

##### Cold Spells ##############################################################
CS_intensity <- CS_df %>% 
  arrange(year) %>% 
  relocate(lon, lat) %>% 
  dplyr::select(-c(total_days)) %>% # Filter data
  pivot_wider(names_from = year,
              values_from = total_icum) # rasterfrom XYZ requires wide data

# Step 1: Create cold spells rasters at 0.25 degree grid
CS_025_grid <- rasterFromXYZ(CS_intensity, 
                              crs = 4326) 

# Step 2: Downsample to base_grid 
CS_001_grid <- disaggregate(CS_025_grid, fact = 25)

# Step 3: Create table that can be joined by lat long
CS_001_table <- rasterToPoints(CS_001_grid) %>% 
  as.data.frame() %>% 
  mutate_all(~replace(., is.na(.), 0)) # Set NAs to zero where there was no cold spells are detected

# Format final data 
final_CS_data <- CS_001_table %>% 
  rename("long" = x, "lat" = y) %>% 
  pivot_longer(cols = X1983:X2021, names_to = "year", values_to = "CS_cummulative") %>% 
  mutate(year = substring(year, 2))

saveRDS(object = final_CS_data,
        file = "Processed_data/SST/CS_cummulative_intensity_1km.rds")

print('marine cold spell data exported')
print('script is finished')
