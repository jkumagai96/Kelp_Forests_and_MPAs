# Date: October 14th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Extract raster data into table 
# BIO 202: Ecological Statistics

##### Set up: packages #########################################################

# Load packages
library(raster)
library(here)
library(sf)
library(tidyverse)

##### Create Pixel ID spatial dataset ##########################################
# load depth data 
depth_r <- raster(here(
  "Processed_data", 
  "kelp_variables",
  "LandsatKelp_Quarterly_depth.tif"
))

# Create points out of depth raster
depth_p <- rasterToPoints(depth_r) %>% 
  as.data.frame() %>% 
  filter(y <= 37.4) %>% 
  mutate(PixelID = 1:nrow(.)) %>% 
  rename(depth = layer)

# export
write.csv(depth_p, 
          "Processed_data/data_tables/PixelID_reference.csv",
          row.names = F
)

##### Create a datatable with all variables ####################################
## Initialize global variables 
years <- 2009:2021                # years of the data we are working with 
PixelID <- 1:nrow(depth_p)        # ID just in case we need to join variables not in order
kelp_list <- list()               # Final list all data will be in 

# Variables we are working in that will go through each iteration of the 
# outer loop (j) 
variables <- c("area", "biomass", "hsmax", "nitrate", "temperature")

# Outer loop goes through variables 
for (j in 1:length(variables)) {
  variable_files <- list.files("Processed_data/kelp_variables", pattern = variables[j], full.names = T)
  variable_all <- c()
  
  # Inner loop goes through each year of data for the variable 
  for (i in 1:length(years)) {
    # load data 
    r <- raster(variable_files[i])
    
    # convert to df
    data_yr <- rasterToPoints(r) %>% 
      as.data.frame() %>% 
      filter(y <= 37.4) %>% 
      mutate(year = years[i],
             PixelID = PixelID) %>% 
      dplyr::select(x, y, PixelID, year, layer)
    
    # Adds all the year data together to the bottom of the data frame 
    variable_all <- rbind(variable_all, data_yr) 
    
  } # end of inner loop 
  
  
  # Assigns the name of the column (variable)
  colnames(variable_all)[5] <- variables[j]     
  
  # Adds each variable into the final list 
  kelp_list[[j]] <- variable_all                

  }


# Now we put each part of the list together so we have a data frame with 
# Latitude, longitude, Pixel ID, year, area, biomass, hsmax, nitrate, and temperature

# sets up the table with area dataset from the list 
area <- left_join(kelp_list[[1]], depth_p) 

# adds each variable to the table 
biomass <- kelp_list[[2]]$biomass
hsmax <- kelp_list[[3]]$hsmax
nitrate <- kelp_list[[4]]$nitrate
temperature <- kelp_list[[5]]$temperature

# Putting the dataset together and formatting 
kelp_all_data <- cbind(area, biomass, hsmax, nitrate, temperature) %>% 
  rename("lon" = "x", 
         "lat" = "y") %>% 
  relocate(depth, .after = PixelID)

##### Export ###################################################################
write.csv(kelp_all_data, "Processed_data/data_tables/kelp_data.csv", row.names = F)

# End of script
