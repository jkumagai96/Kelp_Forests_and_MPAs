# Date: November 28th 2022
# Author: Joy A. Kumagai (kumagaij@stanford.edu)
# Purpose: Summarize the data by year so it is easier to work with 
# BIO 202: Ecological Statistics

##### Set up ###################################################################
# Load packages
library(tidyverse)

# Load Data
all_kelp_data <- read.csv("Processed_data/data_tables/kelp_data_w_mpas_per_quarter.csv")
distances <- read.csv("Processed_data/distances_to_coast.csv") %>% 
  dplyr::select(-depth) %>% 
  rename(long = x, lat = y)

###### Data Manipulation and formating #########################################
# Summarize by year
kelp_data_yr <- all_kelp_data %>% 
  group_by(PixelID, year) %>% 
  summarize(area = mean(area, na.rm = T),
            biomass = mean(biomass, na.rm = T),
            hsmax = mean(hsmax, na.rm = T),
            nitrate = mean(nitrate, na.rm = T),
            temperature = mean(temperature, na.rm = T),
            MHW_intensity = first(MHW_intensity), # same value for all quarters as it was calculated by year
            CS_intensity = first(CS_intensity), # same value for all quarters 
            mpa_status = first(mpa_status), # same value for all quarters 
            mpa_area = first(mpa_area),  # same value for all quarters 
            region = first(region)) # same value for all quarters 

kelp_data_yr$mpa_area[is.na(kelp_data_yr$mpa_area)] <- 0

# Subset the position details from the original data (mpa status, lat, long, Pixel ID, human gravity, and depth)
subset <- all_kelp_data %>% 
  dplyr::select(long, lat, PixelID, depth, gravity, Mpa_ID) %>% 
  unique(.) %>% 
  arrange(PixelID)

# Add this information back in
df <- left_join(kelp_data_yr, subset, by = "PixelID") 

##### Add in distance to coast #################################################

final <- left_join(df, distances, by = c("long", "lat", "PixelID"))

##### Export ###################################################################

write.csv(final, "Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv", 
          row.names = F)
