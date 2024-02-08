# Date: May 2nd 2023
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

##### Checking The Distribution ################################################
plot1 <- all_kelp_data %>% 
  ggplot(aes(x = year, y = area/1000000, color = quarter)) +
  geom_point() +
  labs(x = "Year", y = "Area (km2)")

# Export Figure
png("Figures/Distribution_of_data_by_quarter.png", width = 5, height = 5, units = "in", res = 600)
plot1
dev.off()

11817/241072 #4.9% of the data is missing when filtered for the correct regions

# Check which quarter is the max and min for each pixel, year combination

result <- all_kelp_data %>%
  dplyr::select(PixelID, year, area, quarter) %>% 
  group_by(PixelID, year) %>%
  slice(which.max(area)) %>%
  dplyr::select(PixelID, year, quarter)

# What percentage of the data is the maximum in each quarter 
result %>% filter(quarter == 'Q1') %>% nrow()/nrow(result) # 29%
result %>% filter(quarter == 'Q2') %>% nrow()/nrow(result) # 23%
result %>% filter(quarter == 'Q3') %>% nrow()/nrow(result) # 35%
result %>% filter(quarter == 'Q4') %>% nrow()/nrow(result) # 12% 

###### Removing years with too many NA's #######################################
n_nas_data <- all_kelp_data %>% 
  group_by(PixelID, year) %>% 
  summarize(c_nas = sum(is.na(area))) 

###### Data Manipulation and formating #########################################
# Summarize by year
kelp_data_yr <- all_kelp_data %>% 
  left_join(n_nas_data) %>% 
  filter(c_nas < 2) %>% # Removes pixels with too many NA quarters (removed 5% of data)
  group_by(PixelID, year) %>% 
  summarize(mean_area = mean(area, na.rm = T),
            min_area = min(area, na.rm = T),
            area = max(area, na.rm = T), # Taking the max value of area per year 
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

nrow(df)/(nrow(all_kelp_data)/4) # Filtering data so years with 2 or more quarters with NAs are removed, removed 5% of the data overall

##### Export ###################################################################

write.csv(df, "Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv", 
          row.names = F)
