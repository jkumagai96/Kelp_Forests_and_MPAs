# Date: April 5th 2023
# Author: Joy A. Kumagai (kumagaij@stanford.edu)
# Purpose: Summarize kelp connectivity 
# Kelp Forests, MPAs, and Heat waves Project

##### Set up: packages #########################################################
library(tidyverse)
library(raster)

##### Set up: Load Data ########################################################
connectivity_raw <- data.table::fread("Data/Kelp/connectedness_03022023.csv", header = T)
connectivity <- connectivity_raw

# Load standard .01 grid
base_grid <- read_rds("Data/standard_grid.Rda") 

##### Summarize Data ###########################################################
# Format Data
colnames(connectivity)[1] <- "lat"
colnames(connectivity)[2] <- "lon"

# Tibble of coordinates
coords <- tibble(x = connectivity$lon, y = connectivity$lat)

# Remove coordinates from dataframe
connectivity <- connectivity[,-c(1,2)]

# Rasterize to match base grid 
k <- rasterize(x = coords,                 
               y = base_grid,              
               field = connectivity,    
               fun = mean, # Ask Nur 
               na.rm = T)

connectivity_per_quarter <- rasterToPoints(k) %>% 
  as.data.frame() %>%   
  dplyr::select(-c(X2022.1, X2022.2, X2022.3))#  %>% 
  #filter(y <= 42)

connectivity_quarter_long <- connectivity_per_quarter %>% 
  pivot_longer(X1984.1:X2021.4, names_to = "time", values_to = "connectivity") %>% 
  arrange(time) %>% 
  mutate(year = str_sub(time, -6, -3)) %>% 
  mutate(quarter = str_sub(time, -1)) %>% 
  rename(lon = x, 
         lat = y)
# 33% of the data are NA's

##### Export ###################################################################
write.csv(connectivity_quarter_long, 
          "Processed_data/data_tables/connectivity_per_quarter.csv", 
          row.names = F)
