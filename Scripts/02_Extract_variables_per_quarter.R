# Date: October 25th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Extract NetCDF data into table perquarter
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(here)
library(raster)
library(sf)
library(ncdf4)
library(stringr)

## Load Data 
# Load standard .01 grid
base_grid <- raster(here("Data", "standard_grid.tif"))

# Load all_data NC file
all_data <- nc_open(
  filename = here(
    "Data",
    "Kelp",
    "CAkelpCanopyEnv_2021.nc"),
  suppress_dimvals = F,
  write = F
)

##### Extraction ###############################################################

# Extract area values ----------------------------------------------------------
area <- ncvar_get(
  nc = all_data,
  varid = "area"
)

# Extract biomass values -------------------------------------------------------
biomass <- ncvar_get(
  nc = all_data,
  varid = "biomass"
)

# Extract temperature values ---------------------------------------------------
temperature <- ncvar_get(
  nc = all_data,
  varid = "temperature"
)

# Extract nitrate values -------------------------------------------------------
nitrate <- ncvar_get(
  nc = all_data,
  varid = "nitrate"
)

# Extract hsmax values ---------------------------------------------------------
hsmax <- ncvar_get(
  nc = all_data,
  varid = "hsmax"
)

# Extract lat and long values and depth ----------------------------------------

#Vector of latitudes
lat <- ncvar_get(
  nc = all_data,
  varid = "lat"
)

#vector of longitudes
lon <- ncvar_get(
  nc = all_data,
  varid = "lon"
)

# Extract depth values
depth <- ncvar_get(
  nc = all_data,
  varid = "depth"
)

# Tibble of coordinates
coords <- tibble(x = lon, y = lat)

# Extract time indicators -----------------------------------------------------
year <- ncvar_get(
  nc = all_data,
  varid = "year"
)

quarter <- ncvar_get(
  nc = all_data,
  varid = "quarter"
)

## Create depth raster #########################################################
# Create raster with depth data 
depth_r <- rasterize(x = coords,                 
                     y = base_grid,              
                     field = depth,    
                     fun = mean,
                     na.rm = T)

# Create points out of depth raster
depth_p <- rasterToPoints(depth_r) %>% 
  as.data.frame() %>% 
  filter(y <= 37.4) %>% 
  mutate(PixelID = 1:nrow(.)) %>% 
  rename(depth = layer)

##### Declare Functions ########################################################
# We want to filter out 1km2 cells where more than 25% of the pixels (30mx30m) 
# are NA's, so I came up with a function that returns a 1 if less than 25% of the 
# data are NAs, and NA if more than 25% of the data are NAs, thus this can be put 
# into the rasterize function to create a raster of 1's and NAs to then act as a 
# filter if multiplied by any other raster 

F_filter_NAs <- function(x,...) {
  len <- length(x) # includes NAs
  n_nas <- sum(is.na(x)) # calculate number of NAs
  if (n_nas/len < .25) {return(1)}
  else(return(NA))
}


test <- c(0, 1, 0, 0, 5, 0, 1, NA, NA, NA)
F_filter_NAs(test)

## Process variables ###########################################################
PixelID <- 1:nrow(depth_p)        # ID just in case we need to join variables not in order

study_years <- year
study_quarters <- quarter

variables <- c("biomass", "hsmax", "nitrate", "temperature")
var_list <- list(biomass, hsmax, nitrate, temperature)

data_all <- c()

# First Area 
data_extracted <- area[, year >= 1984] 
# Create a raster brick of quarter-year kelp area 
c <- rasterize(x = coords,                 
               y = base_grid,              
               field = data_extracted,    
               fun = F_filter_NAs,
               na.rm = T)

a <- rasterize(x = coords,                 
               y = base_grid,              
               field = data_extracted,    
               fun = sum, # only works for area 
               na.rm = T)

k <- c*a

names(k) <- paste0("Q", study_quarters, "_", study_years)

data_quarter <- rasterToPoints(k) %>% 
  as.data.frame() %>% 
  filter(y <= 37.4) 

data_quarter_long <- data_quarter %>% 
  pivot_longer(Q1_1984:Q4_2021, names_to = "Time", values_to = "value") %>% 
  mutate(variable = "area") %>% 
  arrange(Time)



data_all <- rbind(data_all, data_quarter_long)

# For loop for the rest of the variables 
for (i in 1:length(variables)) {
  
  var <- var_list[[i]]
  
  data_extracted <- var[, year >= 1984] 
  
  # Create a raster brick of quarter-year kelp variable
  k <- rasterize(x = coords,                 
                 y = base_grid,              
                 field = data_extracted,    
                 fun = mean, 
                 na.rm = T)

  names(k) <- paste0("Q", study_quarters, "_", study_years)
  
  data_quarter <- rasterToPoints(k) %>% 
    as.data.frame() %>% 
    filter(y <= 37.4) 
  
  data_quarter_long <- data_quarter %>% 
    pivot_longer(Q1_1984:Q4_2021, names_to = "Time", values_to = "value") %>% 
    mutate(variable = variables[i]) %>% 
    arrange(Time)
  
  data_all <- rbind(data_all, data_quarter_long)
}

df <- data_all %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  left_join(., depth_p, by = c("x", "y")) %>% 
  rename("lat" = "y", 
         "long" = "x") %>% 
  relocate(depth, .after = area) %>% 
  mutate(year = str_sub(.$Time, start= -4),
         quarter = str_sub(.$Time, 1, 2)) %>% 
  dplyr::select(-Time)

##### Export ###################################################################
write.csv(df, "Processed_data/data_tables/kelp_data_per_quarter.csv", row.names = F)

# End of script