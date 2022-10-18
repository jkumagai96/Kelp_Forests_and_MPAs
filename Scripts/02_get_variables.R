# Date: October 13th 2022
# Author: Juan Carlos' Villaseñor-Derbez (juancvd@stanford.edu) and Joy Kumagai (kumagaij@stanford.edu) 
# Author 2: This code was adapted by Joy, but originally written by Juan Carlos
# Purpose: Create raster's of all data w/ yearly summary's 
# BIO 202: Ecological Statistics

##### Preparing packages and data ##############################################
# Load packages ----------------------------------------------------------------
library(here)
library(ncdf4)
library(raster)
library(stars)
library(tidyverse)

# Load data --------------------------------------------------------------------

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

## Process area and biomass variables ##########################################
# This section of the code is a for loop which goes through a series of steps 
# to extract each variable into a raster format averaged by year. 

# We start with just area, as this is the only variable that is summed to 1km2
# Initialize variables 
v_years <- year[year >= 2009]

# Create a vector of filenames
names_new <-
  paste(
    "LandsatKelp_Quarterly",
    "area",
    unique(v_years),
    sep = "_"
  )

# Create a matrix of years we want, where each column is one quarter of data

data_extracted <- area[, year >= 2009] 

# Create a raster brick of quarter-year kelp area 
k <- rasterize(x = coords,                 
               y = base_grid,              
               field = data_extracted,    
               fun = sum, # only works for area 
               na.rm = T)

means <- stackApply(k, indices=v_years, fun=mean)

# Export one raster per year 
writeRaster(
  x = means,
  bylayer = T,
  format = "GTiff", 
  filename = here::here(
    "Processed_data",
    "kelp_variables", 
    paste0(names_new, ".tif")
  ), 
  overwrite = T
)

# Each loop of the for loop is a different variable. 
# First, we create file names to export the final rasters, then we extract the 
# relevant data within the variable into a matrix, we then create a raster brick 
# of all the data. This raster brick is then averaged by year and exported as a 
# TIFF file.  

## Process biomass, nitrate, hsmax, and temperature  variables #################
# biomass,nitrate, hsmax, and temperature variables are averaged! 
# Initialize variables 
variables <- list(biomass, temperature,  nitrate, hsmax)
var_names <- c("biomass","temperature", "nitrate", "hsmax")

# start for loop
for (i in 1:length(variables)) {
  
  # Create a vector of filenames
  names_new <-
    paste(
      "LandsatKelp_Quarterly",
      var_names[i],
      unique(v_years),
      sep = "_"
    )
  
  # Create a matrix of years we want, where each column is one quarter of data
  
  data_extracted <- variables[[i]][, year >= 2009] 
  
  # Create a raster brick of quarter-year kelp area 
  k <- rasterize(x = coords,                 
                 y = base_grid,              
                 field = data_extracted,    
                 fun = mean, # IMPORTANT LINE OF CODE  
                 na.rm = T)
  
  means <- stackApply(k, indices=v_years, fun=mean)
  
  # Export one raster per year 
  writeRaster(
    x = means,
    bylayer = T,
    format = "GTiff", 
    filename = here::here(
      "Processed_data",
      "kelp_variables", 
      paste0(names_new, ".tif")
    ), 
    overwrite = T
  )
  
}

## Create depth raster #########################################################
# Create raster with depth data 
depth_r <- rasterize(x = coords,                 
               y = base_grid,              
               field = depth,    
               fun = mean,
               na.rm = T)

writeRaster(
  x = depth_r,
  format = "GTiff", 
  filename = here::here(
    "Processed_data",
    "kelp_variables", 
    "LandsatKelp_Quarterly_depth.tif"),
  overwrite = T
  )
### End of script ###
