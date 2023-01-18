# Date: November 23rd 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Download daily SST in a specific bounding box (0.25 degree)
# BIO 202: Ecological Statistics

##### Set up ###################################################################
# Load Packages
library(tidyverse)

##### Process ##################################################################
# Declare variables
year <- 1984:2021

# Download data
for (i in 1:length(year)) {
  website <- paste0("https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21Agg_LonPM180.nc?sst%5B(", year[i], "-01-01T12:00:00Z):1:(", year[i], "-12-31T12:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(42):1:(32.5)%5D%5B(-125):1:(-117)%5D")
  file_name <- paste0("Data/SST/SST_", year[i], ".nc")
  download.file(website, destfile = file_name, quiet = TRUE)
  print(i)
}

# Url created from https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21Agg_LonPM180.html
