# Date: November 22nd 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Combine SST NetCDFs
# BIO 202: Ecological Statistics

##### Set up ###################################################################

# Load Packages 

library(tidyverse)
library(rerddap)
library(doParallel)
library(tidync) # For easily dealing with NetCDF data

# Set up parallel processing 
doParallel::registerDoParallel(5)

#### Load Data ################################################################

SST_load <- function(file_name){
  require(tidyverse)
  OISST_dat <- tidync(file_name) %>%
    hyper_tibble() %>% 
    dplyr::mutate(time_full = as.POSIXct(time, origin = "1970-01-01")) %>% 
    dplyr::mutate(time_final = format(time_full, format="%Y-%m-%d")) %>% 
    dplyr::select(longitude, latitude, time_final, sst) %>% 
    dplyr::rename(t = time_final, temp = sst) %>% 
    rename(lon = longitude, lat = latitude) %>% 
    dplyr::mutate(t = as.Date(t)) %>% 
    na.omit()
  return(OISST_dat)
}

# Locate the files that will be loaded
SST_files <- dir("Data/SST", full.names = T)

# Load the data in parallel
SST_dat <- plyr::ldply(.data = SST_files, .fun = SST_load, .parallel = T)

dir.create("Processed_data/SST")

# Export 
saveRDS(object = SST_dat,
        file = "Processed_data/SST/SST_1983_2021.rds")
