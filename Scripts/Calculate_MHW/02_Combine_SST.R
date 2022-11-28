# Date: November 22nd 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Combine SST NetCDFs
# BIO 202: Ecological Statistics

##### Set up ###################################################################
# Load Packages 
library(tidyverse)
library(tidync) # For easily dealing with NetCDF data
library(rerddap) # For easily downloading subsets of data
library(doParallel) # For parallel processing

# Set up parallel processing 
detectCores()
doParallel::registerDoParallel(3)

#### Load Data ################################################################

SST_load <- function(file_name){
  OISST_dat <- tidync(file_name) %>%
    hyper_tibble() %>% 
    mutate(time_full = as.POSIXct(time, origin = "1970-01-01")) %>% 
    mutate(time_final = format(time_full, format="%Y-%m-%d")) %>% 
    select(longitude, latitude, time_final, sst) %>% 
    dplyr::rename(t = time_final, temp = sst) %>% 
    rename(lon = longitude, lat = latitude) %>% 
    mutate(t = as.Date(t)) %>% 
    na.omit()
  return(OISST_dat)
}

# Locate the files that will be loaded
SST_files <- dir("Data/SST", full.names = T)

# Load the data in parallel
SST_dat <- plyr::ldply(.data = SST_files, .fun = SST_load, .parallel = T)

# Visualize data
plot1 <- SST_dat %>% 
  filter(t == "2019-06-30") %>% 
  ggplot(aes(x = lon, y = lat)) +
  geom_tile(aes(fill = temp)) +
  scale_fill_viridis_c() +
  coord_quickmap(expand = F) +
  labs(x = NULL, y = NULL, fill = "SST (°C)") +
  theme(legend.position = "bottom")

ggsave("Figures/SST_20190630.png")

##### Export ###################################################################
saveRDS(object = SST_dat,
        file = "Processed_data/SST/SST_1984_2021.rds")
