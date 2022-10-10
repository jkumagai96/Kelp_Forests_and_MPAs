# Date: October 9th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu)
# Purpose: Exploring new data and creating a table from NetCDF w/ enviornment data 
# BIO 202: Ecological Statistics

##### Load Packages and Data #####
library(tidyverse)
library(raster)
library(ncdf4)
# exact_extractR

kelp_file <- "Data/Kelp/CAKelpCanopyEnv_2021.nc"
kelp_all <- nc_open(kelp_file)

##### Explore Data #####
kelp_all

# Station variables
lon <- ncvar_get(kelp_all,"lon") 
lat <- ncvar_get(kelp_all, "lat")
depth <- ncvar_get(kelp_all, "station") # does not work!!!! COME BACK HERE

# Time variables 
year <- ncvar_get(kelp_all, "year")
quarter <- ncvar_get(kelp_all, "quarter")

# Station and time variables
biomass <- ncvar_get(kelp_all, "biomass")
temperature <- ncvar_get(kelp_all, "temperature")
area <- ncvar_get(kelp_all, "area")
hsmax <- ncvar_get(kelp_all, "hsmax")
nitrate <- ncvar_get(kelp_all, "nitrate")

dim(biomass) # 152 columns for time and 332,640 rows for each station 

##### Format Data #####
# Create reference tables w/ IDs
station_data <- data.frame(lat, lon, StationID = 1:dim(lat)) %>% 
  filter(lat <= 37.3 & lat >= 32.54) # Just giant kelp in CA 

time_data <- data.frame(year, quarter, TimeID = 1:length(year)) 

# Create usable tables of biomass, temp, hsmax, and nitrate data 
biomass_data <- reshape2::melt(as.matrix(biomass)) 
colnames(biomass_data) <- c("StationID", "TimeID", "biomass")

#biomass_data %>% 
#  left_join(time_data, by = "TimeID") %>% 
#  group_by(StationID, year) %>% 
#  summarise(biomass_mean = mean(biomass, na.rm = T))

temperature_data <- reshape2::melt(as.matrix(temperature))
colnames(temperature_data) <- c("StationID", "TimeID", "temperature")

hsmax_data <- reshape2::melt(as.matrix(hsmax))
colnames(hsmax_data) <- c("StationID", "TimeID", "hsmax")

nitrate_data <- reshape2::melt(as.matrix(nitrate))
colnames(nitrate_data) <- c("StationID", "TimeID", "nitrate")

kelp_data <- cbind(biomass_data, 
                   temperature_data$temperature,
                   hsmax_data$hsmax,
                   nitrate_data$nitrate) %>% 
  left_join(time_data, by = "TimeID") %>% 
  rename(temperature = `temperature_data$temperature`,
         hsmax = `hsmax_data$hsmax`,
         nitrate = `nitrate_data$nitrate`)

#### Format: Summarize data by year ####
start <- Sys.time()
kelp_data_yr <- kelp_data %>%                  
  filter(StationID >= min(station_data$StationID) &  # Filter to CA boundaries of giant kelp
           StationID <= max(station_data$StationID), # Filter to CA boundaries of giant kelp
         year >= 2003) %>%                           # Filter to only look at years of 2003 to present
  arrange(StationID, year) %>% 
  group_by(year, StationID) %>% 
  summarize(biomass_mean = mean(biomass, na.rm = T),
            SST_mean = mean(temperature, na.rm = T),
            hsmax_mean = mean(hsmax, na.rm = T),
            nitrate_mean = mean(nitrate, na.rm = T)) %>% 
  left_join(station_data, by = "StationID")         # Adds in Station location ifo
end <- Sys.time()

end - start

write.csv(kelp_data_yr, "Processed_data/kelp_average_2003_2022.csv") 

##### Create Station Point Dataset w/ MPA #####
# create raster with ID, extract mpa coverage per pixel (ID), 
library(sf) 
station_all <- data.frame(lat, lon, StationID = 1:dim(lat))

station_points = st_as_sf(t, coords = c("lon", "lat"), 
                 crs = 4326, agr = "constant")

# Load in MPA data and filter it 
mpas_all <- st_read("Data/Marine Protected Areas/NOAA_MPAI_2020_IUCN_gdb/NOAA_MPAI_v2020.gdb")
mpas_CA <- mpas_all %>% 
  filter(State == "CA" & Estab_Yr <= 2014) %>% # 2014 is the year the marine heat wave started
  st_transform(crs = 4326)

# export point and mpa data
st_write(station_points, "Processed_data/station_points.shp")
st_write(mpas_CA, "Processed_data/mpas_CA.shp")

# Intersect points with MPAs in CA
mpas_CA <- st_make_valid(mpas_CA)
intersection_list <- st_intersects(station_points, mpas_CA)
points_in_mpas <- lengths(intersection_list) > 0

station_points$MPA <- points_in_mpas 

# There are multiple MPAs per point

##### Join Data and Export #####
kelp_data_yr %>% 
  left_join(st_drop_geometry(station_points), by = "StationID") %>% 
  write.csv("Processed_data/kelp_environ_mpas_2003_2022.csv", row.names = F)

# Goal is to export a table with lat, long, biomass, mpa presence, SST, wave exposure, depth?
# Steps: Reduce years and lat/long, create table w/ ID, biomass, SSt, wave exposure, and depth
# done

# create dataframe with ID, extract mpa presence/absence
# Join the coverage data frame with the dataframe w/ ID, biomass, SST, wave exposure, and depth

