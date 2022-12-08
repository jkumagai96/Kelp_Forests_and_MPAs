# Date: November 28th 2022
# Author: Joy A. Kumagai (kumagaij@stanford.edu)
# Purpose: Create a data table with additional variables including, MPA 
#          protection, human gravity index, and marine heat waves and cold spells 
#          per point with kelp data.
# BIO 202: Ecological Statistics

##### Set up: packages #########################################################
# load packages 
library(sf)
library(tidyverse)

##### Set up: load data ########################################################
mpas_original <- read_sf("Data/Filtered_MPAs/MPAs_CA_Giant_Kelp.shp")
station_data <- read.csv("Processed_data/data_tables/PixelID_reference.csv")
kelp_data <- read.csv("Processed_data/data_tables/kelp_data_per_quarter.csv")
human_gravity <- read.csv("Data/Population/human_gravity_for_kelp_patches.csv") %>% 
  rename(long = lon)

# Marine heatwaves and cold spells 
cold_spells <- readRDS("Processed_data/SST/CS_cummulative_intensity_1km.rds") %>% 
  filter(lat <= 37.4) %>% 
  filter(long >= -122.5)

heat_waves <- readRDS("Processed_data/SST/MHW_cummulative_intensity_1km.rds") %>% 
  filter(lat <= 37.4) %>% 
  filter(long >= -122.5)

##### Formatting ###############################################################
# Select needed attribuets from MPAs
mpas <- mpas_original %>% 
  dplyr::select(Site_ID_12, Estab_Yr_1, AreaMar_12, Status_12)  %>% # I selected Area_Mar12, and not Area? QUESTION
  rename("mpa_status" = "Status_12") %>% 
  # Match projections with kelp data 
  st_transform(crs = 4326) 

# convert station points into spatial data to intersect with mpas 
station_points <- st_as_sf(station_data, 
                           coords = c("x", "y"), 
                           crs = 4326) 

# Make sure the projections match 
raster::crs(station_points) == raster::crs(mpas)


##### Processing ###############################################################
# Spatial intersect between mpas and points to get all points within mpas 
points_in_mpas <- st_intersection(station_points, mpas) %>% 
  dplyr::select(PixelID, mpa_status, Site_ID_12, Estab_Yr_1, AreaMar_12) %>% 
  st_drop_geometry()

# Join data with kelp data 
kelp_w_mpas <- left_join(kelp_data, points_in_mpas, by = "PixelID") 

# Format values of mpa_status so NA is none, no take is full, and uniform multiple use is partial
kelp_w_mpas$mpa_status[is.na(kelp_w_mpas$mpa_status)] <- "None"

not_protected <-kelp_w_mpas %>% 
  filter(mpa_status == "None") %>% 
  arrange(PixelID)

subset <- kelp_w_mpas %>% 
  filter(mpa_status != "None")

for (i in 1:nrow(subset)) {
  status <- subset$mpa_status[i]
  establishment_yr <- subset$Estab_Yr_1[i]
  
  if (subset$year[i] >= establishment_yr) {
    subset$mpa_status[i] <- status
  } else {
    subset$mpa_status[i] <- "None"
    subset$AreaMar_12[i] <- 0
  }
}

subset <- subset %>% arrange(PixelID)

kelp_w_mpas_adjusted <- rbind(not_protected, subset)

# Formatting final data 
final_data <- kelp_w_mpas_adjusted %>% 
  rename("Mpa_ID" = "Site_ID_12") %>% 
  rename("mpa_area" = "AreaMar_12") %>% 
  relocate(Mpa_ID, .after = depth) %>% 
  relocate(mpa_status, .after = Mpa_ID) %>% 
  relocate(mpa_area, .after = mpa_status) %>% 
  select(-Estab_Yr_1)

##### Add in Human Gravity Index ###############################################

final_data <- final_data %>% 
  left_join(human_gravity, by = c("long", "lat")) %>%  # NA's are values > 50km
  mutate(gravity = replace_na(gravity, 0)) # 11% of the data are zero's now 

##### Add in marine heat wave and cold spells ##################################

# Format marine heat wave and cold spell data 
temp_anon <- cbind(heat_waves, cold_spells$CS_cummulative) %>% 
  rename(CS_intensity = "cold_spells$CS_cummulative",
         MHW_intensity = "MHW_cummulative") %>%  
  mutate(year = as.integer(year))

# Join temperature anonmoly data by long, lat and year
final_data <- final_data %>% 
  left_join(temp_anon, by = c("long", "lat", "year"))

# How many cells did not join 
n_na <- final_data %>% filter(is.na(MHW_intensity)) %>% nrow()
total <- nrow(final_data)

n_na/total*100 # 10.8% of the data are NAs now 

##### Export ###################################################################

write.csv(final_data, "Processed_data/data_tables/kelp_data_w_mpas_per_quarter.csv", 
          row.names = F)

# End of script 