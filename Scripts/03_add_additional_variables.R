# Date: December 12th, 2022
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
mpas_original_south <- read_sf("Data/Filtered_MPAs/MPAs_CA_Giant_Kelp.shp") %>% 
  st_transform(crs = 4326) %>% 
  dplyr::select(Site_ID_12, Site_Name_, Estab_Yr_1, AreaMar_12, Status_12) %>% 
  rename("Site_Name" = "Site_Name_")

mpas_original_north <- read_sf("Data/Filtered_MPAs/Nothern_California.shp") %>% 
  dplyr::select(Site_ID, Site_Name, Estab_Yr, AreaMar, Status) %>% 
  rename("Site_ID_12" = "Site_ID",
       "Estab_Yr_1" = "Estab_Yr",
       "AreaMar_12" = "AreaMar",
       "Status_12" = "Status")

# Combine north and south mpas
mpas_original <- rbind(mpas_original_north, mpas_original_south)

station_data <- read.csv("Processed_data/data_tables/PixelID_reference.csv")
kelp_data <- read.csv("Processed_data/data_tables/kelp_data_per_quarter.csv")
human_gravity <- read.csv("Data/Population/human_gravity_for_kelp_patches.csv") %>% 
  rename(long = lon)

# Marine heatwaves and cold spells 
cold_spells <- readRDS("Processed_data/SST/CS_cummulative_intensity_1km.rds") 

heat_waves <- readRDS("Processed_data/SST/MHW_cummulative_intensity_1km.rds") 

# Count of 30x30 pixels within each 1km2 pixel
count_30x30 <- read.csv("Processed_data/data_tables/Count_30x30_pixels.csv")

##### Formatting ###############################################################
# Select needed attributes from MPAs
mpas <- mpas_original %>% 
  rename("mpa_status" = "Status_12") %>% 
  # Match projections with kelp data 
  st_transform(crs = 4326) 

# convert station points into spatial data to intersect with mpas 
station_points <- st_as_sf(station_data, 
                           coords = c("x", "y"), 
                           crs = 4326) 

# Make sure the projections match 
raster::crs(station_points) == raster::crs(mpas)

# Adjust Abalone Cove and Farnsworth MPAs
mpas <- mpas %>% 
  mutate(mpa_status = case_when(Site_Name == "Abalone Cove State Marine Conservation Area" ~ "Full", 
                                Site_Name == "Farnsworth Onshore (Catalina Island) State Marine Conservation Area" ~ "Full",
                                .default = mpa_status)
  ) %>% select(-Site_Name)

# Export MPA data
st_write(mpas, "Processed_data/MPAs.shp", append = FALSE)

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
#  final_data <- kelp_w_mpas %>% 
  arrange(PixelID) %>% 
  rename("Mpa_ID" = "Site_ID_12") %>% 
  rename("mpa_area" = "AreaMar_12") %>% 
  relocate(Mpa_ID, .after = depth) %>% 
  relocate(mpa_status, .after = Mpa_ID) %>% 
  relocate(mpa_area, .after = mpa_status) %>% 
  dplyr::select(-Estab_Yr_1)

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

# Join temperature anomaly data by long, lat and year
data_w_temp <- final_data %>% 
  left_join(temp_anon, by = c("long", "lat", "year")) %>% 
  mutate(has_temp_data = !is.na(MHW_intensity)) 

# How many cells did not join 
n_na <- data_w_temp %>% filter(is.na(MHW_intensity)) %>% nrow()
total <- nrow(data_w_temp)

n_na/total*100 # 11.5% of the data are NAs now 

### Find the nearest heatwave/cold spell data for unique pixels
temp_points <- data_w_temp %>% 
  dplyr::select(lat, long, PixelID, has_temp_data) %>% 
  unique() %>% 
  st_as_sf(coords = c("long", "lat"), remove = FALSE)

na_temp <- temp_points %>% 
  filter(has_temp_data == FALSE)

has_temp <- temp_points %>% 
  filter(has_temp_data == TRUE)

joined_points <- st_join(na_temp, has_temp, join = st_nearest_feature) %>% 
  st_drop_geometry() %>% 
  dplyr::select(PixelID.x, PixelID.y) 

# Assign the blank MHW and coldspell values to the nearest pixel with data 
# Logic is to go through each pixel within the final data that is PixelID.x, then 
# find the corresponding PixelID.y MHW value each year, then assign it to the same one
lookup_table <- data_w_temp %>% 
  dplyr::select(PixelID, year, quarter, MHW_intensity, CS_intensity)

fixed_data <- c()

for (i in joined_points$PixelID.x) {
  
  nearest_pixel <- joined_points[joined_points$PixelID.x == i, 2] 
  
  data_to_be_joined <- lookup_table %>% 
    filter(PixelID == nearest_pixel)
  
  fixed_data_per_pixel <- data_w_temp %>% 
    filter(PixelID == i) %>% 
    dplyr::select(-c(MHW_intensity, CS_intensity)) %>% 
    full_join(data_to_be_joined, by = c("year", "quarter")) 
  
  fixed_data <- rbind(fixed_data, fixed_data_per_pixel)
  
}
# clean up the column names 
fixed_data <- fixed_data %>% 
  dplyr::select(-PixelID.y) %>% 
  rename(PixelID = PixelID.x) 

# Combine the original data with the fixed data. 
a <- data_w_temp %>% 
  filter(has_temp_data == TRUE)

final_data <- rbind(a, fixed_data) %>% 
  arrange(PixelID) %>% 
  dplyr::select(-has_temp_data)

##### Add in regions according to Marine Life Protecion Act Regions ##########
final_data <- final_data %>% 
  mutate(region = ifelse(lat > 34.4486, "Central_Coast", "South_Coast")) %>% 
  mutate(region = ifelse(lat > 37.1819, "North_Central_Coast", region)) %>% 
  mutate(region = ifelse(lat > 39.0044, "North_Coast", region))

##### Remove pixels with 5 or less 30x30 pixels and other regions ##############
final_data <- final_data %>% 
  left_join(count_30x30, by = "PixelID") %>% 
  filter(count > 5) %>% 
  filter(region == "South_Coast" | region == "Central_Coast") # Filter by regions 

region_data <- final_data %>% 
  dplyr::select(PixelID, region) %>% 
  unique()

# Filter poitns in mpas by regions
points_in_mpas <- points_in_mpas %>% 
  left_join(region_data, by = "PixelID") %>% 
  na.omit(region)

##### Export ###################################################################
# Export points_in_mpas so it can be used in the permutation analysis 
write.csv(points_in_mpas, 
          "Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv",
          row.names = F)

# Export final data 
write.csv(final_data, "Processed_data/data_tables/kelp_data_w_mpas_per_quarter.csv", 
          row.names = F)

# End of script 
