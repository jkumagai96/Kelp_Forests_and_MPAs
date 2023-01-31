# Date: January 31st 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Explore kelp per region and mpas per region
# BIO 202: Ecological Statistics

###### Research Question #######################################################
# Do marine protected areas increase the resilience of kelp forests (increase in kelp area)
# following a marine heat wave?

# Is there more kelp area in marine protected areas? Partial vs. fully protected?

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(sf)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")
points_in_mpas <- read.csv("Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv") 

# Mpas 
mpas_original_south <- read_sf("Data/Filtered_MPAs/MPAs_CA_Giant_Kelp.shp") %>% 
  st_transform(crs = 4326) %>% 
  dplyr::select(Site_ID_12, Estab_Yr_1, AreaMar_12, Status_12)

mpas_original_north <- read_sf("Data/Filtered_MPAs/Nothern_California.shp") %>% 
  dplyr::select(Site_ID, Estab_Yr, AreaMar, Status) %>% 
  rename("Site_ID_12" = "Site_ID",
         "Estab_Yr_1" = "Estab_Yr",
         "AreaMar_12" = "AreaMar",
         "Status_12" = "Status")

##### Explore Data #############################################################
# Count and visualize the number of pixels in each region 
kelp_data_all %>% 
  ggplot(aes(x = lat)) +
  geom_histogram(color = "black", fill = "white") +
  geom_vline(aes(xintercept = 34.4486), color = "blue", size = .5) + 
  geom_vline(aes(xintercept = 37.1819), color = "blue", size = .5) +
  geom_vline(aes(xintercept = 39.0044), color = "blue", size = .5) +
  theme_bw()

kelp_data_all %>% count(region) %>% arrange(n)

# Count and visualize the number of pixels that fall within MPAs in each region 
kelp_data_all %>% 
  filter(mpa_status != "None") %>% 
  ggplot(aes(x = lat)) +
  geom_histogram(color = "black", fill = "white") +
  geom_vline(aes(xintercept = 34.4486), color = "blue", size = .5) + 
  geom_vline(aes(xintercept = 37.1819), color = "blue", size = .5) +
  geom_vline(aes(xintercept = 39.0044), color = "blue", size = .5)

kelp_data_all %>% 
  filter(mpa_status != "None") %>% 
  count(region)

# Combine north and south mpas
mpas <- rbind(mpas_original_north, mpas_original_south) %>% 
  rename("mpa_status" = "Status_12") %>% 
  st_transform(crs = 4326) 

### Assign regions to each MPA

# Note: MPA CA265 is in between the central coast and the south coast, 
# but it is mostly in the south coast, so it has been categorized as that. 
# The line which removes number 75 is where I overwrite the region. 
df <- data.frame(st_coordinates(mpas)) %>%
  group_by(L3) %>% 
  mutate(region = ifelse(Y > 34.4486, "Central_Coast", "South_Coast")) %>% 
  mutate(region = ifelse(Y > 37.1819, "North_Central_Coast", region)) %>% 
  mutate(region = ifelse(Y> 39.0044, "North_Coast", region)) %>% 
  distinct(L3, region)

df <- df[-75,]

mpas_w_regions <- mpas %>% 
  mutate(L3 = 1:nrow(mpas)) %>% 
  left_join(df, by = "L3")

plot(mpas_w_regions[,7])

### Count the number of MPAs per region and the total area
mpas_w_regions %>% 
  st_drop_geometry() %>% 
  count(region)
# All regions have similar numbers of MPAs, except for the south which has a lot more

mpas_w_regions %>% 
  st_drop_geometry() %>% 
  group_by(region) %>% 
  summarize(total_area = sum(AreaMar_12)) %>% 
  arrange(total_area)
# South coast has the most area, followed by central coast, north central coast, 
# and then north coast 
