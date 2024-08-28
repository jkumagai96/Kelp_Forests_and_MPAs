# Date:  August 25th 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Clean MLPA data without site_blocks
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(readxl)
library(sf)

# Load MLPA Data
MLPA_master_sites <- read.csv("Data/resourceMap_MLPA_kelpforest_9/data/MLPA_kelpforest_site_table.6.csv")

# Load MLPA duplicated sites that are now corrected
corrected_sites <- read_xlsx("Processed_data/Corrected_MLPA_PISCO_sites.xlsx", sheet = 2) %>% 
  select(-notes)

mpas <- st_read("Processed_data/MPAs.shp")

###### Clean data ##############################################################
MLPA_sites <- MLPA_master_sites

# Remove duplicated rows 
MLPA_sites[duplicated(MLPA_sites),]
MLPA_sites <- distinct(MLPA_sites)

# Remove sites 
## Removed begg rock sites because they are effectively no take 
## Removed Natural bridges and Cabrillo MPA sites because they are not consistently
## sampled within the MPA. Shallow sites are in the MPA and deep sites are outside 
## of it.
Natural_bridges_mpa_sites <- MLPA_sites %>% 
  filter(CA_MPA_Name_Short == "Natural Bridges SMR", 
         site_status == "mpa")

Natural_bridges_mpa_sites <- unique(Natural_bridges_mpa_sites$site)

MLPA_sites <- MLPA_sites %>% 
  filter(!CA_MPA_Name_Short == "Begg Rock SMR") %>% 
  filter(!site == "CABRILLO_NATIONAL_MONUMENT") %>% 
  filter(! site %in% Natural_bridges_mpa_sites) 

# Reduce data 
MLPA_sites_reduced <- MLPA_sites %>% 
  dplyr::select(site, latitude, longitude, CA_MPA_Name_Short) %>% 
  distinct()

# Fix duplication in MLPA sites 
count_of_sites <- MLPA_sites_reduced %>% 
  count(site)

MLPA_sites_reduced <- MLPA_sites_reduced %>% 
  full_join(count_of_sites) 


sites1 <- MLPA_sites_reduced %>% 
  filter(n == 1)

sites2 <- MLPA_sites_reduced %>% 
  filter(n > 1) %>% 
  select(-c(latitude, longitude)) %>% 
  distinct() %>% 
  full_join(corrected_sites, by = "site")

MLPA_sites_corrected <- rbind(sites1, sites2) %>% 
  select(-c(n, CA_MPA_Name_Short))

##### Spatiall join MPAs #######################################################
site_points <- MLPA_sites_corrected %>% 
  mutate(region = ifelse(latitude > 34.4486, "Central_Coast", "South_Coast")) %>% 
  mutate(region = ifelse(latitude > 37.1819, "North_Central_Coast", region)) %>% 
  mutate(region = ifelse(latitude > 39.0044, "North_Coast", region)) %>% 
  filter(region == "South_Coast" | region == "Central_Coast") %>% 
  st_as_sf(., coords = c(x = "longitude", y = "latitude"), crs = 4326)  

# Intersect site data with mpas 
sites_in_mpas <- st_intersection(site_points, mpas) %>% 
  st_drop_geometry() %>% 
  select(-region)

site_points <- site_points %>% 
  left_join(sites_in_mpas, by = c("site")) 

site_points$mpa_status[is.na(site_points$mpa_status)] <- "Reference"

sites_for_joining <- st_drop_geometry(site_points)

##### Export MLPA sites for joining ############################################
st_write(site_points, "Processed_data/MLPAs_sites_with_MPAs.shp", append = FALSE)
write.csv(sites_for_joining, "Processed_data/data_tables/MLPA_sites_with_MPAs.csv", row.names = FALSE)

