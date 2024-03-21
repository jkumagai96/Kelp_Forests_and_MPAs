# Date: February 5th 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Summarize PISCO raw data to site level data so that it can be analyzed
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(readxl)
library(sf)

# Load Data
PISCO_fish_raw <- read.csv("Data/PISCO subtidal data clearinghouse/PISCO_FISH.csv")
PISCO_swath <- read.csv("Data/PISCO subtidal data clearinghouse/PISCO_SWATH.csv")
PISCO_master_species <- read_excel("Data/PISCO subtidal data clearinghouse/master_spp_table.xlsx")
PISCO_master_sites <- read_excel("Data/PISCO subtidal data clearinghouse/PISCO subtidal master site table.xlsx")

mpas <- st_read("Processed_data/MPAs.shp") # mpas

###### Decisions summarized
# Looking at only central and southern california
# We are not including urchin recruits in the analysis
# We assume that any transect which has some data, but not the species of interest
#     the number of species of interest is zero 
# Sites and transects are weighted equally

# We are only looking at data during and after 2002 because the methods are consistent
#     for all species of interest (crowned urchins were counted using a different method)
#     and the large abundances of sheephead (most likely due to the 1998 large El Nino
#     are mostly gone by then 

# Cutoff of 10cm

year_cutoff <- 2002

###### Create PISCO sites dataset with MPA data  ###############################
site_points <- PISCO_master_sites %>% 
  mutate(longitude = as.numeric(LONGITUDE), 
         latitude = as.numeric(LATITUDE)) %>% 
  mutate(region = ifelse(latitude > 34.4486, "Central_Coast", "South_Coast")) %>% 
  mutate(region = ifelse(latitude > 37.1819, "North_Central_Coast", region)) %>% 
  mutate(region = ifelse(latitude > 39.0044, "North_Coast", region)) %>% 
  filter(region == "South_Coast" | region == "Central_Coast") %>% 
  drop_na(longitude) %>% 
  st_as_sf(., coords = c(x = "longitude", y = "latitude"), crs = 4326)  

# Intersect site data with mpas 
sites_in_mpas <- st_intersection(site_points, mpas) %>% 
  st_drop_geometry() %>% 
  dplyr::select(SITE, LATITUDE, LONGITUDE, Site_ID_12, Estab_Yr_1, AreaMar_12, mpa_status)

site_points <- site_points %>% 
  dplyr::select(-REGION) %>% 
  left_join(sites_in_mpas, by = c("SITE", "LATITUDE", "LONGITUDE"))

site_points$mpa_status[is.na(site_points$mpa_status)] <- "Reference"

sites_for_joining_spatial <- site_points %>% 
  dplyr::select(SITE, site_status, Site_ID_12, Estab_Yr_1, AreaMar_12, mpa_status, region)

sites_for_joining <- st_drop_geometry(sites_for_joining_spatial)

# Export PISCO sites for joining 
st_write(sites_for_joining_spatial, "Processed_data/PISCO_sites_with_MPAs.shp", append = FALSE)
write.csv(sites_for_joining, "Processed_data/data_tables/PISCO_sites_with_MPAs.csv", row.names = FALSE)

##### Create unique transects ##################################################
# Count the number of distinct transects 
distinct_transects_fish <- PISCO_fish_raw %>% 
  filter(year >= year_cutoff,
         level == "BOT") %>% 
  distinct(campus, method, year, month, day, site, zone, transect) %>% 
  left_join(sites_for_joining, by = c("site" = "SITE")) %>% 
  filter(region == "South_Coast" |
           region == "Central_Coast")

distinct_transects_inverts <- PISCO_swath %>% 
  filter(year >= year_cutoff) %>% 
  distinct(campus, method, year, month, day, site, zone, transect) %>% 
  left_join(sites_for_joining, by = c("site" = "SITE")) %>% 
  filter(region == "South_Coast" |
           region == "Central_Coast")

##### Look up unique species codes #############################################
test <- PISCO_master_species %>% 
  filter(ScientificName == "Semicossyphus pulcher" |         # California sheephead
         ScientificName == "Mesocentrotus franciscanus" |    # Red urchin
         ScientificName == "Strongylocentrotus purpuratus" | # Purple urchin
         ScientificName == "Centrostephanus coronatus" |     # crowned urchin
         ScientificName == "Panulirus interruptus") %>%      # Spiny lobster
  filter(sample_type == "SWATH"  |
         sample_type == "FISH")

view(test)
# Species codes:
# CENCOR = Crowned Urchin
# MESFRAAD = Red Urchin
# STRPURAD = Purple Urchin
# PANINT = Spiny lobsters
# SPUL = California sheephead 

##### Format PISCO raw fish data ###############################################
PISCO_fish <- PISCO_fish_raw %>% 
  filter(fish_tl >= 10) %>% # IMPORTANT TO FOLLOW PISCO PROCESSING!!!!! 
  filter(year >= year_cutoff,
         level == "BOT") %>% 
  filter(classcode == "SPUL") %>% 
  select(-c(depth, observer, notes, sex, surge))


fish <- PISCO_fish %>% 
  group_by(campus, method, year, month, day, site, zone, transect) %>% 
  mutate(biomass = (count*((0.0144)*fish_tl^3.04))) %>% 
  # total length is in cm and biomass is in g
  summarize(total_count = sum(count), 
            total_biomass = sum(biomass, na.rm = T))

fish_densities <- distinct_transects_fish %>% 
  left_join(fish, by = c("campus", "method", "year", "month", "day", "site", "zone", "transect")) %>% 
  replace_na(list(total_count = 0,
                  total_biomass = 0)) %>% 
  group_by(year, mpa_status, site, region) %>% 
  summarise(SPUL_d = mean(total_count),
            biomass_d = mean(total_biomass))

##### Format PISCO raw inverts data ############################################
inverts <- PISCO_swath %>% 
  filter(year >= year_cutoff) %>% 
  filter(classcode == "CENCOR" |
           classcode == "MESFRAAD" |
           classcode == "STRPURAD" |
           classcode == "PANINT") %>% 
  pivot_wider(names_from = classcode, values_from = count) %>% 
  select(-c(disease, depth, observer, notes, survey_year)) %>% 
  group_by(campus, method, year, month, day, site, zone, transect) %>% 
  summarize(CENCOR = sum(CENCOR, na.rm = T),
            MESFRAAD = sum(MESFRAAD, na.rm = T),
            STRPURAD = sum(STRPURAD, na.rm = T),
            PANINT = sum(PANINT, na.rm = T))

inverts$total_urchins <- rowSums(inverts[,10:12], na.rm = TRUE, dims = 1)

invert_densities <- distinct_transects_inverts %>% 
  left_join(inverts, by = c("campus", "method", "year", "month", "day", "site", "zone", "transect")) %>% 
  replace_na(list(PANINT = 0,
                  CENCOR = 0,
                  MESFRAAD = 0,
                  STRPURAD = 0, 
                  total_urchins = 0)) %>% 
  group_by(year, mpa_status, site, region) %>% 
  summarise(urchin_d = mean(total_urchins),
            PANINT_d = mean(PANINT), 
            CENCOR_d = mean(CENCOR), 
            MESFRAAD_d = mean(MESFRAAD),
            STRPURAD_d = mean(STRPURAD))
  
# Each individual transect is weighted equally due to the distinct transects

# How many sites per region / mpa_status
PISCO_data_summarized %>% 
  ungroup() %>% 
  distinct(site, region, mpa_status) %>% 
  group_by(region, mpa_status) %>% 
  summarize(n = n()) %>% 
  mutate(100*n/sum(n))

PISCO_data_summarized %>% 
  ungroup() %>% 
  distinct(site, mpa_status) %>% 
  group_by(mpa_status) %>% 
  summarize(n = n()) %>% 
  mutate(100*n/sum(n))

PISCO_data_summarized %>% 
  drop_na(urchin_d) %>% 
  ungroup() %>% 
  distinct(site, region) %>%
  group_by(region) %>% 
  summarize(n = n())

##### Combine data and export ##################################################
PISCO_data_summarized <- full_join(fish_densities, invert_densities) %>% 
  mutate(heatwave = case_when(year < 2014 ~ "before", 
                              year > 2016 ~ "after", 
                              .default = "during"),
         heatwave = factor(heatwave, 
                           levels = c("before", "during", "after")))

write.csv(PISCO_data_summarized, 
          "Processed_data/PISCO_data_summarized.csv",
          row.names = F)
