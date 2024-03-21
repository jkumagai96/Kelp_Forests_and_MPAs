# Date: October 25th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Exploring PISCO dataset to see if we can combine the dataset with Reefcheck
# BIO 202: Ecological Statistics

# Load packates 
library(readxl)
library(tidyverse)

# Load Data 
df <- read.csv("Data/PISCO subtidal data clearinghouse/PISCO_kelpforest_combined_data_bysite.csv")
PISCO_master_species <- read_excel("Data/PISCO subtidal data clearinghouse/master_spp_table.xlsx")
PISCO_swath <- read.csv("Data/PISCO subtidal data clearinghouse/PISCO_SWATH.csv")
PISCO_fish <- read.csv("Data/PISCO subtidal data clearinghouse/PISCO_FISH.csv")

# Format Data
invert_species <- PISCO_master_species %>% 
  filter(ScientificName == "Mesocentrotus franciscanus" |      # Red urchin
           ScientificName == "Strongylocentrotus purpuratus" | # Purple urchin
           ScientificName == "Centrostephanus coronatus" |     # crowned urchin
           ScientificName == "Panulirus interruptus") %>%      # Spiny Lobster 
  filter(sample_type == "SWATH")  


test <- invert_species %>% 
  filter(pisco_classcode == "CENCOR" |
           pisco_classcode == "MESFRAAD" |
           pisco_classcode == "STRPURAD" |
           pisco_classcode == "PANINT") %>% 
  filter(campus == "UCSB" |
           campus == "UCSC")


fish_species <- PISCO_master_species %>% 
  filter(ScientificName == "Semicossyphus pulcher")            # Sheephead

species_oi <- rbind(invert_species, fish_species)
unique(species_oi$pisco_classcode) # MESFRAREC and STRPURREC are recruits 

# Data of interest 
df %>% 
  dplyr::select(1:9, fish_SPUL, swath_PANINT, swath_STRPURAD, swath_MESFRAAD, swath_CENCOR)

# Do they survey each site multiple times a year? Yes, sometimes the sites take 
# multiple days of field work 
sites <- unique(PISCO_swath$site)
for (i in 1:length(sites)) {
  site_data <- PISCO_swath %>% 
    filter(site == sites[i]) %>% 
    distinct(year, month, day)
  
  n_records <- nrow(site_data)
  n_years <- length(unique(site_data$year))
  
  n_records == n_years
  
  if (n_records != n_years) {
    print(paste(i, sites[i]))
    print(site_data)
    
  }
  
}

# How did they summarize this data? I can't seem to replicate it
df %>% 
  dplyr::select(1:9, fish_SPUL, swath_PANINT, swath_STRPURAD, swath_MESFRAAD, swath_CENCOR) %>% 
  filter(SITE == "STEWARTS_POINT_MPA_2")

PISCO_swath %>% 
  filter(site == "STEWARTS_POINT_MPA_2") %>% 
  filter(classcode == "STRPURAD" ) %>% 
  group_by(year, site, classcode) %>% 
  summarize(mean_pu = mean(count, na.rm = T), 
            sum_pu = sum(count, na.rm = T))

# How many campuses contributed to this dataset?
unique(PISCO_swath$campus)
unique(PISCO_fish$campus)

# Do multiple campuses survey the same areas? 
for (i in 1:length(sites)) {
  site_data <- PISCO_swath %>% 
    filter(site == sites[i]) %>% 
    distinct(campus)
  
  n_records <- nrow(site_data)
  
  if (n_records > 1) {
    print(paste(i, sites[i]))
    print(site_data)
  }
  
}

# Do they include zeros in the rawish data? No they only include it when there are records 
test_pu <- PISCO_swath %>% 
  filter(site == "ANACAPA_ADMIRALS_E") %>% 
  filter(classcode == "STRPURAD")

test_sh <- PISCO_swath %>% 
  filter(site == "ANACAPA_ADMIRALS_E") %>% 
  filter(classcode == "PANINT")

# How many records include size data for California Sheephead? 
# 99.9% of the data has size records   
PISCO_fish %>% 
  filter(classcode == "SPUL") %>% 
  nrow()

PISCO_fish %>% 
  filter(classcode == "SPUL") %>% 
  filter(fish_tl > 0) %>% 
  nrow()

# which species were looked for?

unique(PISCO_fish$campus)
unique(PISCO_swath$campus)
