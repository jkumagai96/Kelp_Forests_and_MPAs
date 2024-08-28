# Date:  August 25th 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Prepare the MLPA data for analysis without site blocks
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)

# Load MLPA Data
MLPA_fish_raw <- read.csv("Data/resourceMap_MLPA_kelpforest_9/data/MLPA_kelpforest_fish.6.csv")
MLPA_swath <- read.csv("Data/resourceMap_MLPA_kelpforest_9/data/MLPA_kelpforest_swath.7.csv")
MLPA_master_species <- read.csv("Data/resourceMap_MLPA_kelpforest_9/data/MLPA_kelpforest_taxon_table.7.csv")
MLPA_sites <- read.csv("Processed_data/data_tables/MLPA_sites_with_MPAs.csv") 

# Outline of steps: 
# A) Look up species codes
# 1) Adjust transect names
# 2) Create unique transects
# 3) Format Fish data
# 4) Format Inverts data
# 5) Combine
# 6) Descriptive statistics
# 7) Export 

###### Decisions summarized ####
# Looking at only central and southern california
# We are not including urchin recruits in the analysis
# We assume that any transect which has some data, but not the species of interest
#     the number of species of interest is zero 
# Sites and transects are weighted equally

# We are only looking at data during and after 2002 because the methods are consistent
#     for all species of interest (crowned urchins were counted using a different method)
#     and the large abundances of sheephead (most likely due to the 1998 large El Nino
#     are mostly gone by then 

# Cutoff of 15cm

year_cutoff <- 2002

##### A) Look up unique species codes #############################################
test <- MLPA_master_species %>% 
  mutate(ScientificName = paste(Genus, Species)) %>% 
  filter(ScientificName == "Semicossyphus pulcher" |         # California sheephead
           ScientificName == "Mesocentrotus franciscanus" |    # Red urchin
           ScientificName == "Strongylocentrotus purpuratus" | # Purple urchin
           ScientificName == "Panulirus interruptus" |
           ScientificName == "Macrocystis pyrifera" |
           ScientificName == "Nereocystis luetkeana") %>%      # Spiny lobster
  filter(sample_type == "SWATH"  |
           sample_type == "FISH")

view(test)
# Species codes:
# MESFRAAD = Red Urchin
# STRPURAD = Purple Urchin
# PANINT = Spiny lobsters
# SPUL = California sheephead 
# MACPYRAD = Giant kelp
# NERLUE = Bull kelp

# UCSB campus apparently counted all obster found on fish transects, not swath transects
# How will we address this?

##### 1) Adjust transect names #####################################################
# Fish
transect_df <- data.frame(transect = unique(MLPA_fish_raw$transect),
                          correct_t = c(3,4,6,5,7,1,2,8,9,10,12,13,14,11,21,22,23,24,25,26,41,51,42,52,43,53,44,54,61,62,63,64,27,28,30,31,32,29,
                                        71,72,73,74,75,76)) 

mlpa_fish <- MLPA_fish_raw %>% 
  full_join(transect_df, by = "transect") %>% 
  dplyr::select(-transect) %>% 
  rename(transect = correct_t)

# SWATH
transect_df <- data.frame(transect = unique(MLPA_swath$transect),
                          correct_t = c(3,4,5,7,1,8,2,6,21,22,41,42,71,72,61,62,73,51,52,43,53)) 

mlpa_swath <- MLPA_swath %>%
  full_join(transect_df, by = "transect") %>% 
  dplyr::select(-transect) %>% 
  rename(transect = correct_t)

##### 2) Create unique transects ##################################################

# Remove sites we don't want in the analysis! If they are not in the master list they are left out of the analysis 
unique_sites <- unique(MLPA_sites$site) 
mlpa_fish <- mlpa_fish %>% 
  filter(site %in% unique_sites)
mlpa_swath <- mlpa_swath %>% 
  filter(site %in% unique_sites)

# Count the number of distinct transects 
distinct_transects_fish <- mlpa_fish %>% 
  filter(year >= year_cutoff,
         level == "BOT") %>% 
  distinct(campus, method, survey_year, year, month, day, site, zone, transect) %>% 
  left_join(MLPA_sites, by = "site")

distinct_transects_inverts <- mlpa_swath %>% 
  filter(year >= year_cutoff) %>% 
  distinct(campus, method, survey_year, year, month, day, site, zone, transect) %>% 
  left_join(MLPA_sites, by = "site")

##### 3) Format Fish Data ######################################################
fish <- mlpa_fish %>% 
  filter(fish_tl >= 15) %>% # Adjusted to 15cm!!!!!! 
  filter(year >= year_cutoff,
         level == "BOT") %>% 
  filter(classcode == "SPUL") %>% 
  mutate(biomass = (count*((0.0144)*fish_tl^3.04))) %>% 
  group_by(campus, method, survey_year, year, month, day, site, zone, transect) %>% 
  # total length is in cm and biomass is in g
  summarize(total_count = sum(count), 
            total_biomass = sum(biomass, na.rm = T)) %>% 
  ungroup()

fish_densities <- distinct_transects_fish %>% 
  left_join(fish, by = c("campus", "method", "survey_year", "year", "month", "day", "site", "zone", "transect")) %>% 
  replace_na(list(total_count = 0,
                  total_biomass = 0)) %>% 
  group_by(year, mpa_status, site, region) %>% 
  summarise(SPUL_d = mean(total_count),
            biomass_d = mean(total_biomass)) %>% 
  ungroup()

###### 4) Format Inverts Data #####################################################
inverts <- mlpa_swath %>% 
  filter(year >= year_cutoff) %>% 
  filter(classcode == "MESFRAAD" |
           classcode == "STRPURAD" |
           classcode == "PANINT" |
           classcode == "MACPYRAD" |
           classcode == "NERLUE") %>% 
  pivot_wider(names_from = classcode, values_from = count) %>% 
  group_by(campus, method, survey_year, year, month, day, site, zone, transect) %>% 
  summarize(MESFRAAD = sum(MESFRAAD, na.rm = T),
            STRPURAD = sum(STRPURAD, na.rm = T),
            PANINT = sum(PANINT, na.rm = T),
            MACPYRAD = sum(MACPYRAD, na.rm = T),
            NERLUE = sum(NERLUE, na.rm = T),
            kelp = sum(MACPYRAD, NERLUE, na.rm = T))

inverts$total_urchins <- rowSums(inverts[,c("MESFRAAD", "STRPURAD")], na.rm = TRUE, dims = 1)

invert_densities <- distinct_transects_inverts %>% 
  left_join(inverts, by = c("campus", "method", "survey_year", "year", "month", "day", "site", "zone", "transect")) %>% 
  replace_na(list(PANINT = 0,
                  MESFRAAD = 0,
                  STRPURAD = 0, 
                  total_urchins = 0,
                  MACPYRAD = 0,
                  NERLUE = 0,
                  kelp = 0)) %>% 
  group_by(year, mpa_status, site, region) %>% 
  summarise(n_transects = n(),
            urchin_total = sum(total_urchins),
            PANINT_d = mean(PANINT), 
            MESFRAAD_d = mean(MESFRAAD),
            STRPURAD_d = mean(STRPURAD),
            urchins_d = mean(total_urchins),
            MACPYRAD_d = mean(MACPYRAD),
            total_giant_kelp = sum(MACPYRAD), 
            NERLUE_d = mean(NERLUE),
            kelp_d = mean(kelp)) %>% 
  ungroup()

# Each individual transect is weighted equally due to the distinct transects

##### 5) Combine ###############################################################
data_summarized <- full_join(fish_densities, invert_densities, 
                             by = c("year", "mpa_status","site", "region")) %>% 
  mutate(heatwave = case_when(year < 2014 ~ "before", 
                              year > 2016 ~ "after", 
                              .default = "during"),
         heatwave = factor(heatwave, 
                           levels = c("before", "during", "after")))

###### 5) Descriptive Statistics ###############################################
# How many sites per region / mpa_status
data_summarized %>% 
  ungroup() %>% 
  distinct(site, region, mpa_status) %>% 
  group_by(region, mpa_status) %>% 
  summarize(n = n()) %>% 
  mutate(100*n/sum(n))

data_summarized %>% 
  ungroup() %>% 
  distinct(site, mpa_status) %>% 
  group_by(mpa_status) %>% 
  summarize(n = n()) %>% 
  mutate(100*n/sum(n))

data_summarized %>% 
  ungroup() %>% 
  distinct(site, region) %>%
  group_by(region) %>% 
  summarize(n = n())

write.csv(data_summarized, 
          "Processed_data/MLPA_data_summarized_wo_siteblocks.csv",
          row.names = F)

