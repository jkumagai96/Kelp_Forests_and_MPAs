# Date: August 1st 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Clean and combine the reef check invertebrate and fish data 
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(stringr)
library(sf)

# Load Data
reef_check_inverts <- read.csv("Data/Reef Check Data/Inverts.csv")
reef_check_fish <- read.csv("Data/Reef Check Data/Fish.csv")

##### Set up: load data ########################################################
mpas_original <- st_read("Processed_data/MPAs.shp")

##### Filter Data ##############################################################
### Clean Data
inverts <- reef_check_inverts %>% 
  distinct() %>%                                    # Remove duplicates 
  filter(region == "California") %>%                # Filter to just California
  dplyr::select(-region) %>% 
  filter(latitude <= 37.1819) %>%                   # Filter to just southern and central california and add name
  mutate(region = ifelse(latitude > 34.4486, "Central_Coast", "South_Coast")) %>% 
  mutate(year = str_sub(date_start, start= -4)) %>% # Add year to the dataset based on the date of the survey
  mutate(year = as.numeric(year)) %>% 
  dplyr::select(survey_id, site_name, latitude, longitude, # Variables we want to include
         date_start, transect, start_time,
         species, amount, distance,
         year, region) %>% 
  filter(amount >= 0) %>%                             # Amount can't be less than zero or NA
  na.omit(site_name) %>%                              # Site has to be named 
  mutate(density_inverts = amount/(distance*2))

fish <- reef_check_fish %>% 
  distinct() %>%                                    # Remove duplicates 
  filter(region == "California") %>%                # Filter to just California
  dplyr::select(-region) %>% 
  filter(latitude <= 37.1819) %>%                   # Filter to just southern and central california and add name
  mutate(region = ifelse(latitude > 34.4486, "Central_Coast", "South_Coast")) %>% 
  mutate(year = str_sub(date_start, start= -4)) %>% # Add year to the dataset based on the date of the survey
  mutate(year = as.numeric(year)) %>% 
  dplyr::select(survey_id, site_name, latitude, longitude, # Variables we want to include
         date_start, transect, start_time,
         species, amount,region,
         year, min_cm, max_cm, size_category) %>% 
  mutate(area_searched = 60) %>%                    # Setting the area searched to 60m2 
  filter(amount >= 0) %>%                             # Amount can't be less than zero or NA
  na.omit(site_name)                               # Site has to be named 

### Only include survey's that have both fish and invertebrate data 

# Create dataset with one row being survey id, and one row being yes that survey is in inverts, and second row with yes that survey is in fish as well
survey_id_inverts <- data.frame(unique_survey_id = unique(inverts$survey_id), 
                                inverts = 1)

survey_id_fish <- data.frame(unique_survey_id = unique(fish$survey_id), 
                             fish = 1)

survey_id_df <- inner_join(survey_id_inverts, survey_id_fish, by = "unique_survey_id") %>% 
  mutate(include_survey = 1) %>% 
  dplyr::select(unique_survey_id, include_survey)

# Add this as a yes or no to include (left_join), then remove all the other records
inverts <- inverts %>% 
  left_join(survey_id_df, by = c("survey_id" = "unique_survey_id")) %>% 
  filter(include_survey == 1)
  
fish <- fish %>% 
  left_join(survey_id_df, by = c("survey_id" = "unique_survey_id")) %>% 
  filter(include_survey == 1)

# This process removes some of the fish records, but none of the inverts
write.csv(survey_id_df, "Processed_data/data_tables/reef_check_suveyIDs_to_keep.csv", row.names = F)

##### Standardize and Combine data #############################################
### Invertebrates
# Create unique transects
distinct_transects <- inverts %>% 
  distinct(survey_id, site_name, latitude, longitude, date_start, 
           transect, start_time, year, region, distance)

# Check Sunflower star abundances to assess whether an analysis is possible (11 is not enough)
sun_stars <- inverts %>% 
  filter(species == "Sunflower Star" | 
           species == "Sun Star" | 
           species == "Sunflower/Sun Star") %>% 
  filter(amount > 0)

invert_oi_density_per_transect <- inverts %>% 
  filter(distance > 0) %>% 
  filter(species == "Crowned Urchin" |              # Filter to crowned, purple, and red urchins
           species == "Purple Urchin" | 
           species == "Red Urchin" | 
           species == "Red Urchin (<2.5)" |
           species == "California Spiny Lobster") %>% # Include the lobsters
  group_by(survey_id, site_name, latitude, longitude, date_start, 
           transect, start_time, year, region, species) %>% 
  summarise(mean_density = mean(density_inverts)) %>% 
  pivot_wider(names_from = species, 
              values_from = mean_density) %>% 
  replace_na(list('California Spiny Lobster' = 0, 
                 'Red Urchin' = 0, 
                 'Purple Urchin' = 0, 
                 'Crowned Urchin' = 0)) # Replaces the NULLs with zeros!!! 

# This dataset is produced to account the sum of urchins per transect 
urchin_density_per_transect <- inverts %>% 
  filter(species == "Crowned Urchin" |              # Filter to crowned, purple, and red urchins
           species == "Purple Urchin" | 
           species == "Red Urchin" | 
           species == "Red Urchin (<2.5)") %>% 
  # Remove any records where distance is equal to 0
  filter(distance > 0) %>% 

# Need to have a standard amount counted per distance   
  mutate(amount_adjusted = ifelse(distance == 30, 
                                  amount, 
                                  ((amount*30)/distance)
                                  )
         ) %>%
  group_by(survey_id, site_name, latitude, longitude, date_start, 
           transect, start_time, year, region) %>%
  summarise(total_urchins = sum(amount_adjusted)) %>% 
  mutate(urchin_density = total_urchins/60) 
  

inverts_standard <- distinct_transects %>% 
  left_join(invert_oi_density_per_transect) %>% 
  left_join(urchin_density_per_transect) %>% 
  replace_na(list('California Spiny Lobster' = 0, 
                  'Red Urchin' = 0, 
                  'Purple Urchin' = 0, 
                  'Crowned Urchin' = 0,
                  'total_urchins' = 0,
                  'urchin_density' = 0))
# No red urchins less than 2.5 cm included

invert_surveys <- inverts_standard %>% 
  group_by(survey_id, site_name, latitude, longitude, date_start, 
           year, region) %>% 
  summarise(avg_r_urchin = mean(`Red Urchin`),
            avg_p_urchin = mean(`Purple Urchin`),
            avg_c_urchin = mean(`Crowned Urchin`),
            avg_urchins = mean(urchin_density),
            avg_lobster = mean(`California Spiny Lobster`))

### Fish Data
# Create unique transects
distinct_transects <- fish %>% 
  distinct(survey_id, site_name, latitude, longitude, date_start, 
           transect, start_time, year, region, area_searched)

sheephead_density_per_transect <- fish %>% 
  filter(species == "Sheephead (female)" |
           species == "Sheephead (male)" |
           species == "Sheephead (juvenile)") %>%  # Subset fish for sheephead
  group_by(survey_id, site_name, latitude, longitude, date_start, 
           transect, start_time, year, region, area_searched) %>% 
  summarise(total_sheephead = sum(amount)) %>% 
  mutate(density_sheephead = total_sheephead/area_searched)

fish_standard <- distinct_transects %>% 
  left_join(sheephead_density_per_transect) %>% 
  replace_na(list(total_sheephead = 0, 
                  density_sheephead = 0))

fish_surveys <- fish_standard %>% 
  group_by(survey_id, site_name, latitude, longitude, date_start, 
           year, region) %>% 
  summarise(avg_sheephead = mean(density_sheephead))

### Combine data
surveys <- inner_join(invert_surveys, fish_surveys) # 2 records are removed, but there is an incorrect place name 

# There is one survey a year 

##### Join Site Data ###########################################################
# Select needed attributes from MPAs
mpas <- mpas_original

# Sites 
sites <- surveys %>% 
  ungroup %>% 
  distinct(latitude, longitude, site_name, region)

# convert site names into spatial data to intersect with mpas 
site_points <- st_as_sf(sites, 
                        coords = c("longitude", "latitude"), 
                        crs = 4326) 

# export sites 
st_write(site_points, "Processed_data/Reef_check_sites.shp", append = TRUE)

# Intersect sites with MPA data to get all the sites within MPAs 
sites_in_mpas <- st_intersection(site_points, mpas) %>% 
  st_drop_geometry() %>% 
  select(-region)

site_points <- site_points %>% 
  left_join(sites_in_mpas, by = "site_name") 

site_points$mpa_status[is.na(site_points$mpa_status)] <- "None"

# Make sure the projections match 
st_crs(site_points) == st_crs(mpas)

# Drop geometry
sites_final <- site_points %>% 
  st_drop_geometry()

surveys_w_mpa_data <- surveys %>% 
  full_join(sites_final, by = c("site_name", "region")) 

##### Account for establishment year ###########################################

not_protected <- surveys_w_mpa_data %>% 
  filter(mpa_status == "None") 

subset <- surveys_w_mpa_data %>% 
  filter(mpa_status != "None")

for (i in 1:nrow(subset)) {
  status <- subset$mpa_status[i]
  establishment_yr <- subset$Estab_Yr_1[i]
  
  if (subset$year[i] >= establishment_yr) {
    subset$mpa_status[i] <- status
  } else {
    subset$mpa_status[i] <- "None"
  }
}

surveys_final <- rbind(not_protected, subset)

##### Export dataset ###########################################################
write.csv(surveys_final, "Processed_data/data_tables/subtidal_surveys.csv", row.names = F)
