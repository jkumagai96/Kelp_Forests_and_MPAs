# Date: August 2nd 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Create recuritment plot
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(ggridges)

# Load Data
MLPA_fish_raw <- read.csv("Data/resourceMap_MLPA_kelpforest_7/data/MLPA_kelpforest_fish.6.csv")
MLPA_sites <- read.csv("Processed_data/data_tables/MLPA_sites_with_MPAs.csv") 

# Sites for joining


##### Format Data ##############################################################
year_cutoff <- 2002

# Filter for just southern California
MLPA_sites <- MLPA_sites %>% 
  filter(region == "South_Coast")

### 1) Adjust transect names
transect_df <- data.frame(transect = unique(MLPA_fish_raw$transect),
                          correct_t = c(3,4,6,5,7,1,2,8,9,10,12,13,14,11,21,22,23,24,25,26,41,51,42,52,43,53,44,54,61,62,63,64,27,28,30,31,32,29,
                                        71,72,73,74,75,76)) 

mlpa_fish <- MLPA_fish_raw %>% 
  full_join(transect_df, by = "transect") %>% 
  dplyr::select(-transect) %>% 
  rename(transect = correct_t)


### 2) Create unique transects 

# Remove sites we don't want in the analysis! If they are not in the master list they are left out of the analysis 
unique_sites <- unique(MLPA_sites$site) 
mlpa_fish <- mlpa_fish %>% 
  filter(site %in% unique_sites)

distinct_transects_fish <- mlpa_fish %>% 
  filter(year >= year_cutoff,
         level == "BOT") %>% 
  distinct(campus, method, survey_year, year, month, day, site, zone, transect) %>% 
  left_join(MLPA_sites, by = "site")

### 3) Format Fish Data ######################################################
sheephead_by_size <- mlpa_fish %>% 
  filter(fish_tl >= 15) %>% 
  filter(year >= year_cutoff,
         level == "BOT") %>% 
  filter(classcode == "SPUL") %>% 
  mutate(biomass = (count*((0.0144)*fish_tl^3.04))) %>% 
  group_by(campus, method, survey_year, year, month, day, site, zone, transect, fish_tl) %>% 
  # total length is in cm and biomass is in g
  summarize(total_count = sum(count), 
            total_biomass = sum(biomass, na.rm = T)) %>% 
  ungroup()


###### Investigate by size #####################################################
# I would like to make a graph where we have a series of histograms of the size
# classes of sheephead, so x is length, y is frequency, and the group is year 

# Make a dataframe where each row is one observation of one fish
# Include in another dataset the number of unique transects that were run for each year.

sampling_effort <- distinct_transects_fish %>% 
  count(year)
sample_size <- as.character(sampling_effort$n)
sample_size[22] <- paste("n =", sample_size[22])

df <- sheephead_by_size %>% 
  select(year, fish_tl, total_count) %>% 
  uncount(total_count) %>% 
  mutate(year = as.factor(year))

plot_cohort <- ggplot(aes(x = fish_tl, y = year, group = year, fill = year), data = df) +
  geom_density_ridges2(panel_scaling = F, 
                       rel_min_height = 0.001) +
  theme_bw() +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c(rep("#80808090", 12), rep("#ee2c2c90", 3), rep("#80808090", 7))) +
  labs(x = "Sheephead Total Length", y = "Year") +
  geom_vline(xintercept = 15, color = "black") + 
  annotate("text", 
           x = 90, y = 1.25:22.25, 
           size = 3, 
           hjust = 0.6,
           label = sample_size)
plot_cohort

##### Export ###################################################################
png("Figures/Sheephead_cohort.png", width = 5, height = 5, 
    units = "in", res = 600)
plot_cohort
dev.off() 


