# Date: October 18th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Clean the fish data and extract biomass estimates
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(stringr)
library(sf)

# Load Data
reef_check_inverts <- read.csv("Data/Reef Check Data/Inverts.csv")
reef_check_fish <- read.csv("Data/Reef Check Data/Fish.csv")
survey_ids <- read.csv("Processed_data/data_tables/reef_check_suveyIDs_to_keep.csv")
subtidal_data <- read.csv("Processed_data/data_tables/subtidal_surveys.csv") 

# Declare Functions
std <- function(x) sd(x)/sqrt(length(x))

# MPA data 
mpas_original_south <- read_sf("Data/Filtered_MPAs/MPAs_CA_Giant_Kelp.shp") %>% 
  st_transform(crs = 4326) %>% 
  dplyr::select(Site_ID_12, Estab_Yr_1, AreaMar_12, Status_12)

##### Filter Data ##############################################################
### Clean Data
inverts <- reef_check_inverts %>% 
  distinct() %>%                                    # Remove duplicates 
  filter(region == "California") %>%                # Filter to just California
  dplyr::select(-region) %>% 
  filter(latitude <= 34.4486) %>%                   # Filter to just southern california and add name
  mutate(region = "South_Coast")%>% 
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
  filter(latitude <= 34.4486) %>%                   # Filter to just southern california and add name
  mutate(region = "South_Coast")%>% 
  mutate(year = str_sub(date_start, start= -4)) %>% # Add year to the dataset based on the date of the survey
  mutate(year = as.numeric(year)) %>% 
  dplyr::select(survey_id, site_name, latitude, longitude, # Variables we want to include
                date_start, transect, start_time,
                species, amount,region,
                year, min_cm, max_cm, size_category) %>% 
  mutate(area_searched = 60) %>%                    # Setting the area searched to 60m2 
  filter(amount >= 0) %>%                             # Amount can't be less than zero or NA
  na.omit(site_name)                               # Site has to be named 


# Only include the fish data which within the sites we are already analyzing  
fish <- fish %>% 
  left_join(survey_ids, by = c("survey_id" = "unique_survey_id")) %>% 
  filter(include_survey == 1) 

# Seperate analysis of sheephead biomass is needed as 23% of the data has a min value
sheephead <- fish %>% 
  filter(species == "Sheephead (female)" |
           species == "Sheephead (male)" |
           species == "Sheephead (juvenile)") 
nrow(sheephead)

sheephead %>% 
  filter(min_cm > 0) %>% 
  nrow() / nrow(sheephead)

sheephead_biomass <- sheephead %>%
  ungroup() %>% 
  filter(min_cm > 0) %>% 
  mutate(mean_cm = rowMeans(subset(., select = c(min_cm, max_cm)))) %>% 
  mutate(biomass = (amount*((0.0144)*mean_cm^3.04))) %>% 
  mutate(biomass_per_m2 = biomass/area_searched) %>% 
  group_by(survey_id, site_name, latitude, longitude, date_start, 
           transect, start_time, year, region, area_searched) %>% 
  summarise(g_per_m2_per_transect = mean(biomass_per_m2)) 

sheephead_biomass_per_survey <- sheephead_biomass %>% 
  group_by(survey_id, site_name, latitude, longitude, date_start, 
           year, region) %>% 
  summarise(g_per_m2_per_survey = mean(g_per_m2_per_transect))

# Sites are weighted equally?
# Sites are not weighted equally? Does it matter? Let's do it both ways 

##### Join MPA Data per site ###################################################
# Select needed information from subtidal data

site_mpa_data <- subtidal_data %>% 
  select(site_name, mpa_status, year) %>% 
  distinct()

biomass_data <- left_join(sheephead_biomass_per_survey, site_mpa_data, by = c("site_name", "year"))  
biomass_per_transect <- left_join(sheephead_biomass, site_mpa_data, by = c("site_name", "year"))

##### Create Figures ###########################################################
group.colors <- c(Full = "#440154", None = "#FFBA00", Partial ="#21918c")

# Sheephead Biomass averaging across sites 
plot_data <- biomass_data %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_sheephead_biomass = mean(g_per_m2_per_survey),
            se_biomass = std(g_per_m2_per_survey))

plot_biomass <- plot_data %>% 
  ggplot(aes(x = year, y = avg_sheephead_biomass, color = mpa_status)) +
  geom_line(linewidth = 1, show.legend = FALSE) +
  geom_errorbar(aes(ymin = avg_sheephead_biomass-se_biomass, ymax = avg_sheephead_biomass + se_biomass), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 0.4)+
  scale_x_continuous(breaks = seq(2014, 2022, by = 2)) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab(bquote('CA Sheephead biomass')) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  theme_bw() 
plot_biomass

# Quick statistics 
plot_data$mpa_status <- as.factor(plot_data$mpa_status)
res_aov <- aov(avg_sheephead_biomass ~ mpa_status,
               data = plot_data
)
# Are the residuals normal? 
hist(res_aov$residuals)

# QQ-plot
library(car)
qqPlot(res_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)

summary(res_aov)

# Tukey HSD test:
library(multcomp)
post_test <- glht(res_aov,
                  linfct = mcp(mpa_status = "Tukey")
)

summary(post_test)
# 
# 
# # Sheephead Biomass averaging across transects (doesn't make a difference)
# plot_biomass_transect <- biomass_per_transect %>% 
#   group_by(year, mpa_status) %>% 
#   summarise(avg_sheephead_biomass = mean(g_per_m2_per_transect)) %>% 
#   ggplot(aes(x = year, y = avg_sheephead_biomass, color = mpa_status)) +
#   geom_line(linewidth = 1, show.legend = FALSE) +
#   scale_x_continuous(breaks = seq(2014, 2022, by = 2)) +
#   scale_color_manual(values=group.colors, name = "MPA Category") +
#   ylab(bquote('CA Sheephead biomass - g/' ~ m^2)) +
#   xlab("Year") +
#   annotate("rect", fill = "red", alpha = 0.2, 
#            xmin = 2014, xmax = 2016,
#            ymin = -Inf, ymax = Inf) +
#   theme_bw() 
# 
# # Therefore we should go with equally weighting the sites as we do that previously 

##### Export Figure ############################################################
png("Figures/Sheephead_biomass_reefcheck.png", width = 6, height = 4, 
    units = "in", res = 600)
plot_biomass
dev.off() 

# Yay! MPAs work :) 