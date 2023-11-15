# Date: October 30th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Explore the PISCO data for presence of california sheephead and differences between protected and not protected
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(readxl)
library(sf)
library(cowplot)

# Load Data
PISCO_fish_raw <- read.csv("Data/PISCO subtidal data clearinghouse/PISCO_FISH.csv")
PISCO_fish_lookup <- read.csv("Data/PISCO subtidal data clearinghouse/PISCO_taxon_lookup_table.csv")
PISCO_master_species <- read_excel("Data/PISCO subtidal data clearinghouse/master_spp_table.xlsx")
PISCO_master_sites <- read_excel("Data/PISCO subtidal data clearinghouse/PISCO subtidal master site table.xlsx")

# mpas
mpas <- st_read("Processed_data/MPAs.shp")

# Declare Functions
std <- function(x) sd(x)/sqrt(length(x))

###### Format Data #############################################################
# See map of sites 
site_points <- PISCO_master_sites %>% 
  mutate(longitude = as.numeric(LONGITUDE), 
         latitude = as.numeric(LATITUDE)) %>% 
  mutate(region = ifelse(latitude > 34.4486, "Central_Coast", "South_Coast")) %>% 
  mutate(region = ifelse(latitude > 37.1819, "North_Central_Coast", region)) %>% 
  mutate(region = ifelse(latitude > 39.0044, "North_Coast", region)) %>% 
  filter(region == "South_Coast" | region == "Central_Coast") %>% 
  drop_na(longitude) %>% 
  st_as_sf(., coords = c(x = "longitude", y = "latitude"), crs = 4326) 

st_write(site_points, "Processed_data/PISCO_sites.shp", append = TRUE)

# Intersect site data with mpas 
sites_in_mpas <- st_intersection(site_points, mpas) %>% 
  st_drop_geometry() %>% 
  dplyr::select(SITE, LATITUDE, LONGITUDE, Site_ID_12, Estab_Yr_1, AreaMar_12, mpa_status)

site_points <- site_points %>% 
  select(-REGION) %>% 
  left_join(sites_in_mpas, by = c("SITE", "LATITUDE", "LONGITUDE"))

site_points$mpa_status[is.na(site_points$mpa_status)] <- "Reference"

sites_for_joining <- site_points %>% 
  st_drop_geometry() %>% 
  select(SITE, site_status, Site_ID_12, Estab_Yr_1, AreaMar_12, mpa_status, region)

# Export PISCO sites for joining 
write.csv(sites_for_joining, "Processed_data/data_tables/PISCO_sites_with_MPAs.csv", row.names = FALSE)

PISCO_fish_lookup %>% 
  filter(common_name == "California Sheephead") 

sheephead <- PISCO_fish_raw %>% 
  filter(classcode == "SPUL") %>% 
  left_join(sites_for_joining, by = c("site" = "SITE")) %>% 
  filter(region == "South_Coast") %>% # 97% of the data is within southern california 
  group_by(campus, method, year, month, day, site, transect, zone) %>% 
  summarise(n_sheephead = sum(count))

sheephead_by_size <- PISCO_fish_raw %>% 
  filter(classcode == "SPUL") %>% 
  left_join(sites_for_joining, by = c("site" = "SITE")) %>% 
  filter(region == "South_Coast") %>% # 97% of the data is within southern california 
  group_by(campus, method, year, month, day, site, transect, zone, fish_tl) %>% 
  summarise(n_sheephead = sum(count))

# The data is not taking into account when sheephead are not seen, super important due to low numbers
# Figure out how to create a presence absence table...
# Figure out statistics!

# Count the number of distinct transects 
distinct_transects <- PISCO_fish_raw %>% 
  distinct(campus, method, year, month, day, site, zone, transect) %>% 
  left_join(sites_for_joining, by = c("site" = "SITE")) %>% 
  filter(region == "South_Coast") 

sheephead_standardized <- distinct_transects %>% 
  left_join(sheephead, by = c("campus", "method", "year", "month", "day", "site", "transect", "zone")) %>% 
  replace_na(list(n_sheephead = 0)) %>% 
  group_by(year, mpa_status, site) %>% 
  summarise(avg_sheephead_per_site = mean(n_sheephead)) %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_sheephead = mean(avg_sheephead_per_site),
            se_density = std(avg_sheephead_per_site)) %>% 
  filter(year >= 2007)

group.colors <- c(Full = "#440154", Reference = "#FFBA00", Partial ="#21918c")

sheephead_plot <- ggplot(data = sheephead_standardized, aes(x = year, y = avg_sheephead, color = mpa_status)) +
  geom_line(aes(color = mpa_status), linewidth = 1) +
  geom_errorbar(aes(ymin = avg_sheephead-se_density, 
                    ymax = avg_sheephead + se_density), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 0.4)+
  theme_bw() +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  labs(x = "Year", y = "Mean Sheephead Abundance") +
  theme(legend.position = c(0.25, 0.75)) +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) + 
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = .7)

sheephead_plot 

##### Investigate by biomass ###################################################
sheephead_standardized_biomass_per_site <- distinct_transects %>% 
  left_join(sheephead_by_size, by = c("campus", "method", "year", "month", "day", "site", "transect", "zone")) %>% 
  mutate(biomass = (n_sheephead*((0.0144)*fish_tl^3.04))) %>% 
  replace_na(list(biomass = 0)) %>%
  group_by(year, mpa_status, site, Estab_Yr_1) %>% 
  summarise(avg_biomass_per_site = mean(biomass)) 

sheephead_standardized_biomass <- sheephead_standardized_biomass_per_site %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_biomass = mean(avg_biomass_per_site),
            se_biomass = std(avg_biomass_per_site)) %>% 
  filter(year >= 2007)

print("Refer to PISCO figures script to visualize biomass")
# write.csv(sheephead_standardized_biomass, "Processed_data/Data_tables/PISCO_biomass_data.csv", row.names = F)

###### Investigate by size #####################################################
# I would like to make a graph where we have a series of histograms of the size
# classes of sheephead, so x is length, y is frequency, and the group is year 
size_classes_df <- distinct_transects %>% 
  left_join(sheephead_by_size, by = c("campus", "method", "year", "month", "day", "site", "transect", "zone")) %>% 
  replace_na(list(n_sheephead = 0)) %>% 
  replace_na(list(fish_tl = 0)) %>% 
  filter(n_sheephead > 0) %>% 
  filter(year >= 2007) %>% 
  mutate(yr = as.factor(year))


t <- size_classes_df %>% 
  count(year)

sample_size <- as.character(t$n)
sample_size[15] <- paste("n =", sample_size[15])

library(ggridges)
cohort_plot <- ggplot(aes(x = fish_tl, y = yr, group = yr, fill = yr), data = size_classes_df) +
  annotate("rect", fill = "#ADD8E620",
           xmin = 30, xmax = Inf,
           ymin = -Inf, ymax = Inf) +
  geom_density_ridges2(panel_scaling = F, 
                       rel_min_height = 0.001) +
  theme_bw() +
  theme(legend.position = "none") + 
  geom_vline(xintercept = 30, color = "black") + 
  scale_fill_manual(values = c(rep("grey", 7), rep("firebrick2", 3), rep("grey", 5))) +
  annotate("text", 
           x = 90, y = 1.25:15.25, 
           size = 3, 
           hjust = 1,
           label = sample_size) +
  labs(x = "Sheephead Total Length", y = "Year") 

cohort_plot 

ggplot(aes(x = fish_tl, y = yr, group = yr, fill = yr), data = size_classes_df) +
  annotate("rect", fill = "#ADD8E620",
           xmin = 30, xmax = Inf,
           ymin = -Inf, ymax = Inf) +
  geom_ridgeline(height = 4) +
  theme_bw() +
  theme(legend.position = "none") + 
  geom_vline(xintercept = 30, color = "black") + 
  scale_fill_manual(values = c(rep("grey", 7), rep("firebrick2", 3), rep("grey", 5))) +
  labs(x = "Sheephead Total Length", y = "Year") 


l_sheephead_standardized_biomass_per_site <- distinct_transects %>% 
  left_join(sheephead_by_size, by = c("campus", "method", "year", "month", "day", "site", "transect", "zone")) %>% 
  mutate(biomass = (n_sheephead*((0.0144)*fish_tl^3.04))) %>% 
  replace_na(list(biomass = 0)) %>%
  filter(fish_tl >= 30) %>% 
  group_by(year, mpa_status, site, Estab_Yr_1) %>% 
  summarise(avg_biomass_per_site = mean(biomass)) 

l_sheephead_standardized_biomass <- l_sheephead_standardized_biomass_per_site %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_biomass = mean(avg_biomass_per_site),
            se_biomass = std(avg_biomass_per_site)) %>% 
  filter(year >= 2007)

high_biomass_plot <- ggplot(data = l_sheephead_standardized_biomass,
                            aes(x = year,
                                y = avg_biomass,
                                color = mpa_status)) +
  geom_line(aes(color = mpa_status), linewidth = 1) +
  geom_errorbar(aes(ymin = avg_biomass-se_biomass, ymax = avg_biomass + se_biomass),width = 0.2, linewidth = .5, alpha = 0.4)+
  theme_bw() +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  labs(x = "Year", y = "Mean Sheephead Biomass (>= 30cm)") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#ADD8E620")) +
  annotate("rect", fill = "red", alpha = 0.2,
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = .7)

high_biomass_plot

combo_plot <-plot_grid(cohort_plot, 
                       high_biomass_plot, 
                       nrow = 1, 
                       rel_widths = c(1.5,2),
                       labels = "AUTO")
combo_plot

##### Statistics ###############################################################
# Simple ANOVA on biomass (abundances are not significant)
sheephead_anova <- sheephead_standardized_biomass %>% 
  filter(year >= 2012) # after all the MPAs have been implemented

sheephead_anova$mpa_status <- as.factor(sheephead_anova$mpa_status)

res_aov <- aov(avg_biomass ~ mpa_status,
               data = sheephead_anova
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

### Two way fixed effects (look into later)
# effects_data <- sheephead_standardized_biomass_per_site %>% 
#   mutate(yr_si = year - Estab_Yr_1) %>%   # year since implementation
#   mutate(mpa_status = factor(mpa_status, levels = c("Reference", "Partial", "Full")))
# 
# m1_effects <- lm(avg_biomass_per_site ~ mpa_status*Estab_Yr_1 + mpa_status + Estab_Yr_1, 
#    data = effects_data)
# 
# summary(m1_effects)
# plot(m1_effects)
##### Export ###################################################################
png("Figures/Sheephead.png", width = 5, height = 5, 
    units = "in", res = 600)
sheephead_plot 
dev.off() 

png("Figures/Sheephead_cohort_biomass_plot.png", width = 9, height = 5, 
     units = "in", res = 600)
combo_plot
dev.off() 

