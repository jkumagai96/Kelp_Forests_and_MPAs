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

PISCO_fish_lookup %>% 
  filter(common_name == "California Sheephead") 

sheephead <- PISCO_fish_raw %>% 
  filter(classcode == "SPUL") %>% 
  full_join(sites_for_joining, by = c("site" = "SITE")) %>% 
  filter(region == "South_Coast") %>% # 97% of the data is within southern california 
  group_by(campus, method, year, month, day, site, transect, zone) %>% 
  summarise(n_sheephead = sum(count)) %>% 
  filter(!is.na(n_sheephead))

sheephead_by_size <- PISCO_fish_raw %>% 
  filter(classcode == "SPUL") %>% 
  full_join(sites_for_joining, by = c("site" = "SITE")) %>% 
  filter(region == "South_Coast") %>% # 97% of the data is within southern california 
  group_by(campus, method, year, month, day, site, transect, zone, fish_tl) %>% 
  summarise(n_sheephead = sum(count)) %>% 
  filter(!is.na(n_sheephead))

# The data is not taking into account when sheephead are not seen, super important due to low numbers
# Need to create a presence absence table...

# Count the number of distinct transects 
distinct_transects <- PISCO_fish_raw %>% 
  distinct(campus, method, year, month, day, site, zone, transect) %>% 
  full_join(sites_for_joining, by = c("site" = "SITE")) %>% 
  filter(region == "South_Coast") 

sheephead_standardized <- distinct_transects %>% 
  full_join(sheephead, by = c("campus", "method", "year", "month", "day", "site", "transect", "zone")) %>% 
  replace_na(list(n_sheephead = 0)) %>% 
  group_by(year, mpa_status, site) %>% 
  summarise(avg_sheephead_per_site = mean(n_sheephead)) %>% # sites are equally weighted
  group_by(year, mpa_status) %>% 
  summarise(avg_sheephead = mean(avg_sheephead_per_site),
            se_density = std(avg_sheephead_per_site)) %>% 
  filter(year >= 1999)

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
sheephead_standardized_biomass <- distinct_transects %>% 
  left_join(sheephead_by_size, by = c("campus", "method", "year", "month", "day", "site", "transect", "zone")) %>% 
  mutate(biomass = (n_sheephead*((0.0144)*fish_tl^3.04))) %>% 
  replace_na(list(biomass = 0)) %>%
  group_by(year, mpa_status, site) %>% 
  summarise(avg_biomass_per_site = mean(biomass)) %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_biomass = mean(avg_biomass_per_site),
            se_biomass = std(avg_biomass_per_site)) %>% 
  filter(year >= 2007)

print("Refer to PISCO figures script to visualize biomass")
write.csv(sheephead_standardized_biomass, "Processed_data/Data_tables/PISCO_biomass_data.csv", row.names = F)

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


l_sheephead_standardized_biomass_per_site <- distinct_transects %>% 
  left_join(sheephead_by_size, by = c("campus", "method", "year", "month", "day", "site", "transect", "zone")) %>% 
  mutate(biomass = (n_sheephead*((0.0144)*fish_tl^3.04))) %>% 
  replace_na(list(biomass = 0)) %>%
  filter(fish_tl >= 35) %>% 
  group_by(year, mpa_status, site, Estab_Yr_1) %>% 
  summarise(avg_biomass_per_site = mean(biomass)) 

write.csv(l_sheephead_standardized_biomass_per_site, "Processed_data/Data_tables/PISCO_biomass_large_per_site.csv", row.names = F)

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
  labs(x = "Year", y = "Mean Sheephead Biomass (>= 35cm)") +
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


##### Export ###################################################################
png("Figures/Sheephead.png", width = 5, height = 5, 
    units = "in", res = 600)
sheephead_plot 
dev.off() 

png("Figures/Sheephead_cohort_biomass_plot.png", width = 9, height = 5, 
    units = "in", res = 600)
combo_plot
dev.off() 


