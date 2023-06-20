# Date: June 13th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Explore the PISCO data for presence of california sheephead and differences between protected and not protected
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(readxl)
library(sf)

# Load Data
PISCO_fish_raw <- read.csv("Data/PISCO subtidal data clearinghouse/PISCO_FISH.csv")
PISCO_fish_lookup <- read.csv("Data/PISCO subtidal data clearinghouse/PISCO_taxon_lookup_table.csv")
PISCO_master_species <- read_excel("Data/PISCO subtidal data clearinghouse/master_spp_table.xlsx")
PISCO_master_sites <- read_excel("Data/PISCO subtidal data clearinghouse/PISCO subtidal master site table.xlsx")

# mpas
mpas <- st_read("Processed_data/MPAs.shp")
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
  select(SITE, LATITUDE, LONGITUDE, Site_ID_12, Estab_Yr_1, AreaMar_12, mpa_status)

site_points <- site_points %>% 
  select(-REGION) %>% 
  left_join(sites_in_mpas, by = c("SITE", "LATITUDE", "LONGITUDE"))

site_points$mpa_status[is.na(site_points$mpa_status)] <- "Reference"

sites_for_joining <- site_points %>% 
  st_drop_geometry() %>% 
  select(SITE, site_status, Site_ID_12, Estab_Yr_1, AreaMar_12, mpa_status, region)

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
distinct_tansects <- PISCO_fish_raw %>% 
  distinct(campus, method, year, month, day, site, zone, transect) %>% 
  left_join(sites_for_joining, by = c("site" = "SITE")) %>% 
  filter(region == "South_Coast") 

## Account for Establishment year 
not_protected <- distinct_tansects %>% 
  filter(mpa_status == "Reference") 

subset <- distinct_tansects %>% 
  filter(mpa_status != "Reference")

for (i in 1:nrow(subset)) {
  status <- subset$mpa_status[i]
  establishment_yr <- subset$Estab_Yr_1[i]
  
  if (subset$year[i] >= establishment_yr) {
    subset$mpa_status[i] <- status
  } else {
    subset$mpa_status[i] <- "Reference"
  }
}

distinct_transects <- rbind(not_protected, subset)

sheephead_standardized <- distinct_transects %>% 
  left_join(sheephead, by = c("campus", "method", "year", "month", "day", "site", "transect", "zone")) %>% 
  replace_na(list(n_sheephead = 0)) %>% 
  group_by(year, mpa_status, site) %>% 
  summarise(avg_sheephead_per_site = mean(n_sheephead)) %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_sheephead = mean(avg_sheephead_per_site),
            sd_density = sd(avg_sheephead_per_site))

group.colors <- c(Full = "#440154", Reference = "#FFBA00", Partial ="#21918c")

sheephead_plot <- ggplot(data = sheephead_standardized, aes(x = year, y = avg_sheephead, color = mpa_status)) +
  geom_line(aes(color = mpa_status), linewidth = 0.8) +
  geom_errorbar(aes(ymin = avg_sheephead-sd_density, ymax = avg_sheephead + sd_density),width = 0.2, linewidth = .1)+
  theme_bw() +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  labs(x = "Year", y = "Average Sheephead per site") +
  theme(legend.position = c(0.25, 0.85))
sheephead_plot 

# Investigate by size
quantile(sheephead_by_size$fish_tl, na.rm = T)

sheephead_large <- sheephead_by_size %>% 
  filter(fish_tl >= 36)

sheephead_standardized_l <- distinct_transects %>% 
  left_join(sheephead_large, by = c("campus", "method", "year", "month", "day", "site", "transect", "zone")) %>% 
  replace_na(list(n_sheephead = 0)) %>% 
  group_by(year, mpa_status, site) %>% 
  summarise(avg_sheephead_per_site = mean(n_sheephead)) %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_sheephead = mean(avg_sheephead_per_site),
            sd_density = sd(avg_sheephead_per_site))

sheephead_plot_l <- ggplot(data = sheephead_standardized_l, aes(x = year, y = avg_sheephead, color = mpa_status)) +
  geom_line(aes(color = mpa_status), linewidth = 0.8) +
  geom_errorbar(aes(ymin = avg_sheephead-sd_density, ymax = avg_sheephead + sd_density),width = 0.2, linewidth = 0.2)+
  theme_bw() +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  labs(x = "Year", y = "Average large sheephead per site (>36cm)") +
  theme(legend.position = c(0.25, 0.85))
sheephead_plot_l 
  
# small
sheephead_small <- sheephead_by_size %>% 
  filter(fish_tl < 18)

sheephead_standardized_s <- distinct_transects %>% 
  left_join(sheephead_small, by = c("campus", "method", "year", "month", "day", "site", "transect", "zone")) %>% 
  replace_na(list(n_sheephead = 0)) %>% 
  group_by(year, mpa_status, site) %>% 
  summarise(avg_sheephead_per_site = mean(n_sheephead)) %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_sheephead = mean(avg_sheephead_per_site)) 

sheephead_plot_s <- ggplot(data = sheephead_standardized_s, aes(x = year, y = avg_sheephead, group = mpa_status)) +
  geom_line(aes(color = mpa_status), linewidth = 0.8) +
  theme_bw() +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  labs(x = "Year", y = "Average small sheephead per site (< 18 cm)") +
  theme(legend.position = c(0.25, 0.85))
sheephead_plot_s 
##### Export ###################################################################
png("Figures/Sheephead.png", width = 5, height = 5, 
    units = "in", res = 600)
sheephead_plot 
dev.off() 

png("Figures/Sheephead_large.png", width = 5, height = 5, 
    units = "in", res = 600)
sheephead_plot_l 
dev.off() 

png("Figures/Sheephead_small.png", width = 5, height = 5, 
    units = "in", res = 600)
sheephead_plot_s
dev.off() 

site_points <- site_points %>% 
  rename(SITE_CATEGORY = "SITE_CATEGORY__l-t__Project") 
st_write(site_points, "Processed_data/Subtidal/PISCO_site_points.shp", append = F)
