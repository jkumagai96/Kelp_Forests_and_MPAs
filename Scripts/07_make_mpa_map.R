# Date: October 30th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Create a data table with MPA protection per point with kelp data
# Purpose: Create a map of kelp area 
# BIO 202: Ecological Statistics

##### Set up: packages #########################################################
# load packages 
library(tidyverse)
library(sf)
library(tidyverse)
library(ggspatial)
library("rnaturalearth")
library("rnaturalearthdata")
library(viridis)

##### Set up: load data ########################################################
mpas_original <- read_sf("Data/Filtered_MPAs/MPAs_CA_Giant_Kelp.shp")
usa <- read_sf("Data/gadm41_USA_shp/gadm41_USA_0.shp") %>% 
  st_transform(crs = 4326) 

all_kelp_data <- all_kelp_data <- read.csv("Processed_data/data_tables/kelp_data_w_mpas_per_quarter.csv")
  
us <- ne_countries(scale = "large", returnclass = "sf") %>% 
  filter(name == "United States") %>% 
  st_transform(4326) 

##### Formatting ###############################################################
# Select needed attribuets from MPAs
mpas <- mpas_original %>% 
  dplyr::select(Site_ID_12, Estab_Yr_1, AreaMar_12, Status_12)  %>% # I selected Area_Mar12, and not Area? QUESTION
  rename("mpa_status" = "Status_12") %>% 
  # Match projections with kelp data 
  st_transform(crs = 4326) 

# Format kelp area data
all_kelp_data$log_area <- log(all_kelp_data$area + 1)

# Decisions: Focus on 2021
kelp_data <- all_kelp_data %>% 
  filter(year == 2021) %>% 
  filter(quarter == "Q1") %>% 
  # filter(lat <= 37 & lat >= 35) %>% 
  filter(area != 0 )

zeros <- all_kelp_data %>% 
  filter(year == 2021) %>% 
  filter(quarter == "Q1") %>% 
  # filter(lat <= 37 & lat >= 35) %>% 
  filter(area == 0 )
  
points <- st_as_sf(kelp_data, 
                   coords = c("long","lat"), 
                   crs = 4326,
                   remove = FALSE) 

zeros <- st_as_sf(zeros, 
                  coords = c("long", "lat"), 
                  crs = 4326, 
                  remove = FALSE)


poly_cropped <- us %>% 
  st_crop(st_bbox(points) + c(-.1, -.2, .1, .2))

##### Map Attempts #############################################################
ggplot() +
  geom_sf(data = poly_cropped, fill = "#DABC94") +
  geom_sf(size = 1, color = "grey",  data = zeros, alpha = 0.5) +
  #geom_sf(aes(color = log_area), size = points$log_area/3, data = points) +
  scale_color_viridis() +
  theme(panel.grid.major = element_line(color = gray(.5), size = 0.2, linetype = 'dashed'), 
        panel.background = element_rect(fill = 'black'))

ggsave("Santa_cruz_zeros.png", dpi = 600)

ggplot() + 
  geom_sf(data = poly_cropped, fill = "#DABC94") +
  # geom_sf(data = mpas, fill = "grey", color = "black") +
  geom_sf(size = 1, color = "black",  data = zeros, alpha = 0.5) +
  geom_sf(aes(color = log_area), size = 1, data = points) +
  # geom_sf(aes(color = log_area), size = points$log_area/3, data = points) +
  scale_color_viridis() +
  theme(panel.grid.major = element_line(color = gray(.5), size = 0.2, linetype = 'dashed'), 
        panel.background = element_rect(fill = 'white'))

ggsave("Santa_kelp_data.png", dpi = 600)