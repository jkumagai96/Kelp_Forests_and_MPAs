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
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)

##### Set up: load data ########################################################
mpas_original <- read_sf("Data/Filtered_MPAs/MPAs_CA_Giant_Kelp.shp")
usa <- read_sf("Data/gadm41_USA_shp/gadm41_USA_1.shp") %>% 
  st_transform(crs = 4326) 

CA <- usa %>% filter(NAME_1 == "California")

all_kelp_data <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")
  
us <- ne_countries(scale = "medium", returnclass = "sf") %>% 
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

# Decisions: Focus on 2015
kelp_data <- all_kelp_data %>% 
  filter(year == 2015) %>% 
  # filter(lat <= 37 & lat >= 35) %>% 
  filter(area != 0 )

zeros <- all_kelp_data %>% 
  filter(year == 2015) %>% 
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
  st_crop(st_bbox(points) + c(-1, -.2, .1, 1))

##### Map Attempts #############################################################
ggplot() + 
  #geom_sf(data = poly_cropped, fill = "grey80") +
  geom_sf(data = CA, fill = "grey70") +
  geom_sf(data = mpas, color = "red", fill = "red") +
  geom_sf(size = 1, color = "black",  data = zeros) +
  geom_sf(aes(color = log_area), size = 1, data = points) +
  scale_color_viridis() +
  theme(panel.grid.major = element_line(color = gray(.5), size = 0.2, linetype = 'dashed'), 
        panel.background = element_rect(fill = 'white'))

ggsave(filename = "Figures/CA_kelp_2015.png", dpi = 600)


ggplot() + 
  geom_sf(data = poly_cropped, fill = "grey80") +
  geom_sf(data = mpas, color = "red", fill = "red") +
  geom_sf(size = 1, color = "black",  data = zeros) +
  geom_sf(aes(color = log_area), size = 1, data = points) +
  scale_color_viridis() +
  theme(panel.grid.major = element_line(color = gray(.5), size = 0.2, linetype = 'dashed'), 
        panel.background = element_rect(fill = 'white'))

ggsave(filename = "Figures/CA_just_kelp_2015.png", dpi = 600)
