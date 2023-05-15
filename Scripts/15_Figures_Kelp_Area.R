# Date: May 10th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Create Kelp area per year figures 
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")
points_in_mpas <- read.csv("Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv") 

##### Explore MPA distributions ################################################
# Within our dataset, the max MPA establishment year is 2012 
# MPA area in 2021, how is it distributed?

kelp_data_all %>% 
  filter(year == 2021) %>% 
  group_by(region, mpa_status) %>% 
  count(mpa_status)

# So the number of pixels witin MPAs are basically the same, excpet there are more "None" pixels in
# Southern California

# Test actual kelp area 
mpas <- points_in_mpas %>% 
  select(PixelID, Site_ID_12, AreaMar_12, region) %>% 
  unique()
# 66 total MPAs 

mpas %>% 
  group_by(region) %>% 
  summarise(sum(AreaMar_12)) # there's 30% more mpa area in the souterh 

##### Kelp Area Figures ########################################################
# Kelp area through time 
df <- kelp_data_all %>% 
  group_by(year) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area),
            sd_kelp_area = sd(area)) 

ggplot(df, aes(x = year, y = mean_kelp_area)) + 
  geom_line(color = "green") +
  geom_point(color = "darkgreen") +
  theme_bw() 
  
ggplot(df, aes(x = year, y = median_kelp_area)) + 
  geom_line(color = "orange") +
  geom_point(color = "darkorange") +
  theme_bw() 

# Kelp area through time split by region 
df2 <- kelp_data_all %>% 
  group_by(year, region) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area)) 

ggplot(df2, aes(x = year, y = mean_kelp_area, group = region)) + 
  geom_line(aes(linetype=region))+
  geom_point(aes(shape=region)) + 
  theme_bw() 

# Kelp area through time split by mpa status

df3 <- kelp_data_all %>% 
  group_by(year, mpa_status) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area)) 

ggplot(df3, aes(x = year, y = mean_kelp_area, group = mpa_status)) + 
  geom_line(aes(color=mpa_status))+
  geom_point(aes(color=mpa_status)) + 
  theme_bw() 

ggplot(df3, aes(x = year, y = median_kelp_area, group = mpa_status)) + 
  geom_line(aes(color=mpa_status))+
  geom_point(aes(color=mpa_status)) + 
  theme_bw() 

# Just looking at southern california, kelp area through time split by mpa stauts
df4 <- kelp_data_all %>% 
  filter(region == "South_Coast") %>% 
  group_by(year, mpa_status) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area)) 

ggplot(df4, aes(x = year, y = mean_kelp_area, group = mpa_status)) + 
  geom_line(aes(color=mpa_status))+
  geom_point(aes(color=mpa_status)) + 
  theme_bw() 

ggplot(df4, aes(x = year, y = median_kelp_area, group = mpa_status)) + 
  geom_line(aes(color=mpa_status))+
  geom_point(aes(color=mpa_status)) + 
  theme_bw() 

# Just looking at central california, kelp area through time split by mpa stauts
df5 <- kelp_data_all %>% 
  filter(region == "Central_Coast") %>% 
  group_by(year, mpa_status) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area)) 

ggplot(df5, aes(x = year, y = mean_kelp_area, group = mpa_status)) + 
  geom_line(aes(color=mpa_status))+
  geom_point(aes(color=mpa_status)) + 
  theme_bw() 

ggplot(df5, aes(x = year, y = median_kelp_area, group = mpa_status)) + 
  geom_line(aes(color=mpa_status))+
  geom_point(aes(color=mpa_status)) + 
  theme_bw() 
