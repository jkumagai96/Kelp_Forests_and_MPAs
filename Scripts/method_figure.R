# Date: June 8th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Create method graph showing how percent recovery was created
# BIO 202: Ecological Statistics

# Years: 2000 - 2021
# Max: 


##### Set up ###################################################################
# Load packages
library(tidyverse)
library(scales)

# Load Data
all_kelp_data <- read.csv("Processed_data/data_tables/kelp_data_w_mpas_per_quarter.csv")
distances <- read.csv("Processed_data/distances_to_coast.csv") %>% 
  dplyr::select(-depth) %>% 
  rename(long = x, lat = y)

##### Creating mock figure ################################################
# Formatting data
test <- all_kelp_data %>% 
  filter(PixelID == 843) %>% 
  select(year, area, PixelID, region) %>% 
  arrange(year) %>% 
  group_by(year) %>% 
  filter(year >= 2000) %>% 
  mutate(yearly_max = max(area, na.rm = T)) %>% 
  mutate(yearly_min = min(area, na.rm = T))
  
test$is_max_v <- "no"
test$is_max_v[test$area == test$yearly_max] <- "yes"

test$is_max_v[test$year >= 2014] <- "no"
test$is_max_v[test$year >= 2014 & 
                test$year < 2017 & 
                test$area == test$yearly_min] <- "min"

test$is_max_v[test$year > 2016 &
              test$area == test$yearly_max] <- "pp"

test$is_max_v[test$year == 2021] <- "pp"

t <- test %>% 
  filter(is_max_v == "yes") 

v <- mean(t$yearly_max)

group.colors <- c(no = "black", yes = "blue", min = "red", pp = "#FFBA00")

# Create figure
plot1 <- ggplot(data = test, aes(x = year, y = area/1000000, color = is_max_v)) +
  geom_point(size = 2) +
  geom_segment(aes(x = 2000, y = v/1000000, xend = 2013, yend = v/1000000), color = "blue", linewidth = 1.5)+
  geom_segment(aes(x = 2014, y = 0, xend = 2016, yend = 0),color = "red", linewidth = 1.5) +
  scale_color_manual(values=group.colors) + 
  ylab(bquote('Area '(km^2))) +
  labs(x = "Year") +
  theme_classic() +
  theme(legend.position = "none")

plot2 <- ggplot(data = test, aes(x = year, y = area/1000000)) +
  geom_point(size = 2) +
  ylab(bquote('Area '(km^2))) +
  labs(x = "Year") +
  theme_classic() +
  theme(legend.position = "none")
plot2

# Export Figure
png("Figures/Methods_figure_percent_recovery.png", width = 6, height = 4, units = "in", res = 600)
plot1
dev.off()

png("Figures/Methods_figure_percent_recovery_2.png", width = 6, height = 4, units = "in", res = 600)
plot2
dev.off()
