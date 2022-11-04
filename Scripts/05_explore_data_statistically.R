# Date: October 16th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Explore data statistically
# BIO 202: Ecological Statistics

##### Set up ###################################################################
# load packages
library(tidyverse)
library(GGally)
library(ggeffects)

# Load data 
data <- read.csv("Processed_data/data_tables/kelp_data_w_mpas.csv")

##### Explore covariates #######################################################
# Explore each variable with time
data$year <- as.factor(data$year)

## Area
# Across year 
data %>% 
  ggplot(aes(x = year, y = area)) +
  geom_violin() # Lots of low values in the data 

# Mean across year 
data %>% 
  group_by(year) %>% 
  summarize(area_m = mean(area)) %>% 
  ggplot(aes(x = year, y = area_m)) +
  geom_point()

hist(data$area) # Super non-normal distribution, tons of zero's

## Biomass 
# Across year 
data %>% 
  ggplot(aes(x = year, y = biomass)) +
  geom_violin()

data %>% 
  group_by(year) %>% 
  summarize(biomass_m = mean(biomass)) %>% 
  ggplot(aes(x = year, y = biomass_m)) +
  geom_point()

# follows area, we will need to remove biomass 

# Should we use biomass or area?

## Temperature
# Across year 
data %>% 
  ggplot(aes(x = year, y = temperature)) +
  geom_violin()

data %>% 
  group_by(year) %>% 
  summarize(temp_m = mean(temperature)) %>% 
  ggplot(aes(x = year, y = temp_m)) +
  geom_point()

hist(data$temperature) # Not a normal distribution either (?)

## hsmax
data %>% 
  ggplot(aes(x = year, y = hsmax)) +
  geom_violin()

hist(data$hsmax) # non-normal distribution (right side tail)
hist(sqrt(data$hsmax)) # Square root transformation does a good job of normalizing the data 


##### Explore correlations between the variables ###############################

data2 <- data %>% dplyr::select(year, area, biomass, nitrate, temperature, hsmax, depth)
GGally::ggpairs(data2)

# area and biomass are highly correlated (0.665)
# Nitrate and temperature are super highly correlated, use temperature 

# Explore correlation between temperature and nitrate 
data %>% 
  filter(year == 2010) %>% 
  ggplot(aes(x = temperature, y  = nitrate)) +
  geom_point()


# Shorten data to just variables we can use
data_short <- data %>% select(year, area, temperature, hsmax, depth, mpa_status) 

data %>% 
  ggplot(aes(x = year, y = temperature)) +
  geom_point()
