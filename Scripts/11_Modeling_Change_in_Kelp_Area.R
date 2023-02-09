# Date: January 17th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Model change in kelp area from 2008 to 2021. 
# Purpose 2: Model kelp area w/ presence/absence then gaussian w/ log link function but not seperated by year 
# MPAs, heatwaves, and kelp project 

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(mgcv)
library(statmod)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")

###### Data Manipulation and formating #########################################
# Remove rows that are NA for MHW and CS (eventually needs to be corrected)
kelp_data <- kelp_data_all[!is.na(kelp_data_all$MHW_intensity), ]

# Transform area
kelp_data$log_area <- log(kelp_data$area + 1)

# Transform mpa_status into factor
kelp_data$mpa_status <- factor(kelp_data$mpa_status, levels = c("None", "Partial", "Full"))

# Count percentage of zero's for NAs
nrow(kelp_data[kelp_data$area == 0, ])/nrow(kelp_data) # 28.2% of the kelp area are zero's #*** I get 31%?

# Scale the dataset except for area 
kelp_data <- kelp_data %>% 
  mutate(area = area,
         lat = scale(lat, center = TRUE, scale = FALSE),
         area_s = scale(area),
         hsmax = scale(hsmax),
         temperature = scale(temperature),
         MHW_intensity = scale(MHW_intensity),
         CS_intensity = scale(CS_intensity),
         depth = scale(depth), 
         gravity = scale(gravity))

###### Pre-processing ##########################################################
# Calculate mean area per pixel between 1984 to 2007
mean_kelp_area <- kelp_data %>% 
  filter(year <= 2007) %>% 
  group_by(PixelID) %>% 
  summarise(mean_area = mean(area, na.rm = T))

# Join the longstanding average kelp area per pixel, 
# then calculate the change in kelp area per year for 2008 to 2021 (6 years before, 6 years after)

kelp_data_2008 <- kelp_data %>% 
  filter(year > 2007) %>% 
  left_join(mean_kelp_area, by = "PixelID") %>% 
  mutate(change_in_kelp = area - mean_area) 

hist(kelp_data_2008$change_in_kelp)

# Add in factor of before, during, after marine heat wave event
kelp_data_2008 <- kelp_data_2008 %>% 
  mutate(mhw_event = ifelse(year < 2014, "before", "during")) %>% 
  mutate(mhw_event = ifelse(year > 2015, "after", mhw_event)) %>% 
  mutate(mhw_event = factor(mhw_event, levels = c("during", "before", "after"))) 
           ifelse(year < 2014, "before", "during") %>% 
  mutate(mhw_event = ifelse(year > 2015, "after", mhw_event)) %>% 
  mutate(mhw_event = factor(mhw_event, levels = c("during", "before", "after"))) 
#*** Note that case_when, might be easier here? MUCH shorter...
kelp_data_2008 <- kelp_data_2008 %>% 
 mutate(mhw_event = case_when(
   year < 2014 ~ "before",
   year >= 2014 & year < 2015 ~ "during",
   TRUE ~ "after" # Everything that is not before or during is after
   ),
   mhw_event = factor(mhw_event, levels = c("during", "before", "after")))           


##### GAMM Change in Kelp Area #################################################
m1 <- mgcv::gamm(
  change_in_kelp ~
    s(lat, k = 50, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast) +
    s(depth) +
    temperature + hsmax + mhw_event + MHW_intensity + CS_intensity + mpa_status + gravity + mpa_area +
    mpa_status*gravity, 
  data = kelp_data_2008,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m1$gam)
plot(m1$gam, all.terms = T)
summary(m1$lme)
gam.check(m1$gam)


# Removed gravity, added the interaction betweeen mpa status and mhw event 
m2 <- mgcv::gamm(
  change_in_kelp ~
    s(lat, k = 50, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast, k = 50) +
    s(depth) +
    temperature + hsmax + mhw_event + MHW_intensity + CS_intensity + mpa_status + mpa_area +
    mpa_status*mhw_event, 
  data = kelp_data_2008,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m2$gam)

# These models don't really make sense! 

##### Select Datasets for split GAMMs modeling kelp area #######################
# Dataset we will use to model with log gaussian GAMM
data <- kelp_data %>% 
  filter(log_area != 0) %>% 
  mutate(log_area = scale(log_area)) %>% 
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  drop_na()  %>% # some cells with area and temperature are missing
  mutate(mhw_event = ifelse(year < 2014, "before", "during")) %>% 
  mutate(mhw_event = ifelse(year > 2015, "after", mhw_event)) %>% 
  mutate(mhw_event = factor(mhw_event, levels = c("during", "before", "after"))) 


data_from_2008 <- kelp_data %>% 
  filter(log_area != 0) %>% 
  filter(year >= 2008) %>% 
  mutate(log_area = scale(log_area)) %>% 
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  drop_na()  %>% # some cells with area and temperature are missing
  mutate(mhw_event = ifelse(year < 2014, "before", "during")) %>% 
  mutate(mhw_event = ifelse(year > 2015, "after", mhw_event)) %>% 
  mutate(mhw_event = factor(mhw_event, levels = c("during", "before", "after"))) 


# Dataset we will use to model kelp presence absence
data_binary <- kelp_data %>%
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  drop_na() %>%  # some cells with area and temperature are missing
  mutate(kelp_present = ifelse(area == 0, 0, 1)) %>% 
  mutate(mhw_event = ifelse(year < 2014, "before", "during")) %>% 
  mutate(mhw_event = ifelse(year > 2015, "after", mhw_event)) %>% 
  mutate(mhw_event = factor(mhw_event, levels = c("during", "before", "after"))) 

data_binary_from_2008 <- kelp_data %>%
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  filter(year >= 2008) %>% 
  drop_na() %>%  # some cells with area and temperature are missing
  mutate(kelp_present = ifelse(area == 0, 0, 1)) %>% 
  mutate(mhw_event = ifelse(year < 2014, "before", "during")) %>% 
  mutate(mhw_event = ifelse(year > 2015, "after", mhw_event)) %>% 
  mutate(mhw_event = factor(mhw_event, levels = c("during", "before", "after"))) 

##### GAMMS ####################################################################
#*** Just tinkering
m1 <- gam(change_in_kelp ~
            s(year) + s(MHW_intensity) + s(mpa_area) + s(nitrate) + s(gravity) + s(depth) + mpa_status +
            s(PixelID, bs = "re") +
            s(PixelID, year, bs = "re"),
           data = kelp_data_2008)
summary(m1)
plot(m1)




mA <- mgcv::gamm(
  area_s ~
    s(lat, k = 100, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast, k = 20) +
    s(depth) +
    temperature + hsmax + MHW_intensity + CS_intensity +
    mpa_status + gravity + mpa_status*mhw_event + mhw_event + mpa_area, 
  data = data,
  link = "log",
  correlation = corAR1(form = ~ year | PixelID)
)

summary(mA$gam)
plot(mA$gam, all.terms = T)
summary(mA$lme)
gam.check(mA$gam)


model_binomA <- mgcv::gamm(
  kelp_present ~ 
    s(lat, k = 100, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast, k = 20) +
    s(depth) +
    temperature + hsmax + MHW_intensity + CS_intensity +
    mpa_status + gravity + mpa_status*mhw_event + mhw_event + mpa_area, 
  data = data_binary,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)

summary(model_binomA$gam)
summary(model_binomA$lme)
plot(model_binomA$gam, all.terms = T)
gam.check(model_binomA$gam)


##### GAMMS considering just 2008 and after ####################################
mB <- mgcv::gamm(
  area_s ~
    s(lat, k = 100, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast, k = 20) +
    s(depth, k = 5) +
    temperature + hsmax + MHW_intensity + CS_intensity +
    mpa_status + gravity + mpa_status*mhw_event + mhw_event + mpa_area, 
  data = data_from_2008,
  link = "log",
  correlation = corAR1(form = ~ year | PixelID)
)

summary(mB$gam)

model_binomB <- mgcv::gamm(
  kelp_present ~ 
    s(lat, k = 100, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast, k = 20) +
    s(depth, k = 5) +
    temperature + hsmax + MHW_intensity + CS_intensity +
    mpa_status + gravity + mpa_status*mhw_event + mhw_event + mpa_area, 
  data = data_binary_from_2008,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)

summary(model_binomB$gam)
