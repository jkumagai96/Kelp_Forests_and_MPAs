# Date: November 8th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Create model of kelp area from 1984 to 2022
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(cowplot)
library(nlme)

# Load Data
all_kelp_data <- read.csv("Processed_data/data_tables/kelp_data_w_mpas_per_quarter.csv")

###### Data Manipulation and formating #########################################
# Summarize by year
kelp_data_yr <- all_kelp_data %>% 
  group_by(PixelID, year) %>% 
  summarize(area = mean(area, na.rm = T),
         biomass = mean(biomass, na.rm = T),
         hsmax = mean(hsmax, na.rm = T),
         nitrate = mean(nitrate, na.rm = T),
         temperature = mean(temperature, na.rm = T)) 

# Subset the position details from the original data (mpa status, lat, long, Pixel ID, human gravity, and depth)
subset <- all_kelp_data %>% 
  dplyr::select(long, lat, PixelID, depth, gravity, Mpa_ID, mpa_status) %>% 
  slice_head(n = length(unique(all_kelp_data$PixelID)))

# Add this information back in
kelp_data_yr <- left_join(kelp_data_yr, subset, by = "PixelID") %>% 
  ungroup()

# Transform area
kelp_data_yr$log_area <- log(kelp_data_yr$area + 1)

# Create dataset with MPA points removed 
kelp_data_unprotected <- kelp_data_yr %>% 
  filter(mpa_status == "None")

# What percent of the data are points in Mpas? 
1 - nrow(kelp_data_unprotected)/nrow(kelp_data_yr) # 15.57%, Full is 10.66%, 4.91% is partial 
# 84.42% of the data are not within Mpas. 


##### Create figures ###########################################################
# create variable called marine heat wave (mhw) to identify years with 
# marine heat waves in the next graphs
mhw <- rep(0, 38)
mhw[c(14,15, 31, 32, 33)] <- 1

# Number of zero's per time 
plot1 <- kelp_data_yr %>% 
  group_by(year) %>% 
  summarize(n_zeros = sum(area == 0, na.rm = T)) %>% 
  mutate(heatwave = as.factor(mhw)) %>% 
  ggplot(aes(x = year, y = n_zeros, fill = heatwave, color = heatwave)) +
    geom_hline(yintercept = 397) +
    geom_hline(yintercept = 793) +
    geom_col() + 
    theme_bw() +
    theme(legend.position = "none") +
    labs(y = "Number of zeros", x = "Year") +
    scale_color_manual(values=c("#0F4392", "#DD1717")) +
    scale_fill_manual(values=c("#0F4392", "#DD1717"))

# Total area per time 
plot2 <- kelp_data_yr %>% 
  group_by(year) %>% 
  summarize(all_area = sum(area, na.rm = T)/1000000) %>% 
  mutate(heatwave = as.factor(mhw)) %>% 
  ggplot(aes(x = year, y = all_area, fill = heatwave, color = heatwave)) +
  geom_col() + 
  theme_bw() +
  theme(legend.position = c(0.85, 0.9)) +
  labs(y = "Kelp area (square kilometers)",
       x = "Year") +
  scale_color_manual(values=c("#0F4392", "#DD1717")) +
  scale_fill_manual(values=c("#0F4392", "#DD1717"))

plot_grid(plot1, plot2, labels = "auto")

# Average kelp area within full MPAs, partial MPAs, and no MPAs in general, from 2012 to 2021? 
# 2012 is the latest year of establishment where 25 MPAs were established!

# Line Graph - Proportion of Zeros
plot3 <- kelp_data_yr %>% 
  filter(year >= 2012) %>% 
  group_by(year, mpa_status) %>% 
  summarize(n_zeros = sum(area == 0, na.rm = T),
            n_values = sum(area != 0, na.rm = T)) %>% 
  mutate(p_zeros = n_zeros/(n_zeros + n_values)) %>% 
  ggplot(aes(x = year, y = p_zeros, color = mpa_status)) + 
  geom_line(size = 1.5) + 
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Proporation of zeros", x = "Year")

# We see an increase in zero's during the marine heat wave for areas in partial mpas
# and not in MPAs. Between 2014 and 2016 we also see this for full MPAs, but less severe 
# Zero's seem to increase in time from 2014 to 2021.
# There's a higher proportion of zero's outside of MPAs and the samllest proportion 
# of zero's in fully protected MPAs

# Line graph - mean area 
plot4 <- kelp_data_yr %>% 
  filter(year >= 2012) %>% 
  group_by(year, mpa_status) %>% 
  summarize(mean_area = mean(area, na.rm = T)) %>% 
  ggplot(aes(x = year, y = mean_area, color = mpa_status)) + 
  geom_line(size = 1.5) + 
  theme_bw() +
  theme(legend.position = c(0.8, 0.85)) +
  labs(y = "Mean kelp area (m2)", x = "Year")

# Line graph - median area
plot4_median <- kelp_data_yr %>% 
  filter(year >= 2012) %>% 
  group_by(year, mpa_status) %>% 
  summarize(median_area = median(area, na.rm = T)) %>% 
  ggplot(aes(x = year, y = median_area, color = mpa_status)) + 
  geom_line(size = 1.5) + 
  theme_bw() +
  theme(legend.position = c(0.8, 0.8)) +
  labs(y = "Median kelp area (m2)", x = "Year")

# Less area for kelp with no protection, more area on average for paritatially protected MPAs
# which is interesting. 
# Area drastically decreases in 2014 for all categories, partially protected MPAs rebound quickly 

plot_grid(plot3, plot4, labels = "auto")
plot_grid(plot4, plot4_median, labels = "auto")

##### Figures - Box plots per year #############################################
kelp_data_yr$mpa_status <- factor(kelp_data_yr$mpa_status, levels = c("None", "Partial", "Full"))
kelp_data_yr <- kelp_data_yr %>% 
  mutate(hw = if_else(year > 2013, "during", "after")) %>% 
  mutate(hw = if_else(year > 2016, "after", hw)) %>% 
  mutate(hw = if_else(year < 2014, "before", hw)) 

# My favorite plot
plot5 <- kelp_data_yr %>% 
  filter(year >= 2012) %>% 
  group_by(year, mpa_status, hw) %>% 
  summarize(mean_area = mean(area, na.rm = T)) %>%
  ggplot(aes(x = mpa_status, y = mean_area, color = factor(hw, levels = c("before", "during", "after")))) + 
  geom_boxplot(lwd = 1) + 
  theme_bw() +
  scale_color_manual(values = c("#440154FF","#FE6E00", "#2A788EFF")) +
  labs(y = "Mean Kelp Area (m2)", x = "MPA Category", color = "Heatwave") 

##### Models of kelp area without points in  MPAs ##############################
data <- kelp_data_unprotected %>% 
  select(-Mpa_ID, -mpa_status) %>% 
  drop_na()

## Linear Model
Mlm1 <- lm(log_area ~ temperature + hsmax + depth + gravity,
             data = data)

## Mixed effects linear model 
Mlme1 <- lme(log_area ~ temperature + hsmax + depth + gravity,
             random = ~1 + temperature|PixelID,
             method = "ML",
             data = data)

Mlme2 <- lme(log_area ~ temperature + hsmax + gravity,
             random = ~1 + temperature|PixelID,
             method = "ML",
             data = data)

Mlme3 <- lme(log_area ~ temperature + hsmax + depth,
             random = ~1 + temperature|PixelID,
             method = "ML",
             data = data)

summary(Mlme1)
summary(Mlme2)
summary(Mlme3)

AIC(Mlme1, Mlme2, Mlme3) # Mlme3 is better

Mlme3a <- lme(log_area ~ temperature + hsmax + depth,
             random = ~1 + hsmax|PixelID,
             method = "ML",
             data = data)

Mlme3b <- lme(log_area ~ temperature + hsmax + depth,
              random = ~1 + depth|PixelID,
              method = "ML",
              data = data)

summary(Mlme3a) # all terms are significant
summary(Mlme3b) # all terms are significant

AIC(Mlme3, Mlme3a, Mlme3b) #Mlme3a is significantly best 

plot(Mlme3a)

# Remove zeros?
data_no_zeros <- data %>% 
  filter(log_area != 0)

Mlme_new_1 <- lme(log_area ~ temperature + hsmax + depth + gravity,
              random = ~1 + hsmax|PixelID,
              method = "ML",
              data = data_no_zeros)
summary(Mlme_new_1)

Mlme_new_2 <- lme(log_area ~ temperature + hsmax + depth,
                  random = ~1 + hsmax|PixelID,
                  method = "ML",
                  data = data_no_zeros)

summary(Mlme_new_2) # all significant 

AIC(Mlme_new_1, Mlme_new_2) # Mlme_new_2 slightly better

Mlme_new_2a <- lme(log_area ~ temperature + hsmax + depth,
                  random = ~1 + depth|PixelID,
                  method = "ML",
                  data = data_no_zeros)

AIC(Mlme_new_2, Mlme_new_2a)

# Fit the best model with REML
Mlme_no_zeros_best <- lme(log_area ~ temperature + hsmax + depth,
                  random = ~1 + hsmax|PixelID,
                  method = "REML",
                  data = data_no_zeros)
summary(Mlme_no_zeros_best) # all significant
plot(Mlme_no_zeros_best) # residuals are still skewed 
residuals <- resid(Mlme_no_zeros_best)
hist(residuals)

# Try to fit generalized linear mixed model?
library(lme4)
data$presence <- ifelse(data$area == 0, 0, 1)

MGlmm1 <- glmer(presence ~ temperature + hsmax + depth + gravity + 
                  (1 | PixelID),
                family = "binomial", 
                data = data) 


MGlmm2 <- glmer(log_area ~ temperature + hsmax + depth + gravity + 
                  (1 | PixelID),
                family = "gaussian", 
                data = data) 
# This fit a linear mixed effect model? 

# What other variables do you want to throw in there?
# Slope, marine heat wave days, marine cold snap days?


