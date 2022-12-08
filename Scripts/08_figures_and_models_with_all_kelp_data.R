# Date: November 8th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Create model of kelp area from 1984 to 2022
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(cowplot)
library(nlme)
library(mgcv)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")

###### Data Manipulation and formating #########################################
# Remove rows that are NA for MHW and CS (eventually needs to be corrected)
kelp_data <- kelp_data_all[!is.na(kelp_data_all$MHW_intensity), ]

# Transform area
kelp_data$log_area <- log(kelp_data$area + 1)

# Create dataset with MPA points removed 
kelp_data_unprotected <- kelp_data %>% 
  filter(mpa_status == "None")

# What percent of the data are points in Mpas? 
1 - nrow(kelp_data_unprotected)/nrow(kelp_data) 
# 93% of the data are not within Mpas. 


##### Create figures ###########################################################
# create variable called marine heat wave (mhw) to identify years with 
# marine heat waves in the next graphs
mhw <- rep(0, 38)
mhw[c(14,15, 31, 32, 33)] <- 1

# Number of zero's per time 
plot1 <- kelp_data %>% 
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

kelp_data_yr <- kelp_data

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
kelp_data$mpa_status <- factor(kelp_data$mpa_status, levels = c("None", "Partial", "Full"))
kelp_data <- kelp_data %>% 
  mutate(hw = if_else(year > 2013, "2014-2015", "2016-2021")) %>% 
  mutate(hw = if_else(year > 2015, "2016-2021", hw)) %>% 
  mutate(hw = if_else(year < 2014, "2008-2013", hw)) 

# My favorite plot
plot5 <- kelp_data %>% 
  filter(year >= 2008) %>% 
  group_by(year, mpa_status, hw) %>% 
  summarize(mean_area = mean(area, na.rm = T)) %>%
  ggplot(aes(x = mpa_status, y = mean_area, color = factor(hw, levels = c("2008-2013", "2014-2015", "2016-2021")))) + 
  geom_boxplot(lwd = 1) + 
  theme_bw() +
  #theme(legend.position="bottom") +
  scale_color_manual(values = c("#440154FF","#FE6E00", "#2A788EFF")) +
  labs(y = "Mean Kelp Area (m2)", x = "MPA Category", color = "") 

ggsave(last_plot(), filename = "Figures/Before_during_after_heatwave.png",
       dpi = 600,
       units = "in",
       height = 5,
       width = 8)

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

##### Fit Models with GLMs and GLMMs ###########################################
data <- kelp_data_unprotected %>% 
  select(-Mpa_ID, -mpa_status) %>% 
  drop_na()


glm1 <- glm(round(area) ~ temperature + hsmax + depth + gravity + MHW_intensity + CS_intensity,
            data = data,
            family = "quasipoisson")
summary(glm1)

library(MASS)

glm2 <- glm.nb(round(area) ~ temperature + hsmax + depth + gravity + MHW_intensity + CS_intensity,
            data = data)
summary(glm2)

library(lme4)
glme1 <- glmer(round(area) ~ temperature + hsmax + depth + 
                             gravity + MHW_intensity + CS_intensity +
                            (1|PixelID),
               data = data, 
               family = "quasipoisson")
summary(glme1)


# For some reason gamma does not allow zeros?
glme2 <- glmer(area ~ temperature + hsmax + depth + # failing to work
                 gravity + MHW_intensity + CS_intensity +
                 (1|PixelID),
               data = data_no_zeros, 
               family = "Gamma")

glm3 <- glm(log_area ~ temperature + hsmax + depth + 
                 gravity + MHW_intensity + CS_intensity,
               data = data_no_zeros, 
               family = "Gamma")
summary(glm3)



lme1 <- lme(log_area ~ temperature + hsmax + depth + gravity + 
              MHW_intensity + CS_intensity +
              long + lat,
    random = ~1 + hsmax|PixelID,
    method = "ML",
    data = data)

summary(lme1)

# Drop Gravity
lme2 <- lme(log_area ~ temperature + hsmax + depth + 
              MHW_intensity + CS_intensity +
              long + lat,
            random = ~1 + hsmax|PixelID,
            method = "ML",
            data = data)

summary(lme2)

lme3 <- lme(log_area ~ temperature + hsmax +
              MHW_intensity + CS_intensity +
              long + lat,
            random = ~1 + hsmax|PixelID,
            method = "ML",
            data = data)

summary(lme3)

AIC(lme1, lme2, lme3)

residuals <- resid(lme3)
hist(residuals)
plot(lme3)

###### Fit Figures without zeros
# Remove zeros
data_no_zeros <- data %>% 
  filter(log_area != 0)

lme1_nz <- lme(log_area ~ temperature + hsmax + depth + gravity + 
              MHW_intensity + CS_intensity +
              long + lat,
            random = ~1 + hsmax|PixelID,
            method = "ML",
            data = data_no_zeros)
summary(lme1_nz)

lme2_nz <- lme(log_area ~ temperature + hsmax + depth + 
                 MHW_intensity + CS_intensity +
                 long + lat,
               random = ~1 + hsmax|PixelID,
               method = "ML",
               data = data_no_zeros)
summary(lme2_nz)
plot(lme2_nz)

residuals <- resid(lme2_nz)
hist(residuals)

require(MuMIn)
r.squaredGLMM(lme2_nz)

##### Account for time #########################################################
# Need to trouble shoot this 
library(glmmTMB)

m_time1 <- glmmTMB(log_area ~ temperature + hsmax + depth + 
                    MHW_intensity + CS_intensity + long + lat +
                    ar(0 + as.factor(year)|PixelID), 
                  data = data)


##### 
