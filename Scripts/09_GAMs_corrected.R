# Date: December 1st 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Create GAM model of kelp area from 1984 to 2022 w/ spatial and temporal autocorrelation
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(mgcv)
library(statmod)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_per_year.csv")

###### Data Manipulation and formating #########################################
# Remove rows that are NA for MHW and CS (eventually needs to be corrected)
kelp_data <- kelp_data_all[!is.na(kelp_data_all$MHW_intensity), ]

# Transform area
kelp_data$log_area <- log(kelp_data$area + 1)

# Transform mpa_status into factor
kelp_data$mpa_status <- factor(kelp_data$mpa_status, levels = c("None", "Partial", "Full"))

# Count percentage of zero's for NAs
nrow(kelp_data[kelp_data$area == 0, ])/nrow(kelp_data) # 27.7% of the kelp area are zero's

##### Select Datasets ##########################################################

# Dataset we will use to model with log gaussian GAM
data <- kelp_data %>% 
  # filter(year >= 2012) %>% 
  # Explore modeling kelp area with protected areas, if so you need to filter to 2012
  
  filter(mpa_status == "None") %>% 
  filter(log_area != 0) %>% 
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  drop_na() # some cells with area and temperature are missing

data_in_mpas <- kelp_data %>% 
  # filter(year >= 2012) %>% 
  # Explore modeling kelp area with protected areas, if so you need to filter to 2012
  
  filter(mpa_status != "None") %>% 
  # filter(log_area != 0) %>% # we want all the data within mpas
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  drop_na() # some cells with area and temperature are missing


# Dataset we will use to model kelp presence absence
data_binary <- kelp_data %>%
  filter(mpa_status == "None") %>% 
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  drop_na() %>%  # some cells with area and temperature are missing
  mutate(kelp_present = ifelse(area == 0, 0, 1))

##### Run GAMs for Log Gaussian ################################################
## Returns two objects - a GAM and a nlme object with AR1 fit to residuals by pixel id
m1 <- mgcv::gamm(
  log_area ~ s(lat, k = 50, bs = "gp", m = c(3, .1)) +
    temperature + hsmax + depth + gravity + MHW_intensity + CS_intensity, 
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

## Significance of parametric (e.g. intercept, linear terms) and smooth terms
summary(m1$gam)
plot(m1$gam)

## "Phi" is the avg. correlation of pixel value at t and t - 1 for all t
summary(m1$lme)

# Delete cs_intensity
m2 <- mgcv::gamm(
  log_area ~ s(lat, k = 50, bs = "gp", m = c(3, .1)) +
    temperature + hsmax + depth + gravity + MHW_intensity, 
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m2$gam)
plot(m2$gam)
summary(m2$lme)

# Delete gravity instead of cs_intensity
m2.5 <- mgcv::gamm(
  log_area ~ s(lat, k = 50, bs = "gp", m = c(3, .1)) +
    temperature + hsmax + depth + MHW_intensity + CS_intensity, 
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m2.5$gam)
plot(m2.5$gam)
summary(m2.5$lme)

# Delete gravity
m3 <- mgcv::gamm(
  log_area ~ s(lat, k = 50, bs = "gp", m = c(3, .1)) +
    temperature + hsmax + depth + MHW_intensity, 
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m3$gam)
plot(m3$gam)
summary(m3$lme)

# Try to add in long smoother too, added CS_intensity and gravity back in 
m4 <- mgcv::gamm(
  log_area ~ s(lat, k = 50, bs = "gp", m = c(3, .1)) +
    s(long, k = 50, bs = "gp", m = c(3,.1)) +
    temperature + hsmax + MHW_intensity + CS_intensity + depth + gravity, 
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m4$gam)
plot(m4$gam)
summary(m4$lme)

# Remove gravity
m5 <- mgcv::gamm(
  log_area ~ s(lat, k = 50, bs = "gp", m = c(3, .1)) +
    s(long, k = 50, bs = "gp", m = c(3,.1)) +
    temperature + hsmax + MHW_intensity + CS_intensity + depth, 
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m5$gam)
plot(m5$gam)
summary(m5$lme)

# All variables significant, add smoother to temperature 
m6 <- mgcv::gamm(
  log_area ~ s(lat, k = 50, bs = "gp", m = c(3, .1)) +
    s(long, k = 50, bs = "gp", m = c(3,.1)) +
    s(temperature) + 
    hsmax + MHW_intensity + CS_intensity + depth, 
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m6$gam)
plot(m6$gam)
summary(m6$lme)

# A smoother on temperature did not increase R2 by much, and it is very linear 
# so I will not add this to the model 

# Fix the range parameter (based on outcomes lower in the script)
mC <- mgcv::gamm(
  log_area ~ s(lat, k = 50, bs = "gp", m = c(3, .02)) +
    s(long, k = 50, bs = "gp", m = c(3,.02)) +
    temperature + hsmax + MHW_intensity + CS_intensity + depth + mpa_status, 
  data = data_decade,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(mC$gam)
summary(mC$lme)
# Model mC is the best model with mpa_status! 

##### Adjust the range parameter ###############################################
m5_a <- mgcv::gamm(
  log_area ~ s(lat, k = 50, bs = "gp", m = c(3, .05)) +
    s(long, k = 50, bs = "gp", m = c(3,.05)) +
    temperature + hsmax + MHW_intensity + CS_intensity + depth + gravity, 
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m5_a$gam)
plot(m5_a$gam)
summary(m5_a$lme)

m5_b <- mgcv::gamm(
  log_area ~ s(lat, k = 50, bs = "gp", m = c(3, .03)) +
    s(long, k = 50, bs = "gp", m = c(3,.03)) +
    temperature + hsmax + MHW_intensity + CS_intensity + depth + gravity,
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m5_b$gam)
plot(m5_b$gam)
summary(m5_b$lme)

# BEST MODEL SO FAR! :) Adjusted k parameter as well
m5_c <- mgcv::gamm(
  log_area ~ s(lat, k = 100, bs = "gp", m = c(3, .02)) +
    s(long, k = 100, bs = "gp", m = c(3,.02)) +
    temperature + hsmax + MHW_intensity + CS_intensity + depth + gravity, 
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m5_c$gam)
plot(m5_c$gam)
summary(m5_c$lme)

# Adding temperature smoother doesn't make sense, only a .1% increase in predictive power 

##### Run Gamma GAMMs ##########################################################
# Use kelp area instead of log kelp area, that doesn't work... so switch to log_area
# still doesn't work! Adjusting range parameter from .02 to .1 does not fix it 
m1_gamma <- mgcv::gamm(
  log_area ~ s(lat, k = 20, bs = "gp", m = c(3, .1)) + 
    s(long, k = 20, bs = "gp", m = c(3,.1)) + 
    temperature + hsmax + depth + gravity + MHW_intensity + CS_intensity,
  data = data,
  family = Gamma, 
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m1_gamma$gam)
plot(m1_gamma$gam)
summary(m1_gammalme)

##### Run Binomial GAMs ########################################################
model_binom1 <- mgcv::gamm(
  kelp_present ~ s(lat, k = 20, bs = "gp", m = c(3, .1)) + 
    s(long, k = 20, bs = "gp", m = c(3,.1)) + 
    temperature + hsmax + depth + gravity + MHW_intensity + CS_intensity,
  data = data_binary,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)

summary(model_binom1$gam)
summary(model_binom1$lme)
plot(model_binom1$gam)

# Delete gravity
model_binom2 <- mgcv::gamm(
  kelp_present ~ s(lat, k = 20, bs = "gp", m = c(3, .1)) + 
    s(long, k = 20, bs = "gp", m = c(3,.1)) + 
    temperature + hsmax + depth + MHW_intensity + CS_intensity,
  data = data_binary,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)

summary(model_binom2$gam)
summary(model_binom2$lme)
plot(model_binom2$gam)

# Delete temperature
model_binom3 <- mgcv::gamm(
  kelp_present ~ s(lat, k = 20, bs = "gp", m = c(3, .1)) + 
    s(long, k = 20, bs = "gp", m = c(3,.1)) + 
    hsmax + depth + MHW_intensity + CS_intensity,
  data = data_binary,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)

summary(model_binom3$gam)
summary(model_binom3$lme)
plot(model_binom3$gam)

# Adjust range parameter  (not so good)
model_binom4 <- mgcv::gamm(
  kelp_present ~ s(lat, k = 20, bs = "gp", m = c(3, .02)) + 
    s(long, k = 20, bs = "gp", m = c(3,.02)) + 
    hsmax + depth + MHW_intensity + CS_intensity,
  data = data_binary,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)

summary(model_binom4$gam)
summary(model_binom4$lme)
plot(model_binom4$gam)

# Adjust k parameter, this signficantly betters the model! R2 is now 0.186
model_binom5 <- mgcv::gamm(
  kelp_present ~ s(lat, k = 50, bs = "gp", m = c(3, .1)) + 
    s(long, k = 50, bs = "gp", m = c(3,.1)) + 
    hsmax + depth + MHW_intensity + CS_intensity,
  data = data_binary,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)

# adjust range parameter
model_binom6 <- mgcv::gamm(
  kelp_present ~ s(lat, k = 50, bs = "gp", m = c(3, .5)) + 
    s(long, k = 50, bs = "gp", m = c(3,.5)) + 
    hsmax + depth + MHW_intensity + CS_intensity,
  data = data_binary,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)


summary(model_binom6$gam)
summary(model_binom6$lme)
plot(model_binom6$gam)

# adjust range parameter
model_binom7 <- mgcv::gamm(
  kelp_present ~ s(lat, k = 50, bs = "gp", m = c(3, .05)) + 
    s(long, k = 50, bs = "gp", m = c(3,.05)) + 
    hsmax + depth + MHW_intensity + CS_intensity,
  data = data_binary,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)


summary(model_binom7$gam)
summary(model_binom7$lme)
plot(model_binom7$gam)

# adjust range parameter
model_binom8 <- mgcv::gamm(
  kelp_present ~ s(lat, k = 50, bs = "gp", m = c(3, .03)) + 
    s(long, k = 50, bs = "gp", m = c(3,.03)) + 
    hsmax + depth + MHW_intensity + CS_intensity,
  data = data_binary,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)


summary(model_binom8$gam)
summary(model_binom8$lme)
plot(model_binom8$gam)

# adjust range parameter
model_binom9 <- mgcv::gamm(
  kelp_present ~ s(lat, k = 50, bs = "gp", m = c(3, .02)) + 
    s(long, k = 50, bs = "gp", m = c(3,.02)) + 
    hsmax + depth + MHW_intensity + CS_intensity,
  data = data_binary,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)


summary(model_binom9$gam)
summary(model_binom9$lme)
plot(model_binom9$gam)


## model_binom8 is best!!! with 0.03 range parameter! :) gravity and temperature is removed

##### Compare Model Predictions within MPAs in figures #########################

# Obtain model predictions for combined estimates

binomial_prediction <-  predict(model_binom8$gam, newdata = data_in_mpas, type = "response") 
binary_prediction <- ifelse(binomial_prediction > .5, 1, 0)

data_in_mpas$fitted_area <- binary_prediction *
  exp(predict(m5_c$gam, newdata = data_in_mpas) + 1)


# Compare fitted vs. actual area values
data_in_mpas %>% 
  ggplot(aes(x = area, y = fitted_area)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "blue")

# group the data into categories per year (mpa status = partial and full)
colors <- c("#DC3220", "#005AB5")

### Boxplots per year showing difference between fitted and measured values

# Fully protected mpas
data_in_mpas %>% 
  rename(measured_area = area) %>% 
  select(PixelID, year, measured_area, fitted_area, mpa_status) %>% 
  filter(year >= 2012) %>% 
  filter(mpa_status == "Full") %>% 
  mutate(year = as.factor(year)) %>% 
  pivot_longer(c(measured_area,fitted_area),values_to = "kelp_area", names_to = "type") %>% 
  ggplot(aes(x = year, y = log(kelp_area + 1), color = type)) +
  geom_boxplot() +
  scale_color_manual(values = colors, name = "") +
  labs(y = "Kelp Area", x = "Year", title = "No-take MPAs") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(last_plot(), 
       filename = "Figures/comparison_model_measured_full_MPAs.png",
       width = 8,
       height = 5,
       units = "in", 
       dpi = 300)

# Partially protected mpas
data_in_mpas %>% 
  rename(measured_area = area) %>% 
  select(PixelID, year, measured_area, fitted_area, mpa_status) %>% 
  filter(year >= 2012) %>% 
  filter(mpa_status == "Partial") %>% 
  mutate(year = as.factor(year)) %>% 
  pivot_longer(c(measured_area,fitted_area),values_to = "kelp_area", names_to = "type") %>% 
  ggplot(aes(x = year, y = log(kelp_area + 1), color = type)) +
  geom_boxplot() +
  scale_color_manual(values = colors, name = "") +
  labs(y = "Kelp Area", x = "Year", title = "MPAs with Parital Protection") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(last_plot(), 
       filename = "Figures/comparison_model_measured_partial_MPAs.png",
       width = 8,
       height = 5,
       units = "in", 
       dpi = 300)


##### Run Full Models from 2012 to 2021 w/ all data ############################
# Consider how much the model is improved with mpa_status added
# need to filter data to 2012 to do this 
data_decade <- kelp_data %>% 
  filter(year >= 2012) %>% 
  filter(log_area != 0) %>% 
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  drop_na() # some cells with area and temperature are missing

# MPA status is included
mA <- mgcv::gamm(
  log_area ~ s(lat, k = 50, bs = "gp", m = c(3, .02)) +
    s(long, k = 50, bs = "gp", m = c(3,.02)) +
    temperature + hsmax + MHW_intensity + CS_intensity + depth + mpa_status, 
  data = data_decade,
  correlation = corAR1(form = ~ year | PixelID)
)
# removed gravity

summary(mA$gam)
summary(mA$lme)

# adjusted range parameter to 0.03 and 0.04, not better... 


# Adjust K parameter, big improvement!
mB <- mgcv::gamm(
  log_area ~ s(lat, k = 100, bs = "gp", m = c(3, .02)) +
    s(long, k = 100, bs = "gp", m = c(3,.02)) +
    temperature + hsmax + MHW_intensity + CS_intensity + depth + mpa_status, 
  data = data_decade,
  correlation = corAR1(form = ~ year | PixelID)
)
# removed gravity

summary(mB$gam)
summary(mB$lme)
plot(mB$gam)

# mB is best! 


# Binomial Model next!
data_decade_binomial <- kelp_data %>% 
  filter(year >= 2012) %>% 
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  mutate(kelp_present = ifelse(area == 0, 0, 1)) %>% 
  drop_na() # some cells with area and temperature are missing

model_binom_A <- mgcv::gamm(
  kelp_present ~ s(lat, k = 20, bs = "gp", m = c(3, .03)) + 
    s(long, k = 20, bs = "gp", m = c(3,.03)) + 
    temperature + hsmax + depth + gravity + MHW_intensity + CS_intensity + mpa_status,
  data = data_decade_binomial,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)

summary(model_binom_A$gam)
summary(model_binom_A$lme)
plot(model_binom_A$gam)

model_binom_B <- mgcv::gamm(
  kelp_present ~ s(lat, k = 20, bs = "gp", m = c(3, .03)) + 
    s(long, k = 20, bs = "gp", m = c(3,.03)) + 
    temperature + hsmax + depth + MHW_intensity + mpa_status, # removed gravity and CS_intensity
  data = data_decade_binomial,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)

summary(model_binom_B$gam)
summary(model_binom_B$lme)
plot(model_binom_B$gam)

model_binom_C <- mgcv::gamm(
  kelp_present ~ s(lat, k = 50, bs = "gp", m = c(3, .03)) + # adjusted k parameter
    s(long, k = 50, bs = "gp", m = c(3,.03)) + 
    temperature + hsmax + depth + MHW_intensity + mpa_status,
  data = data_decade_binomial,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)

summary(model_binom_C$gam)
summary(model_binom_C$lme)
plot(model_binom_C$gam)

# R2 value is .202

