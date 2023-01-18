# Date: December 5th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Create GAM model of kelp area from 1984 to 2015 w/ spatial and temporal autocorrelation
# BIO 202: Ecological Statistics

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
nrow(kelp_data[kelp_data$area == 0, ])/nrow(kelp_data) # 28.2% of the kelp area are zero's

# Scale the dataset except for area 
kelp_data <- kelp_data %>% 
  mutate(new_area = area,
         lat = scale(lat, center = TRUE, scale = FALSE),
         area = scale(area),
         hsmax = scale(hsmax),
         temperature = scale(temperature),
         MHW_intensity = scale(MHW_intensity),
         CS_intensity = scale(CS_intensity),
         depth = scale(depth), 
         gravity = scale(gravity))

##### Select Datasets ##########################################################
# Dataset we will use to model with log gaussian GAMM
data <- kelp_data %>% 
  filter(log_area != 0) %>% 
  mutate(log_area = scale(log_area)) %>% 
  filter(year <= 2015) %>% 
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  drop_na() # some cells with area and temperature are missing

# Dataset we will use to model kelp presence absence
data_binary <- kelp_data %>%
  filter(year <= 2015) %>% 
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  drop_na() %>%  # some cells with area and temperature are missing
  mutate(kelp_present = ifelse(new_area == 0, 0, 1))

data_after_2015 <- kelp_data %>% 
  filter(log_area != 0) %>% 
  mutate(log_area = scale(log_area)) %>% 
  filter(year > 2015) %>% 
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  drop_na() 
  
data_binary_2015 <- kelp_data %>%
  filter(year > 2015) %>% 
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  drop_na() %>%  # some cells with area and temperature are missing
  mutate(kelp_present = ifelse(new_area == 0, 0, 1))

data_binary_2021 <- kelp_data %>%
  dplyr::select(-c(nitrate, 
                   biomass, 
                   Mpa_ID)) %>% 
  drop_na() %>%  # some cells with area and temperature are missing
  mutate(kelp_present = ifelse(new_area == 0, 0, 1))


##### GAMM Kelp Area ###########################################################
m1 <- mgcv::gamm(
  area ~
    s(lat, k = 50, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast) +
    s(depth) +
  temperature + hsmax + MHW_intensity + CS_intensity +
    mpa_status + gravity, 
  data = data,
  link = "log", 
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m1$gam)
plot(m1$gam, all.terms = T)
summary(m1$lme)
gam.check(m1$gam)

m2 <- mgcv::gamm(
  log_area ~
    s(lat, k = 50, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast) +
    s(depth) +
    temperature + hsmax + MHW_intensity + CS_intensity +
    mpa_status + gravity, 
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m2$gam)
plot(m2$gam, all.terms = T)
summary(m2$lme)
gam.check(m2$gam)

m3 <- mgcv::gamm(
  log_area ~
    s(lat, k = 100, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast) +
    s(depth) +
    temperature + hsmax + MHW_intensity + CS_intensity +
    mpa_status + gravity, 
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m3$gam)
plot(m3$gam, all.terms = T)
summary(m3$lme)
gam.check(m3$gam)

# Adjusted K value to 20
m4 <- mgcv::gamm(
  log_area ~
    s(lat, k = 100, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast, k = 20) +
    s(depth, k = 20) +
    temperature + hsmax + MHW_intensity + CS_intensity +
    mpa_status + gravity, 
  data = data,
  correlation = corAR1(form = ~ year | PixelID)
)

summary(m4$gam)
plot(m4$gam, all.terms = T)
summary(m4$lme)
gam.check(m4$gam)


##### Binomial GAMM for Kelp presence ##########################################
# adjust range parameter
model_binom1 <- mgcv::gamm(
  kelp_present ~ 
    s(lat, k = 100, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast) +
    s(depth) +
    temperature + hsmax + MHW_intensity + CS_intensity + gravity +
    mpa_status,
  data = data_binary,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)


summary(model_binom1$gam)
summary(model_binom1$lme)
plot(model_binom1$gam)
gam.check(model_binom1$gam)


model_binom2 <- mgcv::gamm(
  kelp_present ~ 
    s(lat, k = 100, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast, k = 20) +
    s(depth, k = 20) +
    temperature + hsmax + MHW_intensity + CS_intensity + 
    gravity + mpa_status,
  data = data_binary,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)


summary(model_binom2$gam)
summary(model_binom2$lme)
plot(model_binom2$gam)
gam.check(model_binom2$gam)



##### GAMMS post 2015 ##########################################################
mA <- mgcv::gamm(
  log_area ~
    s(lat, k = 100, bs = "gp", m = c(3, .03)) +
    s(distance_to_coast, k = 20) +
    s(depth, k = 20) +
    temperature + hsmax + MHW_intensity + CS_intensity +
    mpa_status + gravity, 
  data = data_after_2015,
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
    s(depth, k = 20) +
    temperature + hsmax + MHW_intensity + CS_intensity + gravity +
    mpa_status,
  data = data_binary_2015,
  family = binomial(link = "logit"),
  correlation = corAR1(form = ~ year | PixelID)
)

summary(model_binomA$gam)
summary(model_binomA$lme)
plot(model_binomA$gam, all.terms = T)
gam.check(model_binomA$gam)

##### Figures ##################################################################
library(ggplot2)


# Figure 1 heat wave plot
heatwave_years <- c(2014, 2015)
data_binary_2021$mpa_status <- factor(data_binary_2021$mpa_status, levels = c("None", "Partial", "Full"))
data_binary_2021$heat_wave <- ifelse(data_binary_2021$year %in% heatwave_years, 2, 1)
data_binary_2021$kelp_present <- factor(data_binary_2021$kelp_present)

data_binary_2021 %>% 
  ggplot(aes(x = year, y = MHW_intensity, group = year)) +
  geom_point(alpha = 0.01, color = data_binary_2021$heat_wave) +
  geom_boxplot(outlier.alpha = 0) +
  labs(y = "Cummulative Marine Heatwave Intensity", x = "Year") +
  theme_bw()

ggsave(last_plot(), filename = "Figures/Marine_heat_waves_per_year.png",
       units = "in",
       width = 6,
       height = 4,
       dpi = 600)

# Kelp histogram 
data_binary_2021 %>% 
  ggplot(aes(x = log_area, fill = kelp_present)) +
  geom_histogram(bins = 15, color = "black") +
  scale_fill_manual(values = c("#D95F02", "#1B9E77")) +
  labs(y = "Pixel Count", x = "Kelp area (log + 1)") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(last_plot(), filename = "Figures/Kelp_area_histogram.png",
       units = "in",
       width = 3,
       height = 3,
       dpi = 600)

# Picture of the dataframe
data_binary_2021 %>% select(-mpa_area, -gravity, -heat_wave) %>% 
  relocate(c(long, lat, depth), .after = PixelID) %>% 
  sample_n(6) %>% 
  view()

##### Visualize the effects of the Modles

# Log-normal GAMM up to 2015
library(mgcViz)
t <-getViz(m4$gam)

png(filename = "Figures/Log-normal_GAMM_1984_2015_effects.png", 
    res = 600,
    units ="in", 
    width = 8, 
    height = 6)

print(plot(t, allTerms = T), page = T)   
dev.off()

png(filename = "Figures/Log_normal_GAMM_1984_2015_MPAs.png", 
    res = 600,
    units ="in", 
    width = 5, 
    height = 4)

plot(t, select = 8)  
dev.off()

check(t,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))


# Binomial GAMM up to 2015
b <-getViz(model_binom2$gam)

png(filename = "Figures/Binomial_GAMM_1984_2015_effects.png", 
    res = 600,
    units ="in", 
    width = 8, 
    height = 6)

print(plot(b, allTerms = T), page = T)   
dev.off()

png(filename = "Figures/Binomial_GAMM_1984_2015_MPAs.png", 
    res = 600,
    units ="in", 
    width = 5, 
    height = 4)

plot(b, select = 7)  
dev.off()


check(b,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

# Log-normal GAMM from 2016 to 2021

t2 <-getViz(mA$gam)
print(plot(t2, allTerms = T), pages = 1)   
plot(t2, select = 8)  

png(filename = "Figures/Log_normal_GAMM_2016_2021_effects.png", 
    res = 600,
    units ="in", 
    width = 8, 
    height = 6)

print(plot(t2, allTerms = T), page = T)   
dev.off()


png(filename = "Figures/Log_normal_GAMM_2015_2021_MPAs.png", 
    res = 600,
    units ="in", 
    width = 5, 
    height = 4)

plot(t2, select = 8)  
dev.off()


check(t2,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

# Binomial GAMMfrom 2016to 2021
b2 <- getViz(model_binomA$gam)
print(plot(b2, allTerms = T), pages = 1)  
plot(b2, select = 7)

png(filename = "Figures/Binomial_GAMM_2016_2021_effects.png", 
    res = 600,
    units ="in", 
    width = 8, 
    height = 6)

print(plot(b2, allTerms = T), page = T)   
dev.off()


png(filename = "Figures/Binomial_GAMM_2015_2021_MPAs.png", 
    res = 600,
    units ="in", 
    width = 5, 
    height = 4)

plot(b2, select = 7)
dev.off()

plot1a <- plot(t, select = 8) + labs(title = "Log-normal GAMM - 1984 to 2015") 
plot1b <- plot(b, select = 7) + labs(title = "Binomial GAMM - 1984 to 2015") 
plot2a <- plot(t2, select = 8) + labs(title = "Log-normal GAMM 2016 to 2021") 
plot2b <- plot(b2, select = 7) + labs(title = "Binomial GAMM -  2016 to 2021") 

library(cowplot)

plot_grid(plot1b, plot1a, plot2b, plot2a, ncol = 2, labels = TRUE)
