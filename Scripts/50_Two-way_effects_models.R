# Date: November 20th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Two way effects models for abundances of urchins
# BIO 202: Ecological Statistics

##### Load packages and data ###################################################
library(fixest)
library(broom)
library(modelsummary)
library(tidyverse)

# Load Data
PISCO_data <- read.csv("Data/PISCO subtidal data clearinghouse/PISCO_kelpforest_combined_data_bysite.csv")
sites_for_joining <- read.csv("Processed_data/data_tables/PISCO_sites_with_MPAs.csv") 

# Only biomass for sheephead >= 35cm https://royalsocietypublishing.org/doi/full/10.1098/rspb.2016.1936
PISCO_biomass_data <- read.csv("Processed_data/Data_tables/PISCO_biomass_large_per_site.csv")

##### Format Data ##############################################################
data <- PISCO_data %>% 
  left_join(sites_for_joining, by = "SITE") %>% 
  dplyr::select(SITE,
                LATITUDE,
                LONGITUDE,
                region,
                mpa_status,
                Estab_Yr_1,
                year,
                fish_SPUL,
                swath_PANINT,
                swath_CENCOR,
                swath_MESFRAAD,
                swath_STRPURAD
  ) %>% 
  filter(year >= 2007, 
         region == "Central_Coast") %>% 
  #filter(year >= 2007,
  #       region == "South_Coast") %>% 
  mutate(total_urchins = rowSums(across(c(swath_CENCOR,
                                          swath_MESFRAAD,
                                          swath_STRPURAD)))) %>% 
  drop_na(fish_SPUL, swath_PANINT, swath_CENCOR, swath_MESFRAAD, swath_STRPURAD) %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("Reference", "Partial", "Full"))) %>% 
  filter(total_urchins >= 0) %>% 
  mutate(
    before_hw = 1 * (year < 2014),
    during_mhw = 1 * between(year, 2014, 2016),
    after_mhw = 1 * (year >= 2017)
  ) %>% 
  mutate(
    reference = 1 * (mpa_status == "Reference"),
    partial = 1 * (mpa_status == "Partial"), # A dummy for Partial
    full = 1 * (mpa_status == "Full"), # A dummy for Full
    year_imp = case_when(is.na(Estab_Yr_1) ~ 0, .default = year - Estab_Yr_1),
    post = case_when(year_imp > 0 ~ 1, .default = 0)
  ) 

#### Two way effects models ####################################################
# TWFE model
# This regression specification retrieves the dynamic treatment effects.
# It includes fixed effects by site and year
data$log_urchins <- log(ifelse(data$total_urchins == 0, min(data$total_urchins[data$total_urchins > 0]), data$total_urchins))

TWFE_model <- feols(log_urchins ~  i(year_imp, partial, -1) + i(year_imp, full, -1) |
                      SITE + year, # Fixed effects by site and year,
                    panel.id = ~SITE+year, 
                    se = "DK", # Specify we want IID SEs
                    data = data)

#### Interpreting results 
summary(TWFE_model)

# If all we care about are the average treatment effects, then we can use a
# Difference-in-Differences estimation. Basically, take the change in time for
# each treatment (that's the first difference) and then take the difference in across treatments
# that's the second difference
# y = beta_1*X_i + beta_2*post_it:partial_i + beta_3*post_it:full_i + gamma_t + alpha_i epsilon_it
# alpha_i and gamma_t capture site and year-fixed effects
# The average treatment effects are beta_2 and beta_3
DiD <- feols(log_urchins ~ post:partial + post:full |
               SITE + year, # Fixed effects by site and year
             panel.id = ~SITE+year, 
             se = "DK", # Specify IID SEs
             data = data)

modelsummary(DiD,
             stars = T)

##### Build the event effects plot #############################################
# Get coefficients into tabular form
# For the TWFE
TWFE_coefs <- tidy(TWFE_model,
                   conf.int = T) %>%
  filter(str_detect(term, "::")) %>%
  mutate(treatment = str_extract(term, "full|partial"),
         year_imp = as.numeric(str_remove_all(str_extract(term, "::.+:"), ":")))


# For the DiD
DiD_coefs <- tidy(DiD,
                  conf.int = T)

# Plot
ggplot(data = TWFE_coefs,
       aes(x = year_imp, y = estimate, color = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "steelblue") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
  geom_point(x = -1, y = 0, inherit.aes = F, color = "black", size = 4) +
  geom_segment(x = -0.5, xend = 17,
               y = DiD_coefs$estimate[DiD_coefs$term == "post:partial"],
               yend = DiD_coefs$estimate[DiD_coefs$term == "post:partial"],
               color = "steelblue") +
  geom_segment(x = -0.5, xend = 17,
               y = DiD_coefs$estimate[DiD_coefs$term == "post:full"],
               yend = DiD_coefs$estimate[DiD_coefs$term == "post:full"],
               color = "red") +
  geom_pointrange(aes(ymin = conf.low,
                      ymax = conf.high),
                  fatten = 2) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(x = "Time to MPA implementation",
       y = "Urchin Abundances Estimate and 95% Confidence Intervals",
       title = "Event-study plot",
       subtitle = "Points show coefficient estimates (interpreted relative to the event = -1 in the reference treatment)",
       caption = "After accounting for covariate values x, site- and year-level fixed-effects, the difference in y between reference
       and other two treatments is not different from 0 in the pre-period. Then, 'partial' increases by 2 (dashed blue line shows first period)
       and 'full' increases by 5 / yr (dashed red line shows first period). The average treatment effect from the second regression
       (what one might test for with an ANOVA-like analysis) is shown as solid lines. These lines are the mean difference in the post-period relative to the reference level
       after accounting for temoral trends observed across all treatments. It's know as the difference-in-difference estimate.")

TWFE_urchins_plot <- ggplot(data = TWFE_coefs,
         aes(x = year_imp, y = estimate, color = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_point(x = -1, y = 0, inherit.aes = F, color = "black", size = 4) +
  geom_segment(x = -0.5, xend = 17,
               y = DiD_coefs$estimate[DiD_coefs$term == "post:partial"],
               yend = DiD_coefs$estimate[DiD_coefs$term == "post:partial"],
               color = "steelblue") +
  geom_segment(x = -0.5, xend = 17,
               y = DiD_coefs$estimate[DiD_coefs$term == "post:full"],
               yend = DiD_coefs$estimate[DiD_coefs$term == "post:full"],
               color = "red") +
  geom_pointrange(aes(ymin = conf.low,
                      ymax = conf.high),
                  fatten = 2) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(x = "Time since MPA implementation",
       y = "Urchin Abundances Estimate and 95% Confidence Intervals",
       title = "Event-study plot",
       )

TWFE_urchins_plot

# Export figure
# png("Figures/Two_way_effects_model_Urchins.png", width = 9, height = 6, 
#     units = "in", res = 600)
# TWFE_urchins_plot
# dev.off() 


##### Add in MHWs ##############################################################
# Estimation
# Fit a model with appropriate interaction and fixed effects terms. The coefficients of
# interest, in terms of MHW quesitons, are those on the interaction terms between mhw periods and MPA type
MPA_hw_model <- feols(log_urchins ~ partial:post + full:post + 
                        partial:during_mhw + full:during_mhw + 
                        partial:after_mhw + full:after_mhw | 
                        SITE + year,
                      panel.id = ~SITE+year,
                      se = "DK",
                      data = data)

MPA_hw_model

modelsummary(MPA_hw_model,
             stars = T)

TWFE_hw_coefs <- tidy(MPA_hw_model,
                   conf.int = T) %>%
  filter(str_detect(term, ":")) %>%
  mutate(treatment = str_extract(term, "full|partial"))


TWFE_hw_plot <- ggplot(data = TWFE_hw_coefs,
                            aes(x = year_imp, y = estimate, color = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_point(x = -1, y = 0, inherit.aes = F, color = "black", size = 4) +
  geom_segment(x = -0.5, xend = 17,
               y = DiD_coefs$estimate[DiD_coefs$term == "post:partial"],
               yend = DiD_coefs$estimate[DiD_coefs$term == "post:partial"],
               color = "steelblue") +
  geom_segment(x = -0.5, xend = 17,
               y = DiD_coefs$estimate[DiD_coefs$term == "post:full"],
               yend = DiD_coefs$estimate[DiD_coefs$term == "post:full"],
               color = "red") +
  geom_pointrange(aes(ymin = conf.low,
                      ymax = conf.high),
                  fatten = 2) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(x = "Time since MPA implementation",
       y = "Urchin Abundances Estimate and 95% Confidence Intervals",
       title = "Event-study plot",
  )

TWFE_hw_plot


##### Sheephead Biomass Models ###############################################
biomass_data <- data %>% 
  left_join(PISCO_biomass_data, by = c("SITE" = "site", "year", "mpa_status", "Estab_Yr_1")) 

# COME BACK HERE: Removed all NA's as you can't take a log of zero's
# Zero's here mean that there are no large sheephead (>= 35cm seen)
biomass_data <- biomass_data %>% 
  drop_na(avg_biomass_per_site) %>% 
  mutate(log_biomass = log(avg_biomass_per_site))
# 20% of the data are NA's :( 


hist(biomass_data$log_biomass)

TWFE_sheep_m <- feols(log_biomass ~  i(year_imp, partial, -1) + i(year_imp, full, -1) |
                      SITE + year, # Fixed effects by site and year,
                    panel.id = ~SITE+year, 
                    se = "DK", # Specify we want IID SEs
                    data = biomass_data)

#### Interpreting results 
summary(TWFE_sheep_m)

DiD_sheep <- feols(log_biomass ~ post:partial + post:full |
               SITE + year, # Fixed effects by site and year
               panel.id = ~SITE+year, 
               se = "DK", # Specify IID SEs
             data = biomass_data)


modelsummary(DiD_sheep,
             stars = T)

TWFE_coefs <- tidy(TWFE_sheep_m,
                   conf.int = T) %>%
  filter(str_detect(term, "::")) %>%
  mutate(treatment = str_extract(term, "full|partial"),
         year_imp = as.numeric(str_remove_all(str_extract(term, "::.+:"), ":")))


# For the DiD
DiD_coefs <- tidy(DiD_sheep,
                  conf.int = T)


ggplot(data = TWFE_coefs,
       aes(x = year_imp, y = estimate, color = treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_point(x = -1, y = 0, inherit.aes = F, color = "black", size = 4) +
  geom_segment(x = -0.5, xend = 17,
               y = DiD_coefs$estimate[DiD_coefs$term == "post:partial"],
               yend = DiD_coefs$estimate[DiD_coefs$term == "post:partial"],
               color = "steelblue") +
  geom_segment(x = -0.5, xend = 17,
               y = DiD_coefs$estimate[DiD_coefs$term == "post:full"],
               yend = DiD_coefs$estimate[DiD_coefs$term == "post:full"],
               color = "red") +
  geom_pointrange(aes(ymin = conf.low,
                      ymax = conf.high),
                  fatten = 2) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(x = "Time to MPA implementation",
       y = "Estimate and 95% CI of Biomass of large sheepheads (>= 35cm)",
       title = "Sheephead Event-study plot",
       subtitle = "Points show coefficient estimates (interpreted relative to the event = -1 in the reference treatment)",
       caption = "After accounting for covariate values x, site- and year-level fixed-effects, the difference in y between reference
       and other two treatments is not different from 0 in the pre-period. Then, 'partial' increases by 2 (dashed blue line shows first period)
       and 'full' increases by 5 / yr (dashed red line shows first period). The average treatment effect from the second regression
       (what one might test for with an ANOVA-like analysis) is shown as solid lines. These lines are the mean difference in the post-period relative to the reference level
       after accounting for temoral trends observed across all treatments. It's know as the difference-in-difference estimate.")

MPA_hw_model <- feols(log_biomass ~ partial:post + full:post + 
                        partial:during_mhw + full:during_mhw + 
                        partial:after_mhw + full:after_mhw | 
                        SITE + year,
                      panel.id = ~SITE+year,
                      se = "DK",
                      data = biomass_data)

MPA_hw_model

modelsummary(MPA_hw_model,
             stars = T)

##### Lobster Abundance Models ###############################################
### Problem here, you can't use a normal model for count data, as the model 
### should not predict negative lobsters, instead it should be zero lobsters. 
# 75% of the average lobster counts per site are zero 

# # Everything below this should be ignored 
# 
# TWFE_lobster_m <- feols(swath_PANINT ~  i(year_imp, partial, -1) + i(year_imp, full, -1) |
#                         SITE + year, # Fixed effects by site and year,
#                       panel.id = ~SITE+year, 
#                       se = "NW", # Specify we want IID SEs
#                       data = data)
# 
# #### Interpreting results 
# summary(TWFE_lobster_m)
# 
# DiD_lobster <- feols(swath_PANINT ~ post:partial + post:full |
#                      SITE + year, # Fixed effects by site and year
#                    se = "iid", # Specify IID SEs
#                    data = data)
# 
# modelsummary(DiD_lobster,
#              stars = T)
# 
# TWFE_coefs <- tidy(TWFE_lobster_m,
#                    conf.int = T) %>%
#   filter(str_detect(term, "::")) %>%
#   mutate(treatment = str_extract(term, "full|partial"),
#          year_imp = as.numeric(str_remove_all(str_extract(term, "::.+:"), ":")))
# 
# 
# # For the DiD
# DiD_coefs <- tidy(DiD_lobster,
#                   conf.int = T)
# 
# 
# ggplot(data = TWFE_coefs,
#        aes(x = year_imp, y = estimate, color = treatment)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_vline(xintercept = -1, linetype = "dashed") +
#   geom_hline(yintercept = 2, linetype = "dashed", color = "steelblue") +
#   geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
#   geom_point(x = -1, y = 0, inherit.aes = F, color = "black", size = 4) +
#   geom_segment(x = -0.5, xend = 17,
#                y = DiD_coefs$estimate[DiD_coefs$term == "post:partial"],
#                yend = DiD_coefs$estimate[DiD_coefs$term == "post:partial"],
#                color = "steelblue") +
#   geom_segment(x = -0.5, xend = 17,
#                y = DiD_coefs$estimate[DiD_coefs$term == "post:full"],
#                yend = DiD_coefs$estimate[DiD_coefs$term == "post:full"],
#                color = "red") +
#   geom_pointrange(aes(ymin = conf.low,
#                       ymax = conf.high),
#                   fatten = 2) +
#   scale_color_brewer(palette = "Set1") +
#   theme_minimal() +
#   labs(x = "Time to MPA implementation",
#        y = "Estimate and 95% CI of Lobster Abundances",
#        title = "Lobster Event-study plot",
#        subtitle = "Points show coefficient estimates (interpreted relative to the event = -1 in the reference treatment)",
#        caption = "After accounting for covariate values x, site- and year-level fixed-effects, the difference in y between reference
#        and other two treatments is not different from 0 in the pre-period. Then, 'partial' increases by 2 (dashed blue line shows first period)
#        and 'full' increases by 5 / yr (dashed red line shows first period). The average treatment effect from the second regression
#        (what one might test for with an ANOVA-like analysis) is shown as solid lines. These lines are the mean difference in the post-period relative to the reference level
#        after accounting for temoral trends observed across all treatments. It's know as the difference-in-difference estimate.")
# 
# 
# # Can we say that we shouldn't model sheephead biomass and lobster abundances because of 
# # the nature of the data... unless we can figure out a way to implement this type of modeling
# # for other distributions... specifically for zero-inflated data. 