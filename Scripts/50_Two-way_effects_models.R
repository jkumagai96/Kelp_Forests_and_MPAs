# Date: February 19th 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Two way effects models for abundances of urchins
# BIO 202: Ecological Statistics

##### Load packages and data ###################################################
library(fixest)
library(broom)
library(modelsummary)
library(tidyverse)

# Load Data
PISCO_data <- read.csv("Processed_data/PISCO_data_summarized.csv")
sites_for_joining <- read.csv("Processed_data/data_tables/PISCO_sites_with_MPAs.csv") %>% 
  rename(site = SITE)

# Only biomass for sheephead >= 35cm https://royalsocietypublishing.org/doi/full/10.1098/rspb.2016.1936

##### Format Data ##############################################################
data <- PISCO_data %>%  
  rename(total_urchins = urchin_d) %>% 
  left_join(sites_for_joining) %>% 
  filter(region == "South_Coast") %>% 
  drop_na(total_urchins) %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("Reference", "Partial", "Full"))) %>% 
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
                      site + year, # Fixed effects by site and year,
                    panel.id = ~site+year, 
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
               site + year, # Fixed effects by site and year
             panel.id = ~site+year, 
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
  geom_segment(x = -0.5, xend = 20,
               y = DiD_coefs$estimate[DiD_coefs$term == "post:partial"],
               yend = DiD_coefs$estimate[DiD_coefs$term == "post:partial"],
               color = "steelblue") +
  geom_segment(x = -0.5, xend = 20,
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
png("Figures/Two_way_effects_model_Urchins.png", width = 9, height = 6, 
    units = "in", res = 600)
TWFE_urchins_plot
dev.off() 


##### Add in MHWs ##############################################################
# Estimation
# Fit a model with appropriate interaction and fixed effects terms. The coefficients of
# interest, in terms of MHW quesitons, are those on the interaction terms between mhw periods and MPA type
MPA_hw_model <- feols(log_urchins ~ partial:post + full:post + 
                        partial:during_mhw + full:during_mhw + 
                        partial:after_mhw + full:after_mhw | 
                        site + year,
                      panel.id = ~site+year,
                      se = "DK",
                      data = data)

MPA_hw_model

modelsummary(MPA_hw_model,
             stars = T)

TWFE_hw_coefs <- tidy(MPA_hw_model,
                   conf.int = T) %>%
  filter(str_detect(term, ":")) %>%
  mutate(treatment = str_extract(term, "full|partial"))