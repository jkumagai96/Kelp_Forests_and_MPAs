# Date: August 1st 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Testing protection effects on abundances of lobsters, sheephead, and biomass 
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load packages
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(DHARMa)

# Load Data
data <- read.csv("Processed_data/MLPA_data_summarized_wo_siteblocks.csv") %>% 
  filter(region == "South_Coast") %>% 
  mutate(mpa_status = ifelse(mpa_status == "Reference", "Unprotected", mpa_status)) %>% 
  filter(!is.na(urchins_d)) %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("Unprotected", "Partial", "Full")))

# Declare group colors
group.colors <- c(Full = "#440154", Unprotected = "#FFBA00", Partial ="#21918c")

##### Format Data ##############################################################
after_df <- data %>% 
  filter(year >= 2012) %>% 
  mutate(year_fct = as.factor(year))

##### Modeling after ###########################################################
# sheephead
model_sheephead <- glmmTMB(
  SPUL_d ~ mpa_status + 
    (1 | site) + # random intercept 
    ar1(0 + year_fct | site), #  autoregressive function 
  data = after_df, family = tweedie(link = "log")
)

model_sheephead
simulationOutput_sheephead <- simulateResiduals(model_sheephead, plot = F)
plot(simulationOutput_sheephead)

model_sheephead2 <- glmmTMB(
  SPUL_d ~ mpa_status + 
    (1 + year | site), # random intercept and slopes
  data = after_df, family = tweedie(link = "log")
)

model_sheephead2 # convergence error


model_sheephead3 <- glmmTMB(
  SPUL_d ~ mpa_status + 
    (1 | site), # random intercept 
  data = after_df, family = tweedie(link = "log")
)

model_sheephead3
simulationOutput_sheephead3 <- simulateResiduals(model_sheephead3, plot = F)
plot(simulationOutput_sheephead3)

car::Anova(model_sheephead)
plot(emmeans(model_sheephead, specs = pairwise ~ mpa_status, type = "response"))
summary(model_sheephead)
emmeans(model_sheephead, specs = pairwise ~ mpa_status, type = "response")

# sheephead biomass
model_biomass <- glmmTMB(
  biomass_d ~ mpa_status + 
    (1 | site) + # random intercept 
    ar1(0 + year_fct | site), # autoregressive function 
  data = after_df, family = tweedie(link = "log")
)

model_biomass
simulationOutput_biomass <- simulateResiduals(model_biomass, plot = F)
plot(simulationOutput_biomass)

car::Anova(model_biomass)
plot(emmeans(model_biomass, specs = pairwise ~ mpa_status, type = "response"))
emmeans(model_biomass, specs = pairwise ~ mpa_status, type = "response")

# lobsters
model_lobsters <- glmmTMB(
  PANINT_d ~ mpa_status + 
    (1 | site) + # random intercept 
    ar1(0 + year_fct | site), # autoregressive function 
  data = after_df, family = tweedie(link = "log")
)

model_lobsters
simulationOutput_lobsters <- simulateResiduals(model_lobsters, plot = F)
plot(simulationOutput_lobsters)

car::Anova(model_lobsters)
plot(emmeans(model_lobsters, specs = pairwise ~ mpa_status, type = "response"))
emmeans(model_lobsters, specs = pairwise ~ mpa_status, type = "response")


