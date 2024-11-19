# Date: Nov. 11th 2024
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
  #filter(!is.na(urchins_d)) %>% # use all data that is available 
  mutate(mpa_status = factor(mpa_status, levels = c("Unprotected", "Partial", "Full")))

# Declare group colors
group.colors <- c(Full = "#440154", Unprotected = "#FFBA00", Partial ="#21918c")

##### Format Data ##############################################################
after_df <- data %>% 
  filter(year >= 2012) %>% 
  mutate(year_fct = as.factor(year),
         year_std = (year - 2017)/sd(year))

##### Modeling sheephead densities after 2011 ##################################
# Random intercept 
model_sheephead1 <- glmmTMB(
  SPUL_d ~ mpa_status + 
    (1 | site), # random intercept 
  data = after_df, family = tweedie(link = "log")
)
model_sheephead1
simulationOutput_sheephead1 <- simulateResiduals(model_sheephead1, plot = F)
plot(simulationOutput_sheephead1) # good residuals 

# Random intercept and slopes
model_sheephead2 <- glmmTMB(
  SPUL_d ~ mpa_status + 
    (1 + year_std | site), # random intercept and slopes
  data = after_df, family = tweedie(link = "log")
)

model_sheephead2
simulationOutput_sheephead2 <- simulateResiduals(model_sheephead2, plot = F)
plot(simulationOutput_sheephead2)

# Random intercept and slopes and autoregressive function 
model_sheephead3 <- glmmTMB(
  SPUL_d ~ mpa_status + 
    (1 + year_std| site) +
    ar1(0 + year_fct | site),  
  data = after_df, family = tweedie(link = "log")
) # Convergence problem, random slopes are redundant 

model_sheephead3
simulationOutput_sheephead3 <- simulateResiduals(model_sheephead3, plot = F)
plot(simulationOutput_sheephead3)


# Random intercept and autoregressive function 
model_sheephead4 <- glmmTMB(
  SPUL_d ~ mpa_status + 
    (1 | site) + # random intercept 
    ar1(0 + year_fct | site), #  autoregressive function 
  data = after_df, family = tweedie(link = "log")
)

model_sheephead4
simulationOutput_sheephead4 <- simulateResiduals(model_sheephead4, plot = F)
plot(simulationOutput_sheephead4)

###### Model Comparison ########################################################
car::Anova(model_sheephead1)
car::Anova(model_sheephead2)
car::Anova(model_sheephead4)
# all agree that MPA status matters

AIC(model_sheephead1, model_sheephead2, model_sheephead4)
# Model 4 has the lowest AIC and is chosen 

plot(emmeans(model_sheephead4, specs = pairwise ~ mpa_status, type = "response"))
summary(model_sheephead4)
emmeans(model_sheephead4, specs = pairwise ~ mpa_status, type = "response")

##### Sheephead Biomass Models #################################################
# Random intercept 
model_biomass1 <- glmmTMB(
  biomass_d ~ mpa_status + 
    (1 | site), 
  data = after_df, family = tweedie(link = "log")
)
model_biomass1
simulationOutput_biomass1 <- simulateResiduals(model_biomass1, plot = F)
plot(simulationOutput_biomass1) # good residuals

# Random intercept and slopes
model_biomass2 <- glmmTMB(
  biomass_d ~ mpa_status + 
    (1 + year_std| site), 
  data = after_df, family = tweedie(link = "log")
)
model_biomass2
simulationOutput_biomass2 <- simulateResiduals(model_biomass2, plot = F)
plot(simulationOutput_biomass2)

# Random intercept and slopes and autoregressive function 
model_biomass3 <- glmmTMB(
  biomass_d ~ mpa_status + 
    (1 + year_std | site) +
    ar1(0 + year_fct | site), 
  data = after_df, family = tweedie(link = "log")
)
model_biomass3
simulationOutput_biomass3 <- simulateResiduals(model_biomass3, plot = F)
plot(simulationOutput_biomass3)

# Random intercepts and autoregressive function 
model_biomass4 <- glmmTMB(
  biomass_d ~ mpa_status + 
    (1 | site) +
    ar1(0 + year_fct | site), 
  data = after_df, family = tweedie(link = "log")
)
model_biomass4
simulationOutput_biomass4 <- simulateResiduals(model_biomass4, plot = F)
plot(simulationOutput_biomass4) # better residuals

##### Model Comparison Biomass #################################################
AIC(model_biomass1, model_biomass2, model_biomass3, model_biomass4)
# Model 3 with the lowest AIC, but model 4 has better residuals and is consistent with others

car::Anova(model_biomass1)
car::Anova(model_biomass2)
car::Anova(model_biomass3)
car::Anova(model_biomass4) 
# all agree that MPA status is significant

# Lowest AIC is chosen (model 3) 
summary(model_biomass3)
plot(emmeans(model_biomass3, specs = pairwise ~ mpa_status, type = "response"))

##### Lobster density models ###################################################
# Random intercepts
model_lobsters1 <- glmmTMB(
  PANINT_d ~ mpa_status + 
    (1 | site),
  data = after_df, family = tweedie(link = "log")
)
model_lobsters1
simulationOutput_lobsters1 <- simulateResiduals(model_lobsters1, plot = F)
plot(simulationOutput_lobsters1) 


# Random intercepts and slopes
model_lobsters2 <- glmmTMB(
  PANINT_d ~ mpa_status + 
    (1 + year_std | site),
  data = after_df, family = tweedie(link = "log")
)
model_lobsters2
simulationOutput_lobsters2 <- simulateResiduals(model_lobsters2, plot = F)
plot(simulationOutput_lobsters2) 

# Random intercepts, slopes, and autoregressive function
model_lobsters3 <- glmmTMB(
  PANINT_d ~ mpa_status + 
    (1 + year_std | site) + # random intercept 
    ar1(0 + year_fct | site), # autoregressive function 
  data = after_df, family = tweedie(link = "log")
)
model_lobsters3
simulationOutput_lobsters3 <- simulateResiduals(model_lobsters3, plot = F)
plot(simulationOutput_lobsters3) 

# Random intercepts and autoregressive functions 
model_lobsters4 <- glmmTMB(
  PANINT_d ~ mpa_status + 
    (1 | site) + # random intercept 
    ar1(0 + year_fct | site), # autoregressive function 
  data = after_df, family = tweedie(link = "log")
)
model_lobsters4
simulationOutput_lobsters4 <- simulateResiduals(model_lobsters4, plot = F)
plot(simulationOutput_lobsters4) 

##### Model comparison lobsters ################################################
AIC(model_lobsters1, model_lobsters2, model_lobsters3, model_lobsters4)
# Model 3 has lowest residuals 

car::Anova(model_lobsters1)
car::Anova(model_lobsters2)
car::Anova(model_lobsters3)
car::Anova(model_lobsters4)
# all in agreement

summary(model_lobsters3) # Model 3 chosen 
plot(emmeans(model_lobsters3, specs = pairwise ~ mpa_status, type = "response"))


