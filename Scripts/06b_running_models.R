# Date: October 27th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Run models w/ log + 1 transformation
# BIO 202: Ecological Statistics

##### Set up ###################################################################
# load packages
library(tidyverse)
library(GGally)
library(ggeffects)
library(nlme) # Lme4 or LmerTest 
library(pscl)

##### Load Data and Format #####################################################
# Load data 
data <- read.csv("Processed_data/data_tables/kelp_data_w_mpas.csv")

# Define before heat wave and after heat wave time periods
new_data <- data %>% 
  filter(year < 2014 | year > 2016) %>% 
  mutate(hw = if_else(year < 2014, "before", "after"))

new_data$mpa_status <- factor(new_data$mpa_status, labels = c("None", "Partial", "Full"))
new_data$hw <- factor(new_data$hw, labels = c("before", "after"))

new_data$log_area <- log(new_data$area +1)

# Visualize
hist(new_data$area)
hist(new_data$log_area)

##### Running models ###########################################################
### Mixed effects linear models 
# Best MLME before 
Mlme1 <- lme(area ~ hw*mpa_status + temperature + depth,
              random = ~1 + mpa_status|PixelID,
              method = "REML",
              data = new_data)
plot(Mlme1)
plot(effects::allEffects(Mlme1)) 

Mlme2 <- lme(log_area ~ hw*mpa_status + temperature + depth,
             random = ~1 + mpa_status|PixelID,
             method = "REML",
             data = new_data)
plot(Mlme2)
plot(effects::allEffects(Mlme2)) 

### Hurdle Models 

mod.hurdle <- hurdle(round(log_area) ~ temperature + depth, 
                     data = new_data, 
                     dist = "negbin", 
                     zero.dist = "binomial") 

library(countreg)
rootogram(mod.hurdle, max = 15) 

mod.hurdle2 <- hurdle(round(log_area) ~ temperature + depth, 
                     data = new_data, 
                     dist = "poisson", 
                     zero.dist = "binomial") 
rootogram(mod.hurdle2, max = 15) 

AIC(mod.hurdle, mod.hurdle2) # mod.hurdle2 wins

mod.hurdle3 <- hurdle(round(log_area) ~ hw*mpa_status + temperature + depth + hsmax,
                      data = new_data, 
                      dist = "poisson", 
                      zero.dist = "binomial") 

rootogram(mod.hurdle3, max = 15) 
AIC(mod.hurdle2, mod.hurdle3) # mod.hurdle3 is way better 

mod.hurdle4 <- hurdle(round(log_area) ~ hw*mpa_status + temperature + hsmax |
                        hw*mpa_status + temperature + depth + hsmax,
                      data = new_data, 
                      dist = "poisson", 
                      zero.dist = "binomial") 

rootogram(mod.hurdle4, max = 15) 
AIC(mod.hurdle3, mod.hurdle4) # mod.hurdle4 is slightly better 


mod.hurdle5 <- hurdle(round(log_area) ~ hw + mpa_status + temperature + hsmax |
                        hw*mpa_status + temperature + depth + hsmax,
                      data = new_data, 
                      dist = "poisson", 
                      zero.dist = "binomial") 

rootogram(mod.hurdle5, max = 15) 
AIC(mod.hurdle4, mod.hurdle5) # mod.hurdle5 is slightly better 

residuals_data <- residuals(mod.hurdle5, type = "pearson")   
fitted_data <- fitted(mod.hurdle5) 
plot(residuals_data, fitted_data) # Doesn't look good 


## Try a Zero inflated model (negative binomial or poisson)

f1 <- formula(round(log_area) ~ hw*mpa_status + temperature + depth + hsmax)

ZINB <- zeroinfl(f1, dist = "negbin",
                 link = "logit", 
                 data = new_data)

summary(ZINB)
residuals_data <- residuals(ZINB, type = "pearson")   
fitted_data <- fitted(ZINB) 
plot(residuals_data, fitted_data)


# Poisson (slightly better fit)
ZIP <- zeroinfl(f1, dist = "poisson", 
               link = "logit", 
               data = new_data)
summary(ZIP)
residuals_data <- residuals(ZIP, type = "pearson")   
fitted_data <- fitted(ZIP) 
plot(residuals_data, fitted_data)

AIC(mod.hurdle5, ZIP, ZINB)

f2 <- formula(round(log_area) ~ hw*mpa_status + 
                                temperature + 
                                depth + 
                                hsmax + 
                                mpa_status)
ZIP2 <- zeroinfl(f2, dist = "poisson", 
                link = "logit", 
                data = new_data)
summary(ZIP2)
residuals_data <- residuals(ZIP2, type = "pearson")   
fitted_data <- fitted(ZIP2) 
plot(residuals_data, fitted_data)

AIC(mod.hurdle5, ZIP, ZINB, ZIP2) 
# Hurdle model is still best... but not great :( 


### Try a tweedie model
library(statmod)

# Tweedie1's have log transformed area
tweedie1 <- glm(log_area ~ hw*mpa_status + temperature + depth + hsmax,
                family = tweedie(var.power = 2,
                                 link.power = 0),
                data = new_data)
plot(tweedie1)


tweedie1b <- glm(log_area ~ hw*mpa_status + temperature + depth + hsmax,
                family = tweedie(var.power = 1,
                                 link.power = 0),
                data = new_data)
plot(tweedie1b)

# Tweedies 2 - 4 have just the regular area as it looks slightly better
tweedie2 <- glm(area ~ hw*mpa_status + temperature + depth,
                family = tweedie(var.power = 1.2,
                                 link.power = 0),
                data = new_data)
plot(tweedie2)

tweedie3 <- glm(area ~ hw*mpa_status + temperature + depth + hsmax,
                family = tweedie(var.power = 1.2,
                                 link.power = 0),
                data = new_data)
plot(tweedie3)

tweedie4 <- glm(area ~ hw*mpa_status + temperature + depth + hsmax,
                family = tweedie(var.power = 1.5,
                                 link.power = 0),
                data = new_data)
plot(tweedie4)


##### By-hand hurdle model #####################################################
# Presence absence model (logistic glm)

# create presence column 
new_data$presence <- ifelse(new_data$area == 0, 0, 1)

pa_model <- glm(presence ~ hw*mpa_status + temperature + depth + hsmax,
                family = binomial,
                data = new_data)


summary(pa_model)
plot(pa_model)

car::residualPlots(pa_model, fit = FALSE, tests = FALSE)
sqrt(car::vif(pa_model)) # Tests for colinearity of course with an interaction it will be larger 

pa_model2 <- glm(presence ~ hw + mpa_status + temperature + depth + hsmax,
                 family = binomial,
                 data = new_data)
summary(pa_model2)
plot(pa_model2)
sqrt(car::vif(pa_model2))

pa_model3 <- glm(presence ~ hw  + temperature + depth + hsmax,
                 family = binomial,
                 data = new_data)

summary(pa_model3)
plot(pa_model3)
sqrt(car::vif(pa_model3))

# The first pa model was best 
AIC(pa_model, pa_model2, pa_model3)

pa_model4 <- glm(presence ~ hw + mpa_status + temperature + depth + hsmax,
                 data = new_data,
                 family = binomial(link = "cloglog"))

summary(pa_model4)
plot(pa_model4)
sqrt(car::vif(pa_model4))

drop1(pa_model4, test = "Chi")

pa_model5 <- glm(presence ~ hw + mpa_status + temperature + depth,
                 data = new_data,
                 family = binomial(link = "cloglog"))

AIC(pa_model, pa_model2, pa_model3, pa_model4, pa_model5)

pa_model1_b <- glm(presence ~ hw*mpa_status + temperature + depth + hsmax,
                family = binomial(link = "cloglog"),
                data = new_data)

pa_model2_b <- glm(presence ~ hw + mpa_status + temperature + depth + hsmax,
                   family = binomial(link = "cloglog"),
                   data = new_data)

AIC(pa_model, pa_model2, pa_model1_b, pa_model2_b)

# The original PA_model is best (AIC is lower than all)


### Second half of the model
presence_data <- new_data %>% 
  filter(presence == 1)

Mlme3 <- lme(log_area ~ hw*mpa_status + temperature + depth + hsmax,
             random = ~1 + mpa_status|PixelID,
             method = "ML",
             data = presence_data)
summary(Mlme3)
plot(Mlme3)
plot(effects::allEffects(Mlme3)) 

Mlme4 <- lme(log_area ~ hw*mpa_status + temperature + hsmax,
             random = ~1 + mpa_status|PixelID,
             method = "ML",
             data = presence_data)
summary(Mlme4)
plot(Mlme4)
plot(effects::allEffects(Mlme4)) 

Mlme5 <- lme(log_area ~ hw*mpa_status + temperature,
             random = ~1 + mpa_status|PixelID,
             method = "ML",
             data = presence_data)

summary(Mlme5)
plot(Mlme5)
plot(effects::allEffects(Mlme5)) 

AIC(Mlme3, Mlme4, Mlme5) # Mlme4 is best 

# Try to vary the random effects
form_1 <- formula(log_area ~ hw*mpa_status + temperature + hsmax)

Mlme_a <- lme(form_1, random = ~1 + mpa_status|PixelID,
              method = "ML", 
              data = presence_data)

Mlme_b <- lme(form_1, random = ~1 + temperature|PixelID,
              method = "ML", 
              data = presence_data)
# didn't converge 

Mlme_c <- lme(form_1, random = ~1 + hsmax|PixelID,
              method = "ML", 
              data = presence_data)
# didn't converge 

Mlme_d <- lme(form_1, random = ~1 |PixelID,
              method = "ML", 
              data = presence_data)
AIC(Mlme_a, Mlme_d)

form_2 <- formula(log_area ~ hw*mpa_status + temperature + hsmax + depth)

Mlme2_a <- lme(form_2, random = ~1 + mpa_status|PixelID,
              method = "ML", 
              data = presence_data)

Mlme2_b <- lme(form_2, random = ~1 + temperature|PixelID,
              method = "ML", 
              data = presence_data)

Mlme2_c <- lme(form_2, random = ~1 + hsmax|PixelID,
              method = "ML", 
              data = presence_data)
# Didn't converge 

Mlme2_d <- lme(form_2, random = ~1 |PixelID,
              method = "ML", 
              data = presence_data)

AIC(Mlme2_a, Mlme2_b, Mlme2_d)
AIC(Mlme2_b,  Mlme_d)

Best_Mlme <- lme(log_area ~ hw*mpa_status + temperature + hsmax,
             random = ~1 + mpa_status|PixelID,
             method = "ML",
             data = presence_data)
# also confirmed this is the best one with a drop1(Mlme3)

Mlme2_b <- lme(form_2, random = ~1 + temperature|PixelID,
               method = "ML", 
               data = presence_data)

summary(Best_Mlme)
plot(Best_Mlme)
plot(effects::allEffects(Best_Mlme)) 

summary(Mlme2_b)
plot(Mlme2_b)
plot(effects::allEffects(Mlme2_b)) 

AIC(Best_Mlme, Mlme2_b)

Best_Mlme_2 <- lme(log_area ~ hw*mpa_status + temperature + hsmax + depth,
                   random = ~1 + temperature|PixelID,
                   method = "REML", 
                   data = presence_data)
summary(Best_Mlme_2)
plot(Best_Mlme_2)
plot(effects::allEffects(Best_Mlme_2)) 

# Calculate percent increase in full MPAs 
effects_table <- effects::allEffects(Best_Mlme_2)[[4]]
effects_after_hw_none <- exp(7.152034) - 1

effects_after_hw_full <- exp(8.031769) - 1

p_increase <- (effects_after_hw_full - effects_after_hw_none)/effects_after_hw_none*100
print(p_increase)


# GLMs
glm_model <- glm(log_area ~ hw*mpa_status + temperature + depth + hsmax,
                family = Gamma,
                data = presence_data)

plot(glm_model)

glm_model2 <- glm(log_area ~ hw*mpa_status + temperature + hsmax,
                 family = Gamma,
                 data = presence_data)

plot(glm_model2)


glm_model3 <- glm(log_area ~ hw*mpa_status + temperature,
                  family = Gamma,
                  data = presence_data)

summary(glm_model3)

AIC(glm_model, glm_model2, glm_model3) # GLM model 2 is best!


glm_model4 <- glm(log_area ~ hw*mpa_status + temperature + hsmax + year,
                  family = Gamma,
                  data = presence_data)

plot(glm_model4)
# adding year doesn't help. GLM model 2 is best. 




# Glm with mixed effects (not quite working)
library(lme4)
gm_me_1 <- glmer(log_area ~ hw*mpa_status + temperature + hsmax + 
                   (1 | PixelID),
              data = presence_data, 
              family = Gamma)

summary(gm_me_1)
