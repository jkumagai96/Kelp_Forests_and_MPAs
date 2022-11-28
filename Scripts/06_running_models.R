# Date: October 16th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Run models
# BIO 202: Ecological Statistics

##### Set up ###################################################################
# load packages
library(tidyverse)
library(GGally)
library(ggeffects)
library(nlme)

##### Load Data and Format #####################################################
# Load data 
data <- read.csv("Processed_data/data_tables/kelp_data_w_mpas.csv")

# Define before heat wave and after heat wave time periods
new_data <- data %>% 
  filter(year < 2014 | year > 2016) %>% 
  mutate(hw = if_else(year < 2014, "before", "after"))

new_data$mpa_status <- factor(new_data$mpa_status, labels = c("None", "Partial", "Full"))
new_data$hw <- factor(new_data$hw, labels = c("before", "after"))

##### Running models ###########################################################


m1 <- lm(area ~ hw*mpa_status + temperature + hsmax + depth, data = new_data)
summary(m1)
plot(m1) # Obviously something is wrong


test <- ggeffects::ggeffect(m1)
plot(test)

m2 <- lm(area ~ hw*mpa_status + temperature + hsmax, data = new_data)
summary(m2)

test <- ggeffects::ggeffect(m2)
plot(test)

# Try a model with mixed effect
library(nlme)
Mlme1 <- lme(area ~ hw*mpa_status + temperature + hsmax, 
             random = ~1 | PixelID, 
             data = new_data)

summary(Mlme1)
plot(Mlme1) # really bad
plot(effects::allEffects(Mlme1))

Mlme2 <- lme(area ~ hw*mpa_status + temperature + hsmax, 
             random = ~1 + mpa_status | PixelID, 
             method = "ML",
             data = new_data)

summary(Mlme2) # Much better AIC 
plot(Mlme2) # really bad again 
plot(effects::allEffects(Mlme2))

Mlme3 <- lme(area ~ hw*mpa_status + temperature, 
             random = ~1 + mpa_status| PixelID, 
             data = new_data)
summary(Mlme3) # higher AIC
plot(Mlme3)
plot(effects::allEffects(Mlme3))

# what about year?
Mlme4 <- lme(area ~ hw*mpa_status + temperature + year,
             random = ~1 + mpa_status|PixelID,
             data = new_data)
summary(Mlme4)
plot(Mlme4)
plot(effects::allEffects(Mlme4))

Mlme5 <- lme(area ~ hw*mpa_status + temperature + year + hsmax,
             random = ~1 + mpa_status|PixelID,
             data = new_data)
summary(Mlme5)
plot(Mlme5)
plot(effects::allEffects(Mlme5))


# Let me try to go the other way with ML not REML
Mlme6 <- lme(area ~ hw*mpa_status + temperature + year + hsmax + depth,
             random = ~1 + mpa_status|PixelID,
             method = "ML",
             data = new_data)

summary(Mlme6)
plot(Mlme6)
plot(effects::allEffects(Mlme6))
# AIC: 345,927.5 

Mlme7 <- lme(area ~ hw*mpa_status + temperature + hsmax + depth,
             random = ~1 + mpa_status|PixelID,
             method = "ML",
             data = new_data)

summary(Mlme7)
plot(Mlme7)
plot(effects::allEffects(Mlme7))
# AIC: 345,926.2

Mlme8 <- lme(area ~ hw*mpa_status + temperature + depth,
             random = ~1 + mpa_status|PixelID,
             method = "ML",
             data = new_data)

summary(Mlme8)
plot(Mlme8)
plot(effects::allEffects(Mlme8))
# AIC: 345,925

Mlme9 <- lme(area ~ hw*mpa_status + temperature,
             random = ~1 + mpa_status|PixelID,
             method = "ML",
             data = new_data) # not better

# Best
Mlme10 <- lme(area ~ hw*mpa_status + temperature + depth,
              random = ~1 + mpa_status|PixelID,
              method = "REML",
              data = new_data)

summary(Mlme10)
plot(Mlme10)
plot(effects::allEffects(Mlme10)) # stillr really strange

##### Hurdle Model #############################################################
# Try a hurdle model
# https://data.library.virginia.edu/getting-started-with-hurdle-models/

# load package 
library(AER)

# Simplify data at first
new_data <- new_data %>% 
  filter(year == 2021)

hist(new_data$area)
sum(new_data$area == 0) / nrow(new_data) # 34.5% of the data are zeros

new_data$area <- round(new_data$area) %>% as.integer()

mod.hurdle <- hurdle(area ~ temperature + depth, 
                     data = new_data, 
                     dist = "negbin", 
                     zero.dist = "binomial") 

summary(mod.hurdle)

library(countreg)
rootogram(mod.hurdle, max = 100) # fit up to count 162103

short_data <- new_data %>% 
  dplyr::select(PixelID, depth, mpa_status, area, hsmax, temperature)

mod.hurdle2 <- hurdle(area ~ ., 
                     data = short_data, 
                     dist = "negbin", 
                     zero.dist = "binomial") 
summary(mod.hurdle2)
rootogram(mod.hurdle2, max = 100) # fit up to count 162103

AIC(mod.hurdle)
AIC(mod.hurdle2) # Slightly better


mod.hurdle3 <- hurdle(area ~ ., 
                      data = short_data, 
                      dist = "poisson", 
                      zero.dist = "binomial") 
summary(mod.hurdle3)
rootogram(mod.hurdle3, max = 100) # WAY WAY WORSE


mod.hurdle4 <- hurdle(area ~ . | depth + temperature,
                      data = short_data, 
                      dist = "negbin", 
                      zero.dist = "binomial")

summary(mod.hurdle4)
rootogram(mod.hurdle4, max = 100) 

### Practive
# Install the AER package (Applied Econometrics with R)
# install.packages("AER")

# load package 
library(AER)

# load the data
data("NMES1988")         

# select certain columns; Col 1 is number of visits
nmes <- NMES1988[, c(1, 6:8, 13, 15, 18)]  


plot(table(nmes$visits))
sum(nmes$visits < 1)

mod1 <- glm(visits ~ ., data = nmes, family = "poisson")

# predict expected mean count
mu <- predict(mod1, type = "response")

# sum the probabilities of a 0 count for each mean
exp <- sum(dpois(x = 0, lambda = mu))

# predicted number of 0's
round(exp) 

# observed number of 0's
sum(nmes$visits < 1)                      


  
library(pscl)
mod.hurdle <- hurdle(visits ~ ., data = nmes, dist = "poisson", zero.dist = "binomial") 

summary(mod.hurdle)


# Need to install from R-Forge instead of CRAN
# install.packages("countreg", repos="http://R-Forge.R-project.org")
library(countreg)
rootogram(mod.hurdle, max = 80) # fit up to count 80

mod.hurdle.nb <- hurdle(visits ~ ., data = nmes, dist = "negbin")
rootogram(mod.hurdle.nb, max = 80) # fit up to count 80

AIC(mod.hurdle)
AIC(mod.hurdle.nb)
# same as this:
# mod.hurdle <- hurdle(visits ~ ., data = nmes, dist = "poisson", zero.dist = "binomial") 

long <- c()

for (i in 1:nrow(short_data)) {
  print(short_data[i,])
}



##### October 26th #############################################################
# Load Packages
library(tidyverse)
library(pscl)
library(nlme)

# Load Data
data <- read.csv("Processed_data/data_tables/kelp_data_w_mpas.csv")

# Define before heat wave and after heat wave time periods
new_data <- data %>% 
  filter(year < 2014 | year > 2016) %>% 
  mutate(hw = if_else(year < 2014, "before", "after"))

new_data$mpa_status <- factor(new_data$mpa_status, labels = c("None", "Partial", "Full"))
new_data$hw <- factor(new_data$hw, labels = c("before", "after"))

# Explore distributions
par(mfrow = c(1, 2))
hist(data$area)
hist(log(data$area)) # All zero's go to -INF, so similar to a hurdle model if 
                     # we remove the zero's

hist(sqrt(data$area)) 

# Best previous model 
Mlme_a <- lme(area ~ hw*mpa_status + temperature + depth,
              random = ~1 + mpa_status|PixelID,
              method = "REML",
              data = new_data)
plot(Mlme_a)
plot(effects::allEffects(Mlme_a)) # stillr really strange

# Try to use the same model, but change area to sqrt of area
Mlme_b <- lme(sqrt(area) ~ hw*mpa_status + temperature + depth,
              random = ~1 + mpa_status|PixelID,
              method = "REML",
              data = new_data)

plot(Mlme_b) # Slightly better, but still not great... 
plot(effects::allEffects(Mlme_b))
summary(Mlme_b)

# Try a GLM with a gamma distribution?
glm_a <- glm(area ~ hw*mpa_status + temperature + depth,
              #random = ~1 + mpa_status|PixelID,
              family = Gamma(link = "identity"),
              data = new_data)

# Error in eval(family$initialize) : non-positive values not allowed for the 'Gamma' family

glm_a <- glm(round(area) ~ hw*mpa_status + temperature + depth,
             #random = ~1 + mpa_status|PixelID,
             family = poisson,
             data = new_data)

plot(glm_a) # Really bad

# Remove zero's first
test <- new_data[new_data$area != 0, ] 

glm_b <- glm(log(area) ~ hw*mpa_status + temperature + depth,
             #random = ~1 + mpa_status|PixelID,
             family = Gamma(link = "identity"),
             data = test)
summary(glm_b)
plot(glm_b)
plot(effects::allEffects(glm_b)) 

glm_c <- glm(log(area) ~ hw*mpa_status + temperature + depth,
             #random = ~1 + mpa_status|PixelID,
             family = Gamma(link = "log"),
             data = test)

# GLM B is better! 
summary(glm_c)
plot(glm_c)
plot(effects::allEffects(glm_c)) 

glm_d <- glm(log(area) ~ hw*mpa_status + temperature + depth + hsmax,
             #random = ~1 + mpa_status|PixelID,
             family = Gamma(link = "identity"),
             data = test)
summary(glm_d)
plot(glm_d)
plot(effects::allEffects(glm_d)) 
# glm_d has a better AIC value, but the plots don't look better...


glm_e <- glm(log(area) ~ hw*mpa_status + temperature + depth + hsmax + year,
             #random = ~1 + mpa_status|PixelID,
             family = Gamma(link = "identity"),
             data = test)
summary(glm_e)
plot(glm_e)
plot(effects::allEffects(glm_e)) 
# glm_e has a better AIC value, but the plots don't look better...


####### Try a Zero inflated negative binomial ##################################
f1 <- formula(round(area) ~ hw*mpa_status + temperature + depth + hsmax + mpa_status)

ZIHB <- zeroinfl(f1, dist = "negbin",
                   link = "logit", 
                   data = new_data)
summary(ZIHB)
residuals_data <- residuals(ZIHB, type = "pearson")   
fitted_data <- fitted(ZIHB) 

plot(residuals_data, fitted_data)


ZIP <- zeroinfl(f1, dist = "poisson",
                 link = "logit", 
                 data = new_data) # ZIHB is much better

f2 <- formula(round(area) ~ hw*mpa_status + temperature + depth)
ZIHB_2 <- zeroinfl(f2, dist = "negbin",
                 link = "logit", 
                 data = new_data)

residuals_data <- residuals(ZIHB_2, type = "pearson")   
fitted_data <- fitted(ZIHB_2) 

AIC(ZIHB, ZIHB_2)

f3 <- formula(round(area) ~ hw*mpa_status + temperature + hsmax)

ZIHB_3 <- zeroinfl(f3, dist = "negbin",
                 link = "logit", 
                 data = new_data)

# Formula 3 is best. 



######## Try a tweedie model
library(statmod)

y <- rgamma(20,shape=5)
x <- 1:20
# Fit a poisson generalized linear model with identity link
test <- glm(y~x,family=tweedie(var.power=1,link.power=1))

tweedie1 <- glm(area ~ hw*mpa_status + temperature + depth,
                family = tweedie(var.power = 1.2,
                                 link.power = 0),
                data = new_data)
plot(tweedie1)

