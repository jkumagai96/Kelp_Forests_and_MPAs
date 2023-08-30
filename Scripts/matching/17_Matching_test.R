# Date: June 22nd 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Matching for protected vs. not protected to predict percent recovery with the matching package
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(Matching)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")
percent_recovery <- read.csv("Processed_data/data_tables/relative_area.csv") %>% 
  dplyr::select(PixelID, year, mpa_status, percent_recovery)

# Format Data
full_and_none <- kelp_data_all %>% 
  filter(year > 2016) %>% 
  left_join(percent_recovery, by = c("PixelID", "year", "mpa_status")) %>% 
  filter(mpa_status != "Partial") %>%
  dplyr::select(-Mpa_ID) %>% 
  na.omit() %>% 
  mutate(mpa_status_binary = ifelse(mpa_status == "Full", 1, 0)) %>% 
  mutate(region_binary = ifelse(region == "Central_Coast", 1, 0))

Y <- full_and_none$percent_recovery 
Tr <- full_and_none$mpa_status_binary

hist(Y)
hist(log(Y)) # normalizes the data 

# Proposensity score model 
glm1 <- glm(Tr ~ year + hsmax + nitrate + temperature + MHW_intensity + 
              CS_intensity + region_binary + depth + gravity + distance_to_coast, 
            family = binomial, data = full_and_none)

# Matching
rr1 <- Match(Y = Y, Tr = Tr, X = glm1$fitted, estimand = "ATT", M = 3, ties = TRUE, replace = TRUE)


MatchBalance(Tr ~ year + hsmax + nitrate + temperature + MHW_intensity + CS_intensity + 
               region + depth + gravity + distance_to_coast, 
             match.out = rr1, nboots = 1000, data = full_and_none)

X <- cbind(full_and_none$year, full_and_none$hsmax, full_and_none$nitrate, full_and_none$temperature, full_and_none$MHW_intensity, full_and_none$CS_intensity, full_and_none$region_binary, full_and_none$depth, full_and_none$gravity, full_and_none$distance_to_coast) 

# Genetic Matching
gen1 <- GenMatch(Tr = Tr, X = X, weights = NULL, pop.size = 1000)

mgen1 <- Match(Y = Y, Tr = Tr, X = X, Weight.matrix = gen1) 
MatchBalance(Tr ~ year + hsmax + nitrate + temperature + MHW_intensity + CS_intensity + 
               region + depth + gravity + distance_to_coast, 
             data = full_and_none, match.out = mgen1, nboots = 1000) # check balance 

# Balance was not achieved. P value is still very low for gravity and some 
# variables seem to be worse after matching
# Those variables are the following:Distance to coast, gravity

# MatchBalance on distance to coast but do not include it in the genetic matching algorithm 
# Remove gravity?

# Genetic Matching
X2 <- cbind(full_and_none$year, full_and_none$hsmax, full_and_none$nitrate, full_and_none$temperature, full_and_none$MHW_intensity, full_and_none$CS_intensity, full_and_none$region_binary, full_and_none$depth) 

gen2 <- GenMatch(Tr = Tr, X = X2, weights = NULL, pop.size = 1000, wait.generations = 2)

mgen2 <- Match(Y = Y, Tr = Tr, X = X2, Weight.matrix = gen2) 
MatchBalance(Tr ~ year + hsmax + nitrate + temperature + MHW_intensity + CS_intensity + 
               region + depth, 
             data = full_and_none, match.out = mgen2, nboots = 1000) # check balance 

# Removing gravity and distance to coast, leads to semi-okay balance... needs to be improved 
summary(mgen2)

# Results are not significant for percent recovery! 
Y <- full_and_none$area
mgen2 <- Match(Y = Y, Tr = Tr, X = X2, Weight.matrix = gen2) 
# But they are significant for kelp area 


## Genetic Matching for all but distance to coast 
X3 <- cbind(full_and_none$year, full_and_none$hsmax, full_and_none$nitrate, full_and_none$temperature, full_and_none$MHW_intensity, full_and_none$CS_intensity, full_and_none$region_binary, full_and_none$depth, full_and_none$gravity) 

gen3 <- GenMatch(Tr = Tr, X = X3, weights = NULL, pop.size = 1000)

mgen3 <- Match(Y = Y, Tr = Tr, X = X3, Weight.matrix = gen3) 

MatchBalance(Tr ~ year + hsmax + nitrate + temperature + MHW_intensity + CS_intensity + 
               region + depth + gravity, 
             data = full_and_none, match.out = mgen3, nboots = 1000) # check balance 

#summary(mgen3)

# Genetic Matching full dataset, m =2 instead, did not work... 
gen3 <- GenMatch(Tr = Tr, X = X, M = 2, weights = NULL, pop.size = 1000)

mgen3 <- Match(Y = Y, Tr = Tr, X = X, Weight.matrix = gen3) 
MatchBalance(Tr ~ year + hsmax + nitrate + temperature + MHW_intensity + CS_intensity + 
               region + depth + gravity + distance_to_coast, 
             data = full_and_none, match.out = mgen3, nboots = 1000) # check balance 


hist(full_and_none$gravity)
group1 <- full_and_none %>% filter(mpa_status_binary == 1)
group2 <- full_and_none %>% filter(mpa_status_binary == 0)

hist(group1$gravity)
mean(group1$gravity)
hist(group2$gravity)
mean(group2$gravity) # On average, the unprotected kelp pixels are further away from human population centers (may be less impacted)
sd(group1$gravity)
sd(group2$gravity)

# Maybe adjust calipers? 
#v_calipers <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, .25, FALSE) #.25 is too low 
#v_calipers <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, .5, FALSE) #.25 is too low 
#v_calipers <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) #.25 is too low 

# TRY THIS NEXT

gen4 <- GenMatch(Tr = Tr, X = X, weights = NULL, pop.size = 1000, 
                 caliper = v_calipers, wait.generations = 2)

mgen4 <- Match(Y = Y, Tr = Tr, X = X, Weight.matrix = gen4) 
MatchBalance(Tr ~ year + hsmax + nitrate + temperature + MHW_intensity + CS_intensity + 
               region + depth + gravity + distance_to_coast, 
             data = full_and_none, match.out = mgen4, nboots = 1000) # check balance 


##### Genetic Matching starting with gravity ###################################
#X2 <- cbind(full_and_none$year, full_and_none$gravity) 
X2 <- cbind(full_and_none$year, full_and_none$gravity, full_and_none$temperature, 
            full_and_none$hsmax, full_and_none$region_binary)

gen2 <- GenMatch(Tr = Tr, X = X2, weights = NULL, pop.size = 1000, wait.generations = 3)

#Y <- full_and_none$percent_recovery 
Y <- full_and_none$area

mgen2 <- Match(Y = Y, Tr = Tr, X = X2, Weight.matrix = gen2) 

#MatchBalance(Tr ~ year + gravity, 
#             data = full_and_none, match.out = mgen2, nboots = 1000) # check balance 
MatchBalance(Tr ~ year + gravity + temperature + hsmax + region_binary, 
             data = full_and_none, match.out = mgen2, nboots = 1000) # check balance 

summary(mgen2) # significant only for kelp area and not percent recovery/relative area, 
# why is this different than what we find in the permutation analysis?

