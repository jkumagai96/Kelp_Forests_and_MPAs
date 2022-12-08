# Date: December 3rd 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Bootstrap - randomization approach 
# BIO 202: Ecological Statistics

###### Research Question #######################################################
# Do marine protected areas increase the resilience of kelp forests (increase in kelp area)
# following a marine heat wave?

# Is there more kelp area in marine protected areas? Partial vs. fully protected?

###### Set up ##################################################################
# Load Packages
library(tidyverse)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_per_year.csv")

##### Structure Data ###########################################################
kelp_data <- kelp_data_all %>% 
  # filter(year >= 2012) %>% # Latest MPAs were established in 20212
  mutate(mpa_status = factor(mpa_status, levels = c("None", "Partial", "Full"))) %>% 
  select(PixelID, year, area, mpa_status) %>% 
# pivot_wider(names_from = year, 
#              values_from = area) %>% 
  drop_na()

##### Calculate True Values ####################################################
# Change global option to not print out message about group summaries 
options(dplyr.summarise.inform = FALSE)

# Calculate average kelp area per category per year 
true_values <- kelp_data %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_area = mean(area)) %>% 
  pivot_wider(names_from = mpa_status,
              values_from = avg_area) %>% 
  mutate(P_N = Partial - None,
         F_N = Full - None,
         F_P = Full - Partial) %>% 
  select(-c(None, Partial, Full))

##### Bootstrapping ############################################################

# Set up values outside of the vector
years <- 1984:2021
mpa_status_list <- list()

for (i in 1:length(years)) {
  subset <- kelp_data %>% 
    filter(year == years[i]) 
  list[[i]] <- as.vector(subset$mpa_status)
}

new_data <- kelp_data %>% select(-mpa_status)

# So it is repeateable
set.seed(32)

# Where we will save the results
bootstrap_list <- list()

# within the for loop 
start <- Sys.time()
for (i in 1:10000) {
  values <- new_data %>% 
    group_by(year) %>% 
    mutate(mpa_status_random = sample(mpa_status, 
                                      size = length(mpa_status),
                                      replace = F)) %>% 
    ungroup() %>% 
    group_by(year, mpa_status_random) %>% 
    summarise(avg_area = mean(area)) %>% 
    pivot_wider(names_from = mpa_status,
                values_from = avg_area) %>% 
    mutate(P_N = Partial - None,
           F_N = Full - None,
           F_P = Full - Partial) %>% 
    select(-c(None, Partial, Full)) %>% 
    ungroup()
  
  bootstrap_list[[i]] <- values
}
end <- Sys.time()

start - end

#### Processing Results ########################################################
bootstrap_df <- do.call(rbind.data.frame, bootstrap_list)

# Each loop will be a row in the new dataframe
final_results <- data.frame()

for (i in 1:length(years)) {
  
  # Filter bootstrap values by year
  df <- bootstrap_df %>% 
    filter(year == years[i]) 
  
  tv <- true_values %>% 
    filter(year == years[i])
  
  percentile_FP <- ecdf(df$F_P)
  percentile_FN <- ecdf(df$F_N)
  percentile_PN <- ecdf(df$P_N)
  
  t <- data.frame(years[i], 
                  1 - percentile_FN(tv$F_N), 
                  1 - percentile_PN(tv$P_N),
                  1 - percentile_FP(tv$F_P))
  
  colnames(t) <- c("year", "F_N", "P_N", "F_P")
  
  final_results <- rbind(final_results, t)
  
}

0.05/30

# New P Value: 0.001666667

final_results <= 0.001666667
##### Export ###################################################################
final_results

write.csv(final_results, 
          "Processed_data/data_tables/bootstrap_2012_2021.csv", 
          row.names = F)

