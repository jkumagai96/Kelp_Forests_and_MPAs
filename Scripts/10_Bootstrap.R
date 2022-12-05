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
  filter(year >= 2012) %>% # Latest MPAs were established in 20212
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
mpa_vector <- kelp_data$mpa_status
new_data <- kelp_data %>% select(-mpa_status)

# So it is repeateable
set.seed(32)

# Where we will save the results
bootstrap_list <- list()

# within the for loop 
start <- Sys.time()
for (i in 1:10000) {
  values <- new_data %>% 
    mutate(mpa_status = sample(mpa_vector, 
                                size = length(mpa_vector), 
                                replace = F)) %>% 
    group_by(year, mpa_status) %>% 
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

years <- 2012:2021
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

final_results >= 0.001666667
##### Export ###################################################################
final_results

write.csv(final_results, 
          "Processed_data/data_tables/bootstrap_2012_2021.csv", 
          row.names = F)

#### Test 1985 to 1994 #########################################################
# Purpose: Test 1985 to 1994 for points going to be MPAs, but not yet, 
#         to see if before protection, there was still more area. Do they just 
#         protect areas where there is more kelp area/ healthier?

# Load extra packages and data
library(sf)
mpas_original <- read_sf("Data/Filtered_MPAs/MPAs_CA_Giant_Kelp.shp")
mpas <- mpas_original %>% 
  dplyr::select(Site_ID_12, Estab_Yr_1, AreaMar_12, Status_12)  %>% 
  rename("mpa_status" = "Status_12") 


# Identify points to include and not throw out 
# Throw out any mpas established before 1995, so include all at 1995 or later and NAs

sum(mpas$Estab_Yr_1 < 1995) # 12 MPAs which need to be removed

IDs_to_include <- mpas[mpas$Estab_Yr_1 >= 1995,]$Site_ID_12
IDs_to_include <- c(IDs_to_include, NA)

kelp_data <- kelp_data_all %>% 
  filter(Mpa_ID %in% IDs_to_include) %>% 
  filter(year > 1984 & year <= 1994) %>% # Just keep 1985 to 1994
  mutate(mpa_status = factor(mpa_status, levels = c("None", "Partial", "Full"))) %>% 
  select(PixelID, year, area, mpa_status) %>% 
  drop_na()

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

### Bootstrapping
# Set up values outside of the vector
mpa_vector <- kelp_data$mpa_status
new_data <- kelp_data %>% select(-mpa_status)

# So it is repeateable
set.seed(42)

# Where we will save the results
bootstrap_list <- list()

# within the for loop 
start <- Sys.time()
for (i in 1:10000) {
  values <- new_data %>% 
    mutate(mpa_status = sample(mpa_vector, 
                               size = length(mpa_vector), 
                               replace = F)) %>% 
    group_by(year, mpa_status) %>% 
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

bootstrap_df <- do.call(rbind.data.frame, bootstrap_list)

# Each loop will be a row in the new dataframe
final_results_1985_1994 <- data.frame()

years <- 1985:1994
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
  
  final_results_1985_1994 <- rbind(final_results_1985_1994, t)
  
}

final_results_1985_1994
final_results_1985_1994 >= 0.001666667


write.csv(final_results_1985_1994, 
          "Processed_data/data_tables/bootstrap_1985_1994.csv", 
          row.names = F)
