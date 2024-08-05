# Date: July 21st 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Permutation Analysis per region simplified to during and after the heatwave
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(sf)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")
points_in_mpas <- read.csv("Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv") 

##### Format Data ##############################################################
maxes <- kelp_data_all %>% 
  filter(year < 2014) %>% 
  dplyr::select(PixelID, region, year, area) %>% 
  group_by(PixelID) %>% 
  summarize(historic_baseline = mean(area)) %>% 
  filter(historic_baseline > 0)

maxes_south <- kelp_data_all %>% 
  filter(year < 2014,
         region == "South_Coast") %>% 
  dplyr::select(PixelID, region, year, area) %>% 
  group_by(PixelID) %>% 
  summarize(historic_baseline = mean(area)) %>% 
  filter(historic_baseline > 0)

maxes_central <- kelp_data_all %>% 
  filter(year < 2014,
         region == "Central_Coast") %>% 
  dplyr::select(PixelID, region, year, area) %>% 
  group_by(PixelID) %>% 
  summarize(historic_baseline = mean(area)) %>% 
  filter(historic_baseline > 0)


##### Calculate Percent Recovery ###############################################
# Calculate Percent Recovery for each row in the dataset from 2014 to 2021
kelp_data <- kelp_data_all %>% 
  filter(year >= 2014) %>% 
  left_join(maxes, by = "PixelID") %>% 
  select(PixelID, year, area, region, Mpa_ID, mpa_status, historic_baseline) %>% 
  mutate(percent_recovery = area/historic_baseline) %>%
  mutate(timeframe = ifelse(year > 2016, "after", "during")) %>% 
  group_by(PixelID, region, mpa_status, timeframe) %>% 
  summarize(mean_pr = mean(percent_recovery),
            sum_area = sum(area)) %>% 
  ungroup() %>% 
  na.omit(percent_recovery)

##### Structure Data ###########################################################
kelp_data_south_n <- kelp_data %>% 
  filter(region == "South_Coast")

kelp_data_central_n <- kelp_data %>% 
  filter(region == "Central_Coast")

points_in_mpas_south <- points_in_mpas %>% 
  filter(region == "South_Coast")

points_in_mpas_central <- points_in_mpas %>% 
  filter(region == "Central_Coast")

rm(maxes)

##### Start For Loop ###########################################################

percent_cutoffs <- c(0.1, 0.15, 0.2, 0.25, 0.3)

for (i in 1:length(percent_cutoffs)) {
  print(percent_cutoffs[i])

##### Remove lowest Pixels #####################################################
total_area_after2013 <- sum(kelp_data_central_n$sum_area) 

cutoff_percent <- percent_cutoffs[i]

cutoff_area <- quantile(maxes_central$historic_baseline, cutoff_percent) # Based on historical values 

cutoff_table <- maxes_central %>% 
  mutate(keep = ifelse(historic_baseline >= cutoff_area, 1, 0)) 

kelp_data_central <- kelp_data_central_n %>% 
  left_join(cutoff_table, by = c("PixelID")) %>% 
  filter(keep == 1) %>% 
  dplyr::select(-historic_baseline)

area_after2013 <- sum(kelp_data_central$sum_area) 

percent_area_removed <- 100 - 100*(area_after2013/total_area_after2013) # Percent of total removed from the dataset
percent_area_removed

##### Calculate True Values for central region #################################
# Calculate average kelp area per category per year 
true_values <- kelp_data_central %>% 
  group_by(timeframe, mpa_status) %>% 
  summarise(median_pr = median(mean_pr)) %>% 
  pivot_wider(names_from = mpa_status,
              values_from = median_pr) %>% 
  mutate(P_N = Partial - None,
         F_N = Full - None,
         F_P = Full - Partial) %>% 
  select(-c(None, Partial, Full))

# Change global option to not print out message about group summaries 
options(dplyr.summarise.inform = FALSE)

##### Bootstrapping ############################################################
## Set up:
set.seed(20) # So the results are repeatable

# Set up variables outside of the loops
bootstrap_list <- list()

# Remove original mpa_status from kelp_data
kelp_data_r <- kelp_data_central %>% select(-mpa_status)

all_pixels <- kelp_data_central$PixelID %>% unique()

# within the for loop 
start <- Sys.time()
for (j in 1:10000) {
  
  # Sample random pixels that will overlap w/ MPAs 
  r_pixels <- sample(all_pixels, size = nrow(points_in_mpas_central), replace = F) # number of pixels originally overlapping with MPAs 
  
  # Assign the random pixels to the points in the mpas instead of the original values
  points_in_mpas_central$PixelID <- r_pixels
  
  # Join randomized mpa data with kelp data 
  kelp_w_mpas_r <- left_join(kelp_data_r, points_in_mpas_central, by = "PixelID") 
  
  # Ensure that the mpa status accounts for the year of establishment
  kelp_w_mpas_r$mpa_status[is.na(kelp_w_mpas_r$mpa_status)] <- "None"
  
  new_data <- kelp_w_mpas_r
  
  # Calculate table 
  values <- new_data %>% 
    group_by(timeframe, mpa_status) %>% 
    summarise(median_pr = median(mean_pr)) %>% 
    pivot_wider(names_from = mpa_status,
                values_from = median_pr) %>% 
    mutate(P_N = Partial - None,
           F_N = Full - None,
           F_P = Full - Partial) %>% 
    select(-c(None, Partial, Full)) %>% 
    ungroup()
  
  bootstrap_list[[j]] <- values
}
end <- Sys.time()

start - end

#### Processing Results ########################################################
bootstrap_df <- do.call(rbind.data.frame, bootstrap_list)

# Each loop will be a row in the new dataframe
final_results_central <- data.frame()

# During 
df <- bootstrap_df %>% filter(timeframe == "during")
tv <- true_values %>% filter(timeframe == "during")

sig_F_N <- 1 - ecdf(df$F_N)(tv$F_N)
sig_P_N <- 1 - ecdf(df$P_N)(tv$P_N)
sig_F_P <- 1 - ecdf(df$F_P)(tv$F_P)

during_pvalues <- data.frame("2014-2016",
                             sig_F_N,
                             sig_P_N,
                             sig_F_P)

colnames(during_pvalues) <- c("timeframe","F_N", "P_N", "F_P")

# After 
df <- bootstrap_df %>% filter(timeframe == "after")
tv <- true_values %>% filter(timeframe == "after")

sig_F_N <- 1 - ecdf(df$F_N)(tv$F_N)
sig_P_N <- 1 - ecdf(df$P_N)(tv$P_N)
sig_F_P <- 1 - ecdf(df$F_P)(tv$F_P)

after_pvalues <- data.frame("2017-2021",
                            sig_F_N,
                            sig_P_N,
                            sig_F_P)

colnames(after_pvalues) <- c("timeframe","F_N", "P_N", "F_P")

final_results_central <- rbind(during_pvalues, after_pvalues) %>% 
  mutate(F_N = p.adjust(F_N, method = "bonferroni", n = 6),
         P_N = p.adjust(P_N, method = "bonferroni", n = 6),
         F_P = p.adjust(F_P, method = "bonferroni", n = 6))
print(final_results_central)

filename <- paste0("Processed_data/data_tables/percent_recovery/pr_central_", cutoff_percent, ".csv")

write.csv(final_results_central, filename, row.names = F)

rm(kelp_data_central, total_area_after2013, area_after2013, cutoff_area, cutoff_table)
rm(true_values, kelp_data_r, all_pixels, bootstrap_list)

##### Calculate True Values for south region ###################################

### Remove lowest Pixels 
total_area_after2013 <- sum(kelp_data_south_n$sum_area) 

cutoff_area <- quantile(maxes_south$historic_baseline, cutoff_percent) # Based on historical values 

cutoff_table <- maxes_south %>% 
  mutate(keep = ifelse(historic_baseline >= cutoff_area, 1, 0)) 

kelp_data_south <- kelp_data_south_n %>% 
  left_join(cutoff_table, by = c("PixelID")) %>% 
  filter(keep == 1) %>% 
  dplyr::select(-historic_baseline)

area_after2013 <- sum(kelp_data_south$sum_area) 

percent_area_removed <- 100 - 100*(area_after2013/total_area_after2013) # Percent of total removed from the dataset
percent_area_removed

# Calculate average kelp area per category per year 
true_values <- kelp_data_south %>% 
  group_by(timeframe, mpa_status) %>% 
  summarise(median_pr = median(mean_pr)) %>% 
  pivot_wider(names_from = mpa_status,
              values_from = median_pr) %>% 
  mutate(P_N = Partial - None,
         F_N = Full - None,
         F_P = Full - Partial) %>% 
  select(-c(None, Partial, Full))

##### Bootstrapping ############################################################
## Set up:
set.seed(20) # So the results are repeatable

# Set up variables outside of the loops
bootstrap_list <- list()

# Remove original mpa_status from kelp_data
kelp_data_r <- kelp_data_south %>% select(-mpa_status)

all_pixels <- kelp_data_south$PixelID %>% unique()

# within the for loop 
start <- Sys.time()
for (j in 1:10000) {
  
  # Sample random pixels that will overlap w/ MPAs 
  r_pixels <- sample(all_pixels, size = nrow(points_in_mpas_south), replace = F) # number of pixels originally overlapping with MPAs 
  
  # Assign the random pixels to the points in the mpas instead of the original values
  points_in_mpas_south$PixelID <- r_pixels
  
  # Join randomized mpa data with kelp data 
  kelp_w_mpas_r <- left_join(kelp_data_r, points_in_mpas_south, by = "PixelID") 
  
  # Ensure that the mpa status accounts for the year of establishment
  kelp_w_mpas_r$mpa_status[is.na(kelp_w_mpas_r$mpa_status)] <- "None"
  
  new_data <- kelp_w_mpas_r
  
  # Calculate table 
  values <- new_data %>% 
    group_by(timeframe, mpa_status) %>% 
    summarise(median_pr = median(mean_pr)) %>% 
    pivot_wider(names_from = mpa_status,
                values_from = median_pr) %>% 
    mutate(P_N = Partial - None,
           F_N = Full - None,
           F_P = Full - Partial) %>% 
    select(-c(None, Partial, Full)) %>% 
    ungroup()
  
  bootstrap_list[[j]] <- values
}
end <- Sys.time()

start - end

#### Processing Results ########################################################
bootstrap_df <- do.call(rbind.data.frame, bootstrap_list)

# Each loop will be a row in the new dataframe
final_results_south <- data.frame()

# During 
df <- bootstrap_df %>% filter(timeframe == "during")
tv <- true_values %>% filter(timeframe == "during")

sig_F_N <- 1 - ecdf(df$F_N)(tv$F_N)
sig_P_N <- 1 - ecdf(df$P_N)(tv$P_N)
sig_F_P <- 1 - ecdf(df$F_P)(tv$F_P)

during_pvalues <- data.frame("2014-2016",
                             sig_F_N,
                             sig_P_N,
                             sig_F_P)

colnames(during_pvalues) <- c("timeframe","F_N", "P_N", "F_P")

# After 
df <- bootstrap_df %>% filter(timeframe == "after")
tv <- true_values %>% filter(timeframe == "after")

sig_F_N <- 1 - ecdf(df$F_N)(tv$F_N)
sig_P_N <- 1 - ecdf(df$P_N)(tv$P_N)
sig_F_P <- 1 - ecdf(df$F_P)(tv$F_P)

after_pvalues <- data.frame("2017-2021",
                            sig_F_N,
                            sig_P_N,
                            sig_F_P)

colnames(after_pvalues) <- c("timeframe","F_N", "P_N", "F_P")

final_results_south <- rbind(during_pvalues, after_pvalues) %>% 
  mutate(F_N = p.adjust(F_N, method = "bonferroni", n = 6),
         P_N = p.adjust(P_N, method = "bonferroni", n = 6),
         F_P = p.adjust(F_P, method = "bonferroni", n = 6))

print(final_results_south)

filename <- paste0("Processed_data/data_tables/percent_recovery/pr_south_", cutoff_percent, ".csv")

write.csv(final_results_south, filename, row.names = F)

}