# Date: July 21st 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Permutation Analysis simplified to during and after the heatwave
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
  dplyr::select(PixelID, year, area) %>% 
  group_by(PixelID) %>% 
  summarize(historic_baseline = mean(area)) %>% 
  filter(historic_baseline > 0)

##### Calculate Percent Recovery ###############################################
# Calculate Percent Recovery for each row in the dataset from 2014 to 2021
kelp_data_individual <- kelp_data_all %>% 
  filter(year >= 2014) %>% 
  left_join(maxes, by = "PixelID") %>% 
  select(PixelID, year, area, region, Mpa_ID, mpa_status, historic_baseline) %>% 
  mutate(percent_recovery = area/historic_baseline) %>% 
  mutate(timeframe = ifelse(year > 2016, "after", "during"))
  
kelp_data <- kelp_data_individual %>%   
  group_by(PixelID, region, mpa_status, timeframe) %>% 
  summarize(mean_pr = mean(percent_recovery)) %>% 
  ungroup() %>% 
  na.omit(percent_recovery)

# Calculate total declines 
kelp_data_individual %>% 
  group_by(year, timeframe) %>% 
  summarize(sum_area = sum(area), 
            baseline = sum(historic_baseline, na.rm = T)) %>% 
  group_by(timeframe) %>% 
  summarize(mean_total_area = mean(sum_area), 
            baseline = mean(baseline)) %>% 
  mutate(p_loss = (1 - (mean_total_area / baseline))*100)

kelp_data_individual %>% 
  group_by(year, timeframe, region) %>% 
  summarize(sum_area = sum(area), 
            sum_baseline = sum(historic_baseline, na.rm = T)) %>% 
  group_by(timeframe, region) %>% 
  summarize(mean_total_area = mean(sum_area), 
            baseline = mean(sum_baseline)) %>% 
  mutate(p_loss = (1 - (mean_total_area / baseline))*100)

# SI Table 5, relative area means and medians  
kelp_data %>% 
  group_by(timeframe, mpa_status) %>% 
  summarise(median = median(mean_pr),
    mean_pr = mean(mean_pr)) 

kelp_data %>% 
  group_by(timeframe, mpa_status, region) %>% 
  summarise(median = median(mean_pr),
            mean_pr = mean(mean_pr)) %>% 
  arrange(region, timeframe)

##### Calculate True Values  ###################################################
# Change global option to not print out message about group summaries 
options(dplyr.summarise.inform = FALSE)

# Calculate average kelp area per category per year 
true_values <- kelp_data %>% 
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
kelp_data_r <- kelp_data %>% select(-mpa_status)

all_pixels <- kelp_data$PixelID %>% unique()

# within the for loop 
start <- Sys.time()
for (j in 1:10000) {
  
  # Sample random pixels that will overlap w/ MPAs 
  r_pixels <- sample(all_pixels, size = nrow(points_in_mpas), replace = F) # number of pixels originally overlapping with MPAs 
  
  # Assign the random pixels to the points in the mpas instead of the original values
  points_in_mpas$PixelID <- r_pixels
  
  # Join randomized mpa data with kelp data 
  kelp_w_mpas_r <- left_join(kelp_data_r, points_in_mpas, by = "PixelID") 
  
  # Set MPA status to none where no join occured
  kelp_w_mpas_r$mpa_status[is.na(kelp_w_mpas_r$mpa_status)] <- "None"
  
  # Calculate table 
  values <- kelp_w_mpas_r %>% 
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
final_results <- data.frame()

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

final_results <- rbind(during_pvalues, after_pvalues)

final_results

write.csv(final_results, "Processed_data/data_tables/percent_recovery/pr_all.csv", row.names = F)

##### Figures ##################################################################
# Format Data into long format
results_long <- final_results %>% 
  pivot_longer(names_to = "Comparison", values_to = "pvalues", F_N:F_P) 

# Plot 
plot1 <- results_long %>% 
  ggplot(aes(x = timeframe, y = -log10(pvalues), group = Comparison)) +
  geom_point(aes(color = Comparison, shape = Comparison), size = 2) +
  geom_hline(yintercept = -log10(0.05/(6))) +
  scale_color_manual(values=c('#FF5C00', '#999999','#000EDD')) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(y = "-log10(P values)", x = "Time frame") 
plot1

# Export 
png("Figures/Percent_recovery/All.png", width = 4, height = 3, units = "in", res = 600)
plot1
dev.off() 
