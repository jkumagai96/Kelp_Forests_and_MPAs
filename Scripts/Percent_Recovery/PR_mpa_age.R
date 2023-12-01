# Date: July 21st 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Permutation Analysis simplified to during and after the heatwave (testing age of the MPAs)
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(sf)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")
points_in_mpas <- read.csv("Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv") 

##### Format Data ##############################################################
# Set cutoff year for delineation between old and new 
cutoff_yr <- 2007

maxes <- kelp_data_all %>% 
  filter(year < 2014) %>% 
  dplyr::select(PixelID, year, area) %>% 
  group_by(PixelID) %>% 
  summarize(historic_baseline = mean(area)) %>% 
  filter(historic_baseline > 0)

##### Calculate Percent Recovery ###############################################
pixels_to_remove <- points_in_mpas %>% 
  filter(Estab_Yr_1 <= cutoff_yr) %>%  # Marking all pixels within MPAs before or during cutoff to be removed
  select(PixelID) %>% 
  mutate(remove_me = 1)

# Calculate Percent Recovery for each row in the dataset from 2014 to 2021
kelp_data <- kelp_data_all %>% 
  filter(year >= 2014) %>% 
  left_join(maxes, by = "PixelID") %>% 
  select(PixelID, year, area, region, Mpa_ID, mpa_status, historic_baseline) %>% 
  mutate(percent_recovery = area/historic_baseline) %>%
  mutate(timeframe = ifelse(year > 2016, "after", "during")) %>% 
  group_by(PixelID, region, mpa_status, timeframe) %>% 
  summarize(mean_pr = mean(percent_recovery)) %>% 
  ungroup() %>% 
  na.omit(percent_recovery) %>% 
  left_join(pixels_to_remove, by = "PixelID") %>% 
  filter(is.na(remove_me)) 

points_in_mpas <- points_in_mpas %>% # Marking all pixels within MPAs to stay if established after cutoff year
  filter(Estab_Yr_1 > cutoff_yr)

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
write.csv(final_results, "Processed_data/data_tables/percent_recovery/pr_new_mpas.csv", row.names = F)

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
filename <- paste0("Figures/Percent_recovery/New_", cutoff_yr, "_mpas.png")
png(filename, width = 4, height = 3, units = "in", res = 600)
plot1
dev.off() 

kelp_data_new <- kelp_data %>% mutate(age = "new")
##### Old MPAs #################################################################
rm(list=setdiff(ls(), c("kelp_data_new", "cutoff_yr")))

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")
points_in_mpas <- read.csv("Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv") 

maxes <- kelp_data_all %>% 
  filter(year < 2014) %>% 
  dplyr::select(PixelID, year, area) %>% 
  group_by(PixelID) %>% 
  summarize(historic_baseline = mean(area)) %>% 
  filter(historic_baseline > 0)

pixels_to_remove <- points_in_mpas %>% 
  filter(Estab_Yr_1 > cutoff_yr) %>%  # Marking all pixels within MPAs before or during cutoff year to be removed
  select(PixelID) %>% 
  mutate(remove_me = 1)

# Calculate Percent Recovery for each row in the dataset from 2014 to 2021
kelp_data <- kelp_data_all %>% 
  filter(year >= 2014) %>% 
  left_join(maxes, by = "PixelID") %>% 
  select(PixelID, year, area, region, Mpa_ID, mpa_status, historic_baseline) %>% 
  mutate(percent_recovery = area/historic_baseline) %>%
  mutate(timeframe = ifelse(year > 2016, "after", "during")) %>% 
  group_by(PixelID, region, mpa_status, timeframe) %>% 
  summarize(mean_pr = mean(percent_recovery)) %>% 
  ungroup() %>% 
  na.omit(percent_recovery) %>% 
  left_join(pixels_to_remove, by = "PixelID") %>% 
  filter(is.na(remove_me)) 

points_in_mpas <- points_in_mpas %>% # Marking all pixels within MPAs to stay if established after cutoff
  filter(Estab_Yr_1 <= cutoff_yr)

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
write.csv(final_results, "Processed_data/data_tables/percent_recovery/pr_old_mpas.csv", row.names = F)

# Format Data into long format
results_long <- final_results %>% 
  pivot_longer(names_to = "Comparison", values_to = "pvalues", F_N:F_P) 

# Plot 
plot2 <- results_long %>% 
  ggplot(aes(x = timeframe, y = -log10(pvalues), group = Comparison)) +
  geom_point(aes(color = Comparison, shape = Comparison), size = 2) +
  geom_hline(yintercept = -log10(0.05/(6))) +
  scale_color_manual(values=c('#FF5C00', '#999999','#000EDD')) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(y = "-log10(P values)", x = "Time frame") 
plot2

# Export 
filename2 <- paste0("Figures/Percent_recovery/Old_", cutoff_yr, "_mpas.png")
png(filename2, width = 4, height = 3, units = "in", res = 600)
plot2
dev.off() 

kelp_data_old <- kelp_data %>% mutate(age = "old")


##### Percent Recovery Figure ##################################################
rm(kelp_data)

library(ggpubr)
library(rstatix)

# Significance Data
final_results_new<- read.csv("Processed_data/data_tables/percent_recovery/pr_new_mpas.csv") %>% 
  mutate(age = "new")
final_results_old <- read.csv("Processed_data/data_tables/percent_recovery/pr_old_mpas.csv") %>% 
  mutate(age = "old")

kelp_data <- rbind(kelp_data_old, kelp_data_new) %>% 
  mutate(timeframe = ifelse(timeframe == "after", "2017-2021", "2014-2016")) %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("Full", "Partial", "None")))

# Combine the final results and format 
final_results_all <- rbind(final_results_old, final_results_new)
results_long <- final_results_all %>% 
  pivot_longer(names_to = "Comparison", values_to = "pvalues", F_N:F_P) 

group.colors <- c(Full = "#440154", None = "#FFBA00", Partial ="#21918c")

base <- ggboxplot(kelp_data, 
          x = "mpa_status", 
          y = "mean_pr", 
          fill = "mpa_status",
          facet.by = c("age", "timeframe"),
          palette = group.colors, outlier.shape = NA) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     limits = c(0, 3)) +
  labs(y = "Relative Area", x = "") +
  theme(legend.position = "none")


# Using the final results above make the stat.test tibble that the figure reads
stat.test <- results_long %>% 
  mutate(group1 = c("Full", "Partial", "Full", "Full", "Partial", "Full",
                    "Full", "Partial", "Full", "Full", "Partial", "Full")) %>% 
  mutate(group2 = c("None", "None", "Partial", "None", "None", "Partial", 
                    "None", "None", "Partial", "None", "None", "Partial")) %>%
  mutate(p.adj = pvalues*6) %>% 
  select(timeframe, age, group1, group2, p.adj) %>% 
  mutate(y.position = c(2.85, 2.6, 2.6, 2.85, 2.6, 2.6,
                        2.85, 2.6, 2.6, 2.85, 2.6, 2.6)) %>% 
  mutate("p.adj.signif" = c("***", "**", "ns","ns", "ns", "ns",
                            "**", "ns", "**", "***", "ns", "***")) 
# create figure with significance bars 
boxplot_pr_w_sig <- base + 
  stat_pvalue_manual(stat.test, 
                     label = "p.adj.signif", 
                     tip.length = .00003, 
                     bracket.shorten = .05) +
  theme(strip.background =element_rect(fill="grey"))

boxplot_pr_w_sig 

png("Figures/Percent_recovery/pr_boxplot_w_sig_MPA_age.png", 
    width = 6, 
    height = 6, 
    units = "in", 
    res = 600)
boxplot_pr_w_sig 
dev.off() 
