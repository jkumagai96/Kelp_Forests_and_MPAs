# Date: July 25th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Figure 2 - visualizing percent recovery with significance levels determined by permutation analysis 
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(sf)
library(ggpubr)
library(rstatix)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")

# Significance Data
final_results_south <- read.csv("Processed_data/data_tables/percent_recovery/pr_south.csv") %>% 
  mutate(region = "Southern")
final_results_central <- read.csv("Processed_data/data_tables/percent_recovery/pr_central.csv") %>% 
  mutate(region = "Central")
final_results <- read.csv("Processed_data/data_tables/percent_recovery/pr_all.csv") %>% 
  mutate(region = "All")

##### Format Data ##############################################################
# Combind the final results and format 
final_results_all <- rbind(final_results, final_results_central, final_results_south)

results_long <- final_results_all %>% 
  pivot_longer(names_to = "Comparison", values_to = "pvalues", F_N:F_P) 

##### Calculate Percent Recovery ###############################################
maxes <- kelp_data_all %>% 
  filter(year < 2014) %>% 
  dplyr::select(PixelID, year, area) %>% 
  group_by(PixelID) %>% 
  summarize(historic_baseline = mean(area)) %>% 
  filter(historic_baseline > 0)

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
  na.omit(percent_recovery)# Making Figure 2 for real this time 

##### Format Data for Figure ###################################################
kelp_data <- kelp_data %>% 
  mutate(timeframe = ifelse(timeframe == "after", "2017-2021", "2014-2016")) %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("Full", "Partial", "None"))) %>% 
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern"))

kelp_data2 <- kelp_data %>% 
  mutate(region = "All")

kelp_data_combo <- rbind(kelp_data, kelp_data2)

##### Create FIgure ############################################################
# Pallete 
group.colors <- c(Full = "#440154", None = "#FFBA00", Partial ="#21918c")

# Create base plot
base <- ggboxplot(kelp_data_combo, 
                  x = "mpa_status", 
                  y = "mean_pr", 
                  fill = "mpa_status",
                  facet.by = c("region", "timeframe"),
                  palette = group.colors, outlier.shape = NA) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     limits = c(0, 3)) +
  labs(y = "Relative Area", x = "") +
  theme(legend.position = "none")

base

# Using the final results above make the stat.test tibble that the figure reads
stat.test <- results_long %>% 
  mutate(group1 = c("Full", "Partial", "Full", "Full", "Partial", "Full",
                    "Full", "Partial", "Full", "Full", "Partial", "Full", 
                    "Full", "Partial", "Full", "Full", "Partial", "Full")) %>% 
  mutate(group2 = c("None", "None", "Partial", "None", "None", "Partial", 
                    "None", "None", "Partial", "None", "None", "Partial",
                    "None", "None", "Partial", "None", "None", "Partial")) %>%
  mutate(p.adj = pvalues*6) %>% 
  select(timeframe, region, group1, group2, p.adj) %>% 
  mutate(y.position = c(2.85, 2.6, 2.6, 2.85, 2.6, 2.6,
                        2.85, 2.6, 2.6, 2.85, 2.6, 2.6, 
                        2.85, 2.6, 2.6, 2.85, 2.6, 2.6)) %>% 
  mutate("p.adj.signif" = c("***", "ns", "ns", "**", "ns", "ns",
                            "ns", "*", "ns", "ns", "ns", "ns", 
                            "***", "ns", "**","*", "ns", "*")) 
# create figure with significance bars 
boxplot_pr_w_sig <- base + 
  stat_pvalue_manual(stat.test, 
                     label = "p.adj.signif", 
                     tip.length = .00003, 
                     bracket.shorten = .05) +
  theme(strip.background =element_rect(fill="grey"))

boxplot_pr_w_sig 

##### Export ###################################################################
png("Figures/Percent_recovery/pr_boxplot_w_sig.png", 
    width = 6, 
    height = 8, 
    units = "in", 
    res = 600)
boxplot_pr_w_sig 
dev.off() 

