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
  na.omit(percent_recovery) 

##### Format Data for Figure ###################################################
kelp_data <- kelp_data %>% 
  mutate(timeframe = ifelse(timeframe == "after", "2017-2021", "2014-2016")) %>% 
  mutate(mpa_status = ifelse(mpa_status == "None", "Unprotected", mpa_status)) %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("Full", "Partial", "Unprotected"))) %>% 
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern"))

kelp_data2 <- kelp_data %>% 
  mutate(region = "All")

kelp_data_combo <- rbind(kelp_data, kelp_data2)

##### Create Figure ############################################################
# Pallete 
group.colors <- c(Full = "#440154", Unprotected = "#FFBA00", Partial ="#21918c")

# Create base plot 
p <- ggplot(data = kelp_data_combo) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_boxplot(aes(x = mpa_status, y = mean_pr, fill = mpa_status), outlier.shape = NA, color = "grey30") +
  scale_fill_manual(values=group.colors) +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(y = "Relative Area", x = "") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

base <- facet(p, facet.by = c("region", "timeframe")) +
  theme(legend.position = "none") +
  stat_summary(aes(x = mpa_status, y = mean_pr), fun = "mean", 
               geom = "point", color = "black", pch = 21, fill = "white", size = 2.5)
base
# base <- ggboxplot(kelp_data_combo, 
#                   x = "mpa_status", 
#                   y = "mean_pr", 
#                   fill = "mpa_status",
#                   facet.by = c("region", "timeframe"),
#                   palette = group.colors, outlier.shape = NA) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
#                      limits = c(0, 3)) +
#   labs(y = "Relative Area", x = "") +
#   theme(legend.position = "none") 

dat_text <- data.frame(
  label = c("A", "D", "E", "B", "C", "F"),
  region   = c("All", "Central", "Southern"),
  timeframe = c("2014-2016", "2017-2021")
)

# Using the final results above make the stat.test tibble that the figure reads
stat.test <- results_long %>% 
  mutate(group1 = c("Full", "Partial", "Full", "Full", "Partial", "Full",
                    "Full", "Partial", "Full", "Full", "Partial", "Full", 
                    "Full", "Partial", "Full", "Full", "Partial", "Full")) %>% 
  mutate(group2 = c("Unprotected", "Unprotected", "Partial", "Unprotected", "Unprotected", "Partial", 
                    "Unprotected", "Unprotected", "Partial", "Unprotected", "Unprotected", "Partial",
                    "Unprotected", "Unprotected", "Partial", "Unprotected", "Unprotected", "Partial")) %>%
  mutate(p.adj = pvalues*6) %>% 
  select(timeframe, region, group1, group2, p.adj) %>% 
  mutate(y.position = c(4.25, 3.9, 3.9, 4.25, 3.9, 3.9,
                        4.25, 3.9, 3.9, 4.25, 3.9, 3.9, 
                        4.25, 3.9, 3.9, 4.25, 3.9, 3.9)) %>% 
  mutate("p.adj.signif" = c("***", "ns", "ns", "**", "ns", "ns",
                            "ns", "*", "ns", "ns", "ns", "ns", 
                            "***", "ns", "**","*", "ns", "**")) 
# create figure with significance bars 
boxplot_pr_w_sig <- base + 
  stat_pvalue_manual(stat.test, 
                     label = "p.adj.signif", 
                     tip.length = .00003, 
                     bracket.shorten = .05) +
  theme(strip.background =element_rect(fill="grey")) +
  geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -0.5,
    vjust   = -20,
    fontface = 2
  ) 

  

boxplot_pr_w_sig 

##### Export ###################################################################
png("Figures/Percent_recovery/pr_boxplot_w_sig.png", 
    width = 6, 
    height = 8, 
    units = "in", 
    res = 600)
boxplot_pr_w_sig 
dev.off() 

