# Date: August 10th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Produce sheephead, lobster, and urchin figures 
# BIO 202: Ecological Statistics

##### Load Data and Packages ###################################################
# load packages
library(tidyverse)
library(cowplot)
library(png)
library(patchwork)

# Load data
data <- read.csv("Processed_data/data_tables/subtidal_surveys.csv") 

##### Format Data ##############################################################
data <- data %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("None", "Partial", "Full")))

# Add in little black outlines of organisms
urchins <- readPNG("Figures/Organisms/Urchin_Figure_t.png",native = TRUE)
sheephead <- readPNG("Figures/Organisms/California_Sheephead.png", native = TRUE)
spiny_lobster <- readPNG("Figures/Organisms/Spiny_Lobster_Figure_t.png", native = TRUE)

##### Figures ##################################################################
group.colors <- c(Full = "#440154", None = "#FFBA00", Partial ="#21918c")

# Urchin Densities in Southern California 
plot_urchins <- data %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_urchins_all = mean(avg_urchins)) %>% 
  ggplot(aes(x = year, y = avg_urchins_all, color = mpa_status)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab(bquote('Urchins per' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  theme_bw() +
  theme(legend.position = c(.85, .7)) +
  theme(legend.background = element_rect(colour = 'black', 
                                         fill = 'white', 
                                         linetype='solid', 
                                         linewidth = 0.3)) +
  inset_element(p = urchins,
                             left = 0.1,
                             bottom = 0.7,
                             right = 0.24,
                             top = 0.95)

# Sheephead Densities 
plot_sheephead <- data %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_sheephead_all = mean(avg_sheephead)) %>% 
  ggplot(aes(x = year, y = avg_sheephead_all*60, color = mpa_status)) +
  geom_line(linewidth = 1, show.legend = FALSE) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab(bquote('CA Sheepheads per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  theme_bw() +
  inset_element(p = sheephead,
                            left = 0.1,
                            bottom = 0.7,
                            right = 0.24,
                            top = 0.95)
  

# Lobster Densities 
plot_lobsters <- data %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_lobster_all = mean(avg_lobster)) %>% 
  ggplot(aes(x = year, y = avg_lobster_all*60, color = mpa_status)) +
  geom_line(linewidth = 1, show.legend = FALSE) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab(bquote('CA Lobsters per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  theme_bw() +
  inset_element(p = spiny_lobster,
                left = 0.07,
                bottom = 0.7,
                right = 0.27,
                top = 0.95)

combo_plot <- plot_grid(plot_urchins, plot_sheephead, plot_lobsters,
                        labels = "AUTO", nrow = 3)

##### Export ###################################################################
png("Figures/Urchins_sheephead_lobsters.png", width = 6, height = 9, 
    units = "in", res = 600)
combo_plot
dev.off() 



