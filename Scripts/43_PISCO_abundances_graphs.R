# Date: February 6th 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Visualizing PISCO data (abundance graphs)
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load packages
library(tidyverse)
library(cowplot)
library(png)
library(patchwork)

# Load Data
data <- read.csv("Processed_data/PISCO_data_summarized.csv")

# Declare Functions 
std <- function(x) sd(x)/sqrt(length(x))

# Declare group colors
group.colors <- c(Full = "#440154", Reference = "#FFBA00", Partial ="#21918c")

##### Load organism photos #####################################################
urchins <- readPNG("Figures/Organisms/Urchin_Figure_t.png",native = TRUE)
sheephead <- readPNG("Figures/Organisms/California_Sheephead.png", native = TRUE)
sheephead_3 <- readPNG("Figures/Organisms/Sheephead_3.png", native = TRUE)
spiny_lobster <- readPNG("Figures/Organisms/Spiny_Lobster_Figure_t.png", native = TRUE)

##### Format Data ##############################################################
# How many sites and how are they distributed across MPA grous?
unique(data$site) %>% length()

data %>% 
  distinct(site, mpa_status) %>% 
  count(mpa_status) %>% 
  mutate(p_n = 100*n/sum(n))

inverts_data <- data %>% filter(!is.na(urchin_d))
fish_data <- data %>% filter(!is.na(biomass_d))

##### Urchins through time #####################################################
Urchins_per_region <- inverts_data %>% 
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern")) %>% 
  group_by(year, mpa_status, region) %>% 
  summarise(total = mean(urchin_d),
            total_se = std(urchin_d)) %>% 
  ggplot(aes(x = year, y = total, color = mpa_status)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = total - total_se, 
                    ymax = total + total_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  facet_grid(cols = vars(region)) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab("Number of Urchnis") +
  ylab(bquote('Urchins per 60' ~ m^2)) +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = c(.1, .75),
        legend.background = element_rect(linewidth = 0.2, colour = 1),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Urchins_per_region

red_urchins <- inverts_data %>% 
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern")) %>% 
  group_by(year, mpa_status, region) %>% 
  summarise(Red = mean(MESFRAAD_d), 
            Red_se = std(MESFRAAD_d)) %>% 
  ggplot(aes(x = year, y = Red, color = mpa_status)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = Red - Red_se, 
                    ymax = Red + Red_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  facet_grid(vars(region)) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab(bquote('Red urchins per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
red_urchins

purple_urchins <- inverts_data %>% 
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern")) %>% 
  group_by(year, mpa_status, region) %>% 
  summarise(Purple = mean(STRPURAD_d), 
            Purple_se = std(STRPURAD_d)) %>% 
  ggplot(aes(x = year, y = Purple, color = mpa_status)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = Purple - Purple_se, 
                    ymax = Purple + Purple_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  facet_grid(vars(region)) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab(bquote('Purple urchins per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = c(.2, .85),
        legend.background = element_rect(linewidth = 0.2, colour = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

purple_urchins

combo_urchin_plot <- cowplot::plot_grid(purple_urchins, red_urchins)

##### Lobster through time #####################################################
lobster_plot <- inverts_data %>% 
  filter(region == "South_Coast") %>% 
  group_by(year, mpa_status) %>% 
  summarise(lobster = mean(PANINT_d), 
            lobster_se = std(PANINT_d)) %>% 
  ggplot(aes(x = year, y = lobster, color = mpa_status)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = lobster - lobster_se, 
                    ymax = lobster + lobster_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab(bquote('Spiny lobsters per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = c(.85, .85),
        legend.background = element_rect(linewidth = 0.2, colour = 1),
        panel.grid.minor = element_blank()) +
  inset_element(p = spiny_lobster,
                left = 0.03,
                bottom = 0.75,
                right = 0.2,
                top = 0.95) 
lobster_plot

##### Sheephead through time ###################################################
sheephead_plot <- fish_data %>% 
  filter(region == "South_Coast") %>% 
  group_by(year, mpa_status, region) %>% 
  summarise(sheepheads = mean(SPUL_d), 
            sheepheads_se = std(SPUL_d)) %>% 
  ggplot(aes(x = year, y = sheepheads, color = mpa_status)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = sheepheads - sheepheads_se, 
                    ymax = sheepheads + sheepheads_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab(bquote('California sheephead per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = "none",
    panel.grid.major = element_blank()) +
  inset_element(p = sheephead_3,
                left = 0.04,
                bottom = 0.7,
                right = 0.24,
                top = 0.95)

sheephead_plot


biomass_plot <- fish_data %>% 
  filter(region == "South_Coast") %>% 
  group_by(year, mpa_status, region) %>% 
  summarise(biomass = mean(biomass_d/1000), 
            biomass_se = std(biomass_d/1000)) %>% 
  ggplot(aes(x = year, y = biomass, color = mpa_status)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = biomass - biomass_se, 
                    ymax = biomass + biomass_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab(bquote('Sheephead biomass (kg) per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  inset_element(p = sheephead,
                left = 0.04,
                bottom = 0.7,
                right = 0.24,
                top = 0.95) 

biomass_plot

combo_plot <- plot_grid(lobster_plot, sheephead_plot, biomass_plot, 
                        labels = c("A", "B", "C"),
                        ncol = 1)
lobster_plot / sheephead_plot / biomass_plot

##### Export ###################################################################
png("Figures/Urchins_PISCO.png", width = 9, height = 6, 
    units = "in", res = 600)
combo_urchin_plot
dev.off()

png("Figures/Urchins_per_region_PISCO.png", width = 9, height = 5, 
    units = "in", res = 600)
Urchins_per_region
dev.off()

