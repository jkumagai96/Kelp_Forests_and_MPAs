# Date: Feb. 24th 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Create Kelp area per year figures 
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(scales)
library(cowplot)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv") %>% 
  mutate(mpa_status = ifelse(mpa_status == "None", "Unprotected", mpa_status))
points_in_mpas <- read.csv("Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv") 

##### Format Data ##############################################################
maxes <- kelp_data_all %>% 
  filter(year < 2014) %>% 
  dplyr::select(PixelID, year, region, area) %>% 
  group_by(PixelID, region) %>% 
  summarize(historic_baseline = mean(area)) %>% 
  group_by(region) %>% 
  summarize(baseline = sum(historic_baseline, na.rm = T)) 

# On average, how much more kelp was within central vs. southern during and after the heatwave
kelp_data_all %>% 
  filter(year > 2013) %>% 
  group_by(year, region) %>% 
  summarize(mean_area = mean(area)) %>% 
  pivot_wider(names_from = region, values_from = mean_area) %>% 
  mutate(Difference = Central_Coast/South_Coast) %>% view()

##### Total kelp per region #####
plot_central_total <- kelp_data_all %>% 
  filter(region == "Central_Coast") %>% 
  group_by(year) %>% 
  summarise(sum_kelp_area = sum(area)) %>% 
  ggplot(aes(x = year, y = sum_kelp_area/1e6)) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 1997, xmax = 1998,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "Central", fontface = 2, 
           x = 1984.9, y = 35) +
  geom_hline(yintercept = maxes$baseline[1]/1e6, color = "grey40", lty = "dashed", linewidth = 1) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(label = comma) + 
  theme_bw() +
  theme(legend.position = "none") + 
  ylab(bquote('Total Kelp Area '(km^2))) +
  xlab("Year")

plot_south_total <- kelp_data_all %>% 
  filter(region == "South_Coast") %>% 
  group_by(year) %>% 
  summarise(sum_kelp_area = sum(area)) %>% 
  ggplot(aes(x = year, y = sum_kelp_area/1e6)) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 1997, xmax = 1998,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "Southern", fontface = 2, 
           x = 1985.7, y = 35) + 
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = maxes$baseline[2]/1e6, color = "grey40", lty = "dashed", linewidth = 1) +
  scale_y_continuous(label = comma) + 
  theme_bw() +
  theme(legend.position = "none") + 
  ylab(bquote('Total Kelp Area '(km^2))) +
  xlab("Year")



##### Average kelp per mpa status ######
group.colors <- c(Full = "#440154", Unprotected = "#FFBA00", Partial ="#21918c")

plot_central_mpa <- kelp_data_all %>% 
  filter(region == "Central_Coast") %>% 
  filter(year >= 2012) %>% 
  group_by(year, mpa_status) %>% 
  summarise(mean_kelp_area = mean(area)) %>% 
  ggplot(aes(x = year, y = mean_kelp_area, group = mpa_status)) + 
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "Central", fontface = 2, 
           x = 2012.3, y = 79000) + 
  geom_line(aes(color=mpa_status), linewidth = 0.8)+
  scale_color_manual(values=group.colors, name = "Protection Status") +
  scale_y_continuous(label = comma) + 
  theme_bw() +
  theme(legend.position = c(.84,.8), legend.title = element_text(size = 10)) + 
  ylab(bquote('Mean Kelp Area '(m^2))) +
  xlab("Year")

plot_south_mpa <- kelp_data_all %>% 
  filter(region == "South_Coast") %>% 
  filter(year >= 2012) %>% 
  group_by(year, mpa_status) %>% 
  summarise(mean_kelp_area = mean(area)) %>% 
  ggplot(aes(x = year, y = mean_kelp_area, group = mpa_status)) + 
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "Southern", fontface = 2, 
           x = 2012.5, y = 29000) + 
  geom_line(aes(color=mpa_status), linewidth = 0.8)+
  scale_color_manual(values=group.colors, name = "Protection Status") +
  scale_y_continuous(label = comma) + 
  theme_bw() +
  theme(legend.position = "none") + 
  ylab(bquote('Mean Kelp Area '(m^2))) +
  xlab("Year")


##### Export 
combo_plot <- cowplot::plot_grid(plot_central_total, 
                   plot_central_mpa,
                   plot_south_total, 
                   plot_south_mpa, labels = "AUTO")

png("Figures/Kelp_area_central_south.png", width = 10, height = 7, 
    units = "in", res = 600)
combo_plot
dev.off() 

