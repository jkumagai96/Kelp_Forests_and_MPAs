# Date: May 10th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Create Kelp area per year figures 
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(scales)
library(cowplot)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")
points_in_mpas <- read.csv("Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv") 

##### Explore MPA distributions ################################################
# Within our dataset, the max MPA establishment year is 2012 
# MPA area in 2021, how is it distributed?

kelp_data_all %>% 
  filter(year == 2021) %>% 
  group_by(region, mpa_status) %>% 
  count(mpa_status)

# So the number of pixels witin MPAs are basically the same, excpet there are more "None" pixels in
# Southern California

# Test actual kelp area 
mpas <- points_in_mpas %>% 
  select(PixelID, Site_ID_12, AreaMar_12, region) %>% 
  unique()
# 66 total MPAs 

mpas %>% 
  group_by(region) %>% 
  summarise(sum(AreaMar_12)) # there's 30% more mpa area in the souterh 

##### Kelp Area Figures ########################################################
# Kelp area through time 
df <- kelp_data_all %>% 
  group_by(year) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area),
            sum_kelp_area = sum(area)) 

# Mean kelp area through time 
ggplot(df, aes(x = year, y = mean_kelp_area)) + 
  geom_line(color = "green") +
  geom_point(color = "darkgreen") +
  theme_bw() 

# Median kelp area through time 
ggplot(df, aes(x = year, y = median_kelp_area)) + 
  geom_line(color = "orange") +
  geom_point(color = "darkorange") +
  theme_bw() 

# Total kelp area through time 
ggplot(df, aes(x = year, y = sum_kelp_area)) + 
  geom_line(color = "black") +
  geom_point(color = "grey50") +
  theme_bw() 

# Combo plot of mean and median through time
plot_all_sum <- ggplot(df) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 1997, xmax = 1998,
          ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "All", fontface = 2, 
           x = 2018, y = 55) + 
  geom_line(aes(x = year, y = sum_kelp_area/1e6), linewidth = 0.8) +
  scale_y_continuous(label = comma) + 
  theme_bw() +
  theme(legend.position = "none") + 
  ylab(bquote('Total Kelp Area '(km^2))) +
  xlab("Year")


# Kelp area through time split by region 
df2 <- kelp_data_all %>% 
  group_by(year, region) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area)) 

ggplot(df2, aes(x = year, y = mean_kelp_area, group = region)) + 
  geom_line(aes(linetype=region))+
  geom_point(aes(shape=region)) + 
  theme_bw() 

# Kelp area through time split by mpa status

df3 <- kelp_data_all %>% 
  group_by(year, mpa_status) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area)) 

group.colors <- c(Full = "#440154", None = "#FFBA00", Partial ="#21918c")

plot_all_split_median <- ggplot(df3, aes(x = year, y = median_kelp_area, group = mpa_status)) + 
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 1997, xmax = 1998,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "All", fontface = 2, 
           x = 2018, y = 45000) + 
  geom_line(aes(color=mpa_status))+
  scale_color_manual(values=group.colors, name = "MPA Category") +
  scale_y_continuous(label = comma) + 
  theme_bw() +
  theme(legend.position = "none") + 
  ylab(bquote('Median Kelp Area '(m^2))) +
  xlab("Year")

plot_all_split_avg <- ggplot(df3, aes(x = year, y = mean_kelp_area, group = mpa_status)) + 
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 1997, xmax = 1998,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "All", fontface = 2, 
           x = 2018, y = 105000) + 
  geom_line(aes(color=mpa_status), linewidth = 0.8)+
  scale_color_manual(values=group.colors, name = "MPA Category") +
  scale_y_continuous(label = comma) + 
  theme_bw() +
  theme(legend.position = "none") + 
  ylab(bquote('Mean Kelp Area '(m^2))) +
  xlab("Year")

# Just looking at southern california, kelp area through time split by mpa stauts
df4 <- kelp_data_all %>% 
  filter(region == "South_Coast") %>% 
  group_by(year, mpa_status) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area),
            sum_kelp_area = sum(area))

df4_b <- kelp_data_all %>% 
  filter(region == "South_Coast") %>% 
  group_by(year) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area),
            sum_kelp_area = sum(area))

plot_south_all_sum <- ggplot(df4_b, aes(x = year, y = sum_kelp_area/1e6)) + 
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 1997, xmax = 1998,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "Southern", fontface = 2, 
           x = 2018, y = 35) + 
  geom_line(linewidth = 0.8)+
  scale_y_continuous(label = comma) + 
  theme_bw() +
  ylab(bquote('Total Kelp Area '(km^2))) +
  xlab("Year")

plot_south_split_median <- ggplot(df4, aes(x = year, y = median_kelp_area, group = mpa_status)) + 
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 1997, xmax = 1998,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "Southern", fontface = 2, 
           x = 2018, y = 40000) + 
  geom_line(aes(color=mpa_status))+
  scale_color_manual(values=group.colors, name = "MPA Category") +
  scale_y_continuous(label = comma) + 
  theme_bw() +
  theme(legend.position = "none") + 
  ylab(bquote('Median Kelp Area '(m^2))) +
  xlab("Year")

plot_south_split_avg <- ggplot(df4, aes(x = year, y = mean_kelp_area, group = mpa_status)) + 
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 1997, xmax = 1998,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "Southern", fontface = 2, 
           x = 2018, y = 40000) + 
  geom_line(aes(color=mpa_status), linewidth = 0.8)+
  scale_color_manual(values=group.colors, name = "MPA Category") +
  scale_y_continuous(label = comma) + 
  theme_bw() +
  theme(legend.position = "none") + 
  ylab(bquote('Mean Kelp Area '(m^2))) +
  xlab("Year")

# Just looking at central california, kelp area through time split by mpa stauts
df5 <- kelp_data_all %>% 
  filter(region == "Central_Coast") %>% 
  group_by(year, mpa_status) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area)) 

df5_b <- kelp_data_all %>% 
  filter(region == "Central_Coast") %>% 
  group_by(year) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area),
            sum_kelp_area = sum(area))

plot_central_all_sum <- ggplot(df5_b, aes(x = year, y = sum_kelp_area/1e6)) + 
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 1997, xmax = 1998,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "Central", fontface = 2, 
           x = 2018, y = 35) + 
  geom_line(linewidth = 0.8)+
  scale_y_continuous(label = comma) + 
  theme_bw() +
  ylab(bquote('Total Kelp Area '(km^2))) +
  xlab("Year")

plot_central_split_median <- ggplot(df5, aes(x = year, y = median_kelp_area, group = mpa_status)) + 
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 1997, xmax = 1998,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "Central", fontface = 2, 
           x = 2018, y = 200000) + 
  geom_line(aes(color=mpa_status), linewidth = 0.8)+
  scale_color_manual(values=group.colors, name = "MPA Category") +
  scale_y_continuous(label = comma) + 
  theme_bw() +
  theme(legend.position = "none") + 
  ylab(bquote('Median Kelp Area '(m^2))) +
  xlab("Year")

plot_central_split_avg <- ggplot(df5, aes(x = year, y = mean_kelp_area, group = mpa_status)) + 
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 1997, xmax = 1998,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "Central", fontface = 2, 
           x = 2018, y = 200000) + 
  geom_line(aes(color=mpa_status), linewidth = 0.8)+
  scale_color_manual(values=group.colors, name = "MPA Category") +
  scale_y_continuous(label = comma) + 
  theme_bw() +
  theme(legend.position = "none") + 
  ylab(bquote('Mean Kelp Area '(m^2))) +
  xlab("Year")

plot_medians <- plot_grid(plot_all_sum, 
          plot_all_split_median, 
          plot_central_all_sum,
          plot_central_split_median,
          plot_south_all_sum,
          plot_south_split_median,
          ncol = 2,
          labels = c('A', 'B', 'C', 'D', 'E', 'F'), 
          label_size = 12)

plot_averages <- plot_grid(plot_all_sum, 
          plot_all_split_avg, 
          plot_central_all_sum,
          plot_central_split_avg,
          plot_south_all_sum,
          plot_south_split_avg,
          ncol = 2,
          labels = c('A', 'B', 'C', 'D', 'E', 'F'), 
          label_size = 12)

##### Kelp total and averages per category from 2012 ###########################
kelp_data_all %>% 
  filter(region == "Central_Coast") %>% 
  #filter(year >= 2012) %>% 
  group_by(year) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area),
            sum_kelp_area = sum(area)) %>% 
  ggplot(aes(x = year, y = sum_kelp_area/1e6)) + 
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "Central", fontface = 2, 
           x = 2018, y = 25) + 
  geom_line(linewidth = 0.8)+
  scale_y_continuous(label = comma) + 
  theme_bw() +
  ylab(bquote('Total Kelp Area '(km^2))) +
  xlab("Year")

kelp_data_all %>% 
  filter(region == "South_Coast") %>% 
  #filter(year >= 2012) %>% 
  group_by(year) %>% 
  summarise(mean_kelp_area = mean(area),
            median_kelp_area = median(area),
            sum_kelp_area = sum(area)) %>% 
  ggplot(aes(x = year, y = sum_kelp_area/1e6)) + 
  annotate("rect", fill = "red", alpha = 0.4, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  annotate("text", label = "Southern", fontface = 2, 
           x = 2018, y = 25) + 
  geom_line(linewidth = 0.8)+
  scale_y_continuous(label = comma) + 
  theme_bw() +
  ylab(bquote('Total Kelp Area '(km^2))) +
  xlab("Year")

##### Export ###################################################################
png("Figures/Kelp_area_medians.png", width = 10, height = 10, 
    units = "in", res = 600)
plot_medians
dev.off() 

png("Figures/Kelp_area_averages.png", width = 10, height = 10, 
    units = "in", res = 600)
plot_averages
dev.off() 
