# Date: December 13th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: PCA analysis of environmental variables pre, during, and post heatwave - Figure 3
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(factoextra)
library(cowplot)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")

##### Format Data ##############################################################
min_non_zero <-  sort(unique(kelp_data_all$gravity)) %>% head(2)  
min_non_zero <- min_non_zero[2]

kelp_data_all$gravity[kelp_data_all$gravity == 0] <- min_non_zero

# 2013
data <- kelp_data_all %>% 
  filter(year == 2013) %>%
  mutate(log_gravity = log(gravity)) %>% 
  select(PixelID, mpa_status, region, hsmax, nitrate, temperature, MHW_intensity, CS_intensity,
         depth, log_gravity) %>% 
  rename(temp = temperature, 
         MHW = MHW_intensity, 
         CS = CS_intensity, 
         gravity = log_gravity)

data.pca <- prcomp(data[,-c(1,2,3)], center = TRUE, scale = TRUE)
summary(data.pca)

screeplot(data.pca, type = "line", main = "Scree plot")

group.colors <- c(Full = "#440154", None = "#FFBA00", Partial ="#21918c")
region.colors <- c(Central_Coast = "#BF2C23", South_Coast = "#2F67B1")

plot_points_2013 <- fviz_pca_ind(data.pca, habillage=data$mpa_status, label = "var",
                palette = group.colors,
                ggtheme = theme_minimal(), 
                addEllipses = TRUE, 
                ellipse.level = 0.95,
                title = "Pre Heatwave - 2013") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

data %>% 
  mutate(PC1 = data.pca$x[,1], PC2 = data.pca$x[,2]) %>% 
  ggplot(aes(PC1, PC2, color = mpa_status)) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = region)) + 
  stat_ellipse(aes(fill = mpa_status), geom = "polygon", alpha = 0.1) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = group.colors) + 
  scale_fill_manual(values = group.colors) + 
  theme_minimal()


plot_points_2013$data %>% 
  ggplot(aes(x, y, color = Groups))

plot_points_2013

region_plot_2013 <- fviz_pca_ind(data.pca, habillage=data$region, label = "var",
                               palette = region.colors,
                               ggtheme = theme_minimal(), 
                               addEllipses = TRUE, 
                               ellipse.level = 0.95,
                               title = " ") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

region_plot_2013

var_plot_2013 <- fviz_pca_var(data.pca, 
                              col.var = "black", 
                              col.circle = "white", title = " ") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())





# 2015
data <- kelp_data_all %>% 
  filter(year == 2015) %>%
  mutate(log_gravity = log(gravity)) %>% 
  select(PixelID, mpa_status, region, hsmax, nitrate, temperature, MHW_intensity, CS_intensity,
         depth, log_gravity) %>% 
  rename(temp = temperature, 
         MHW = MHW_intensity, 
         CS = CS_intensity, 
         gravity = log_gravity) #%>% 
  #select(-CS)

data.pca <- prcomp(data[,-c(1,2,3)], center = TRUE, scale = TRUE)
summary(data.pca)

screeplot(data.pca, type = "line", main = "Scree plot")

group.colors <- c(Full = "#440154", None = "#FFBA00", Partial ="#21918c")
region.colors <- c(Central_Coast = "#BF2C23", South_Coast = "#2F67B1")

plot_points_2015 <- fviz_pca_ind(data.pca, habillage=data$mpa_status, label = "var",
                                    palette = group.colors,
                                    ggtheme = theme_minimal(), 
                                    addEllipses = TRUE, 
                                    ellipse.level = 0.95,
                                    title = "Pre Heatwave - 2015") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_points_2015

region_plot_2015 <- fviz_pca_ind(data.pca, habillage=data$region, label = "var",
                                    palette = region.colors,
                                    ggtheme = theme_minimal(), 
                                    addEllipses = TRUE, 
                                    ellipse.level = 0.95,
                                 title = " ") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

region_plot_2015

var_plot_2015 <- fviz_pca_var(data.pca, 
                              col.var = "black", 
                              col.circle = "white", title = " ") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# 2019
data <- kelp_data_all %>% 
  filter(year == 2019) %>%
  mutate(log_gravity = log(gravity)) %>% 
  select(PixelID, mpa_status, region, hsmax, nitrate, temperature, MHW_intensity, CS_intensity,
         depth, log_gravity) %>% 
  rename(temp = temperature, 
         MHW = MHW_intensity, 
         CS = CS_intensity, 
         gravity = log_gravity)

data.pca <- prcomp(data[,-c(1,2,3)], center = TRUE, scale = TRUE)
summary(data.pca)

screeplot(data.pca, type = "line", main = "Scree plot")

group.colors <- c(Full = "#440154", None = "#FFBA00", Partial ="#21918c")
region.colors <- c(Central_Coast = "#BF2C23", South_Coast = "#2F67B1")

plot_points_2019 <- fviz_pca_ind(data.pca, habillage=data$mpa_status, label = "var",
                                    palette = group.colors,
                                    ggtheme = theme_minimal(), 
                                    addEllipses = TRUE, 
                                    ellipse.level = 0.95,
                                    title = "Post Heatwave - 2019") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_points_2019

region_plot_2019 <- fviz_pca_ind(data.pca, habillage=data$region, label = "var",
                                    palette = region.colors,
                                    ggtheme = theme_minimal(), 
                                    addEllipses = TRUE, 
                                    ellipse.level = 0.95, 
                                 title = " ") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

region_plot_2019

var_plot_2019 <- fviz_pca_var(data.pca, 
                              col.var = "black", 
                              col.circle = "white", title = " ") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##### Combine plots ############################################################
leg <- get_legend(plot_points_2019)
leg2 <- get_legend(region_plot_2019)

combo_plot <- plot_grid(plot_points_2013,
                        plot_points_2015,
                        plot_points_2019 + theme(legend.position = "None"),
                        leg,
                        rel_heights = c(1, 1, 1, 0.2),
          ncol = 1,
          labels = c("A", "D", "G"))
combo_plot

combo_plot2 <- plot_grid(region_plot_2013,
                         region_plot_2015,
                         region_plot_2019 + theme(legend.position = "None"),
                        leg2,
                        rel_heights = c(1, 1, 1, 0.2),
                        ncol = 1,
                        labels = c("B", "E", "H"))
combo_plot2


combo_plot3 <- plot_grid(var_plot_2013,
                         var_plot_2015,
                         var_plot_2019,
                         leg2,
                         rel_heights = c(1, 1, 1, 0.2),
                         ncol = 1,
                         labels = c("C", "F", "I"))
combo_plot3

final_plot <- plot_grid(combo_plot, combo_plot2, combo_plot3, ncol = 3)
final_plot
##### Export ###################################################################
png("Figures/Environmental_PCA_points_span_hw.png", width = 12, height = 12, 
    units = "in", res = 600)
final_plot
dev.off() 


