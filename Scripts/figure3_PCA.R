# Date: July 25th 2023
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
data <- kelp_data_all %>% 
  filter(year == 2019) %>% 
  select(PixelID, mpa_status, hsmax, nitrate, temperature, MHW_intensity, CS_intensity,
         depth, gravity, distance_to_coast) 

data <- na.omit(data) # Removes 12% of the data 

data.pca <- prcomp(data[,-c(1,2)], center = TRUE, scale = TRUE)
summary(data.pca, loadings=TRUE)

screeplot(data.pca, type = "line", main = "Scree plot")

group.colors <- c(Full = "#440154", None = "#FFBA00", Partial ="#21918c")


plot_points <- fviz_pca_ind(data.pca, label="Mpa Category", habillage=data$mpa_status,
                            addEllipses=TRUE, ellipse.level=0.95, palette = group.colors, 
                            title = "Post Heatwave - 2019")
plot_points

#plot_points_2013 <- plot_points
#plot_points_2015 <- plot_points
#plot_points_2019 <- plot_points

combo_plot <- plot_grid(plot_points_2013, plot_points_2015, plot_points_2019, 
          ncol = 1, labels = c("A", "B", "C"))

##### Export ###################################################################
png("Figures/Environmental_PCA_points_span_hw.png", width = 5, height = 8, 
    units = "in", res = 600)
combo_plot
dev.off() 


