# Date: May 22nd 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: PCA analysis of environmental variables to determine whether the 
# environments of kelp within protected areas are experience are different than 
# the environments of kelp within non protected areas
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(factoextra)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")

##### Format Data ##############################################################
# Only for 2014, and then 2016, 2017, 2019, 2021?

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
             addEllipses=TRUE, ellipse.level=0.95, palette = group.colors)
plot_points
fviz_pca_var(data.pca, col.var = "steelblue")

plot_biplot <- fviz_pca_biplot(data.pca, habillage=data$mpa_status, label = "var",
                palette = group.colors,
                ggtheme = theme_minimal())

##### Export ###################################################################
png("Figures/Environmental_PCA_points_2019.png", width = 7, height = 5, 
    units = "in", res = 600)
plot_points
dev.off() 

png("Figures/Environmental_PCA_biplot_2019.png", width = 7, height = 5, 
    units = "in", res = 600)
plot_biplot 
dev.off() 

##### Explore MPA age in the PCA! 



print("Script is finished")