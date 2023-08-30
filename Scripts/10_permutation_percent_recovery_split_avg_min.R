# Date: July 8th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Bootstrap - Permutation approach based on difference in percent recovery instead of kelp area 
# Spitting the data by region instead of treating it all the same
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(sf)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")
points_in_mpas <- read.csv("Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv") 

##### Declare Functions ########################################################
Calculate_percent_recovery <- function(x, min_v, max_v) {
  value <- 100*(x - min_v)/(max_v - min_v)
  return(value)
}

##### Format Data ##############################################################
data_area_long <- kelp_data_all %>% 
  filter(year < 2014) %>% 
  dplyr::select(PixelID, year, area, min_area) 

maxes <- data_area_long %>% 
  dplyr::select(-min_area) %>% 
  pivot_wider(names_from = year, values_from = area) 

mins <- kelp_data_all %>% 
  filter(year >= 2014 & year < 2017) %>% 
  dplyr::select(PixelID, year, min_area) %>% 
  pivot_wider(names_from = year, values_from = min_area) 

# Create new dataframe to store values in 
df <- maxes %>% 
  dplyr::select(PixelID)

# Calculate Max and Min value 
df$min <- apply(mins, 1, FUN = mean, na.rm = T) # Changed to mean of min's 
df$max <- apply(maxes, 1, FUN = mean, na.rm = T)

# Check to see whether the mean of maxes is lower or equal to the "min" value
df[df$max <= df$min,] # 118 records

df$remove <- df$max <= df$min

##### Calculate Percent Recovery ###############################################
# Calculate Percent Recovery for each row in the dataset from 2014 to 2021
df_percent_recovery <- kelp_data_all %>% 
  filter(year >= 2014) %>% 
  dplyr::select(PixelID, year, area) %>% 
  left_join(df, by = "PixelID") %>% 
  mutate(percent_recovery = NA)

for (i in 1:nrow(df_percent_recovery )) {
  area_i <- df_percent_recovery$area[i]
  min_i <- df_percent_recovery$min[i]
  max_i <- df_percent_recovery$max[i]
  
  df_percent_recovery$percent_recovery[i] <- Calculate_percent_recovery(x = area_i,
                                                                        min_v = min_i,
                                                                        max_v = max_i)
}

# Join data back 
kelp_data <- kelp_data_all %>% 
  filter(year >= 2014) %>% 
  select(PixelID, year, mpa_status, area, region) %>% 
  left_join(df_percent_recovery, by = c("PixelID", "year", "area")) %>% 
  filter(remove == FALSE)

##### Structure Data ###########################################################
kelp_data_south <- kelp_data %>% 
  filter(region == "South_Coast")

kelp_data_central <- kelp_data %>% 
  filter(region == "Central_Coast")

# Add region onto the points dataset
a <- kelp_data %>% select(region, PixelID) %>% distinct()

points_in_mpas_south <- points_in_mpas %>% 
  #left_join(a, by = "PixelID") %>% 
  filter(region == "South_Coast")

points_in_mpas_central <- points_in_mpas %>% 
  #left_join(a, by = "PixelID") %>% 
  filter(region == "Central_Coast")

##### Calculate True Values for south region ###################################
# Change global option to not print out message about group summaries 
options(dplyr.summarise.inform = FALSE)

# Calculate average kelp area per category per year 
true_values <- kelp_data_south %>% 
  group_by(year, mpa_status) %>% 
  summarise(median_pr = median(percent_recovery)) %>% 
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
years <- 2014:2021

# Remove original mpa_status from kelp_data
kelp_data_r <- kelp_data_south %>% select(-mpa_status)

all_pixels <- kelp_data_south$PixelID %>% unique()

# within the for loop 
start <- Sys.time()
for (j in 1:10000) {
  
  # Sample random pixels that will overlap w/ MPAs 
  r_pixels <- sample(all_pixels, size = nrow(points_in_mpas_south), replace = F) # number of pixels originally overlapping with MPAs 
  
  # Assign the random pixels to the points in the mpas instead of the original values
  points_in_mpas_south$PixelID <- r_pixels
  
  # Join randomized mpa data with kelp data 
  kelp_w_mpas_r <- left_join(kelp_data_r, points_in_mpas_south, by = "PixelID") 
  
  # Ensure that the mpa status accounts for the year of establishment
  kelp_w_mpas_r$mpa_status[is.na(kelp_w_mpas_r$mpa_status)] <- "None"
  
  not_protected <- kelp_w_mpas_r %>% 
    filter(mpa_status == "None") %>% 
    arrange(PixelID)
  
  subset <- kelp_w_mpas_r %>% 
    filter(mpa_status != "None")
  
  for (i in 1:nrow(subset)) {
    status <- subset$mpa_status[i]
    establishment_yr <- subset$Estab_Yr_1[i]
    
    if (subset$year[i] >= establishment_yr) {
      subset$mpa_status[i] <- status
    } else {
      subset$mpa_status[i] <- "None"
      subset$AreaMar_12[i] <- 0
    }
  }
  
  subset <- subset %>% arrange(PixelID)
  new_data <- rbind(not_protected, subset)
  
  # Calculate table 
  values <- new_data %>% 
    group_by(year, mpa_status) %>% 
    summarise(median_pr = median(percent_recovery)) %>% 
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
final_results_south <- data.frame()

for (i in 1:length(years)) {
  
  # Filter bootstrap values by year
  df <- bootstrap_df %>% 
    filter(year == years[i]) 
  
  tv <- true_values %>% 
    filter(year == years[i])
  
  sig_F_N <- 1 - ecdf(df$F_N)(tv$F_N)
  sig_P_N <- 1 - ecdf(df$P_N)(tv$P_N)
  sig_F_P <- 1 - ecdf(df$F_P)(tv$F_P)
  
  t <- data.frame(years[i],
                  sig_F_N,
                  sig_P_N,
                  sig_F_P)
  
  colnames(t) <- c("year", "F_N", "P_N", "F_P")
  
  final_results_south <- rbind(final_results_south, t)
  
}


##### Export ###################################################################
final_results_south

write.csv(final_results_south, 
          "Processed_data/data_tables/permutation_sensitivity/bootstrap_avg_min_south.csv", 
          row.names = F)


rm(kelp_data_south)
rm(true_values, kelp_data_r, all_pixels, bootstrap_list)

##### Calculate True Values for central region #################################
# Calculate average kelp area per category per year 
true_values <- kelp_data_central %>% 
  group_by(year, mpa_status) %>% 
  summarise(median_pr = median(percent_recovery)) %>% 
  pivot_wider(names_from = mpa_status,
              values_from = median_pr) %>% 
  mutate(P_N = Partial - None,
         F_N = Full - None,
         F_P = Full - Partial) %>% 
  select(-c(None, Partial, Full))

## Set up:
# Set up variables outside of the loops
bootstrap_list <- list()
years <- 2014:2021

# Remove original mpa_status from kelp_data
kelp_data_r <- kelp_data_central %>% select(-mpa_status)

all_pixels <- kelp_data_central$PixelID %>% unique()

# within the for loop 
start <- Sys.time()
for (j in 1:10000) {
  
  # Sample random pixels that will overlap w/ MPAs 
  r_pixels <- sample(all_pixels, size = nrow(points_in_mpas_central), replace = F) # number of pixels originally overlapping with MPAs 
  
  # Assign the random pixels to the points in the mpas instead of the original values
  points_in_mpas_central$PixelID <- r_pixels
  
  # Join randomized mpa data with kelp data 
  kelp_w_mpas_r <- left_join(kelp_data_r, points_in_mpas_central, by = "PixelID") 
  
  # Ensure that the mpa status accounts for the year of establishment
  kelp_w_mpas_r$mpa_status[is.na(kelp_w_mpas_r$mpa_status)] <- "None"
  
  not_protected <- kelp_w_mpas_r %>% 
    filter(mpa_status == "None") %>% 
    arrange(PixelID)
  
  subset <- kelp_w_mpas_r %>% 
    filter(mpa_status != "None")
  
  for (i in 1:nrow(subset)) {
    status <- subset$mpa_status[i]
    establishment_yr <- subset$Estab_Yr_1[i]
    
    if (subset$year[i] >= establishment_yr) {
      subset$mpa_status[i] <- status
    } else {
      subset$mpa_status[i] <- "None"
      subset$AreaMar_12[i] <- 0
    }
  }
  
  subset <- subset %>% arrange(PixelID)
  new_data <- rbind(not_protected, subset)
  
  # Calculate table 
  values <- new_data %>% 
    group_by(year, mpa_status) %>% 
    summarise(median_pr = median(percent_recovery)) %>% 
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
final_results_central <- data.frame()

for (i in 1:length(years)) {
  
  # Filter bootstrap values by year
  df <- bootstrap_df %>% 
    filter(year == years[i]) 
  
  tv <- true_values %>% 
    filter(year == years[i])
  
  sig_F_N <- 1 - ecdf(df$F_N)(tv$F_N)
  sig_P_N <- 1 - ecdf(df$P_N)(tv$P_N)
  sig_F_P <- 1 - ecdf(df$F_P)(tv$F_P)
  
  t <- data.frame(years[i],
                  sig_F_N,
                  sig_P_N,
                  sig_F_P)
  
  colnames(t) <- c("year", "F_N", "P_N", "F_P")
  
  final_results_central <- rbind(final_results_central, t)
  
}

final_results_central
write.csv(final_results_central, 
          "Processed_data/data_tables/permutation_sensitivity/bootstrap_avg_min_central.csv", 
          row.names = F)


##### Figures ##################################################################
##### Plot for southern california 
final_results_south <- read.csv("Processed_data/data_tables/permutation_sensitivity/bootstrap_avg_min_south.csv")

# Format Data into long format
results_long <- final_results_south %>% 
  pivot_longer(names_to = "Comparison", values_to = "pvalues", F_N:F_P) 

results_long$pvalues[results_long$pvalues == 0] <- .000001

# Plot 
plot1 <- results_long %>% 
  ggplot(aes(x = year, y = -log10(pvalues), group = Comparison)) +
  geom_point(aes(color = Comparison, shape = Comparison), size = 2) +
  geom_hline(yintercept = -log10(0.05/(18))) +
  geom_hline(yintercept = -log10(0.05/3), linetype = "dashed") +
  scale_color_manual(values=c('#FF5C00', '#999999','#000EDD')) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021)) +
  labs(y = "-log10(P values)", x = "Year") 

# Export 
png("Figures/Permutation_sensitivity/Bootstrap_avg_min_south.png", width = 5, height = 3, units = "in", res = 600)
plot1
dev.off() 


##### Plot for central california
final_results_central <- read.csv("Processed_data/data_tables/permutation_sensitivity/bootstrap_avg_min_central.csv")

# Format Data into long format
results_long <- final_results_central %>% 
  pivot_longer(names_to = "Comparison", values_to = "pvalues", F_N:F_P) 

results_long$pvalues[results_long$pvalues == 0] <- .000001

# Plot 
plot2 <- results_long %>% 
  ggplot(aes(x = year, y = -log10(pvalues), group = Comparison)) +
  geom_point(aes(color = Comparison, shape = Comparison), size = 2) +
  geom_hline(yintercept = -log10(0.05/(18))) +
  geom_hline(yintercept = -log10(0.05/3), linetype = "dashed") +
  scale_color_manual(values=c('#FF5C00', '#999999','#000EDD')) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylim(0, 6) +
  scale_x_continuous(breaks = c(2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021)) +
  labs(y = "-log10(P values)", x = "Year") 

# Export
png("Figures/Permutation_sensitivity/Bootstrap_avg_min_central.png", width = 5, height = 3, units = "in", res = 600)
plot2
dev.off() 


#### change figures 
