# Date: May 10th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Bootstrap - Permutation approach based on difference in percent recovery 
#          Looking at the time period from 1984 to 2002
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(sf)

# Load Data
kelp_raw <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")
points_in_mpas <- read.csv("Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv") 

##### Declare Functions ########################################################
Calculate_percent_recovery <- function(x, min_v, max_v) {
  value <- 100*(x - min_v)/(max_v - min_v)
  return(value)
}

##### Remove Pixels with MPAs established before 2003 ##########################
points_in_mpas_alt <- points_in_mpas %>% 
  select(PixelID, Estab_Yr_1) %>% 
  filter(is.na(Estab_Yr_1) == TRUE | Estab_Yr_1 > 2002)

points_in_mpas_loop <- points_in_mpas %>% 
  select(PixelID, Estab_Yr_1, mpa_status) %>% 
  filter(is.na(Estab_Yr_1) == TRUE | Estab_Yr_1 > 2002)

kelp_data_all <- kelp_raw %>% 
  left_join(points_in_mpas_alt, by = "PixelID") %>% 
  dplyr::select(-Estab_Yr_1) %>% 
  filter(year < 2003) # 2002 is the last year so we are analyzing 1999, 2000, 2001, and 2002 recovery data

##### Format Data ##############################################################
maxes <- kelp_data_all %>% 
  filter(year < 1997) %>% 
  dplyr::select(PixelID, year, area) %>% 
  pivot_wider(names_from = year, values_from = area) 

mins <- kelp_data_all %>% 
  filter(year >= 1997 & year < 1999) %>% 
  dplyr::select(PixelID, year, min_area) %>% 
  pivot_wider(names_from = year, values_from = min_area) 

# Create new dataframe to store values in 
df <- maxes %>% 
  dplyr::select(PixelID)

# Calculate Max and Min value 
mins$mins <- apply(mins, 1, FUN = min, na.rm = T) 
df$max <- apply(maxes, 1, FUN = mean, na.rm = T)

df <- df %>% left_join(mins[,c(1,4)], by = "PixelID")
df <- na.omit(df) # 2% of the data removed where there are no mins for the year (not enough data that year I guess)

##### Calculate Percent Recovery ###############################################
# Calculate Percent Recovery for each row in the dataset from 1997 to 2002
df_percent_recovery <- kelp_data_all %>% 
  filter(year >= 1997) %>% 
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

#### Assign MPA status of the future (2021) to 1984 to 2002 ###
df_mpa_staus_2021 <- kelp_raw %>% 
  filter(year == 2021) %>% 
  dplyr::select(PixelID, mpa_status)

# Join data back 
kelp_data <- kelp_data_all %>% 
  filter(year >= 1997) %>% 
  select(PixelID, year, area) %>% 
  left_join(df_percent_recovery, by = c("PixelID", "year", "area")) %>% 
  left_join(df_mpa_staus_2021, by = "PixelID") 

##### Calculate True Values ####################################################
# Start here, I need to understand how the mpa categories are assigned here.
# Need to check the points in mpas to make sure points with an establishent year before 2002 are not included
# Then, I need to make sure that I assign the mpa categories correctly and they stay consistent through the iteration 


# Change global option to not print out message about group summaries 
options(dplyr.summarise.inform = FALSE)

# Calculate median percent kelp recovery per category per year for the real data
true_values <- kelp_data %>% 
  group_by(year, mpa_status) %>% 
  summarise(median_pr = median(percent_recovery, na.rm = T)) %>% 
  pivot_wider(names_from = mpa_status,
              values_from = median_pr) %>% 
  mutate(P_N = Partial - None,
         F_N = Full - None,
         F_P = Full - Partial) %>% 
  select(-c(None, Partial, Full))

##### Bootstrapping ############################################################
## Set up:
set.seed(456) # So the results are repeatable

# Set up variables outside of the loops
years <- 1997:2002

# Remove original mpa_status from kelp_data
kelp_data_r <- kelp_data %>% select(-mpa_status)

all_pixels <- kelp_data$PixelID %>% unique()

bootstrap_list <- list()
# within the for loop 
start <- Sys.time()
for (j in 1:10000) {
  
  # Sample random pixels that will overlap w/ MPAs 
  r_pixels <- sample(all_pixels, size = nrow(points_in_mpas_loop), replace = F) # number of pixels originally overlapping with MPAs 
  
  # Assign the random pixels to the points in the mpas instead of the original values
  points_in_mpas_loop$PixelID <- r_pixels
  
  # Join randomized mpa data with kelp data 
  kelp_w_mpas_r <- dplyr::left_join(kelp_data_r, points_in_mpas_loop, by = "PixelID") 
  
  # Ensure that the mpa status accounts for the year of establishment
  kelp_w_mpas_r$mpa_status[is.na(kelp_w_mpas_r$mpa_status)] <- "None"
  
  # Calculate table, changes here too!
  values <- kelp_w_mpas_r %>% 
    group_by(year, mpa_status) %>% 
    summarise(median_pr = median(percent_recovery, na.rm = T)) %>% 
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

#### Processing Results ########################################################
bootstrap_df <- do.call(rbind.data.frame, bootstrap_list)

# Each loop will be a row in the new dataframe
final_results <- data.frame()

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
  
  final_results <- rbind(final_results, t)
  
}

##### Export ###################################################################
final_results

write.csv(final_results, 
          "Processed_data/data_tables/bootstrap_previous_heatwave_analysis.csv", 
          row.names = F)

print('bootstrap results saved')
##### Figures ##################################################################
# Format Data into long format
results_long <- final_results %>% 
  pivot_longer(names_to = "Comparison", values_to = "pvalues", F_N:F_P) 

results_long$pvalues[results_long$pvalues == 0] <- .000001

print('now attempting to make figure')
# Plot 
plot1 <- results_long %>% 
  ggplot(aes(x = year, y = -log10(pvalues), group = Comparison)) +
  geom_point(aes(color = Comparison, shape = Comparison), size = 2) +
  geom_hline(yintercept = -log10(0.05/(3*6))) +
  geom_hline(yintercept = -log10(0.05/3), linetype = "dashed") +
  scale_color_manual(values=c('#FF5C00', '#999999','#000EDD')) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(1997, 1998, 1999, 2000, 2001, 2002)) +
  labs(y = "-log10(P values)", x = "Year") 

# Export
print('attempting export')

png("Figures/bootstrap_previous_heatwave_analysis.png", width = 5, height = 3, units = "in", res = 600)
plot1
dev.off() 

