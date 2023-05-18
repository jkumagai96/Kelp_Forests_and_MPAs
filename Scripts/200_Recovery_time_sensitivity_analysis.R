# Date: May 17th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Try to calculate recovery time 
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(sf)
library(viridis)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")

##### Calculate baseline area per pixel ########################################
baseline_data <- kelp_data_all %>% 
  dplyr::select(PixelID, year, area) %>% 
  filter(year < 2014) %>% # 30 year baseline starting in 1984 to 2013
  group_by(PixelID) %>% 
  summarize(mean_area = mean(area, na.rm = T),
            sd_area = sd(area, na.rm = T), #)
            one_sd_low = mean_area - sd_area) #%>% 
#filter(one_sd_low > 0)

# Create a dataset with the post heat wave data (2017-2021)
post_hw_data <- kelp_data_all %>% 
  dplyr::select(PixelID, year, area) %>% 
  filter(year > 2016)

##### Remove lowest Pixels #####################################################
# ADJUST THIS SECTION!
total_area <- sum(post_hw_data$area) 

maxes <- kelp_data_all %>% 
  filter(year < 2014) %>% 
  dplyr::select(PixelID, year, area) %>%  
  pivot_wider(names_from = year, values_from = area) 

# Create new dataframe to store values in 
df <- maxes %>% 
  dplyr::select(PixelID)

# Calculate mean of maxes value 
df$max <- apply(maxes, 1, FUN = mean, na.rm = T)

cutoff_percent <- .5
cutoff_area <- quantile(df$max, cutoff_percent) # Based on historical values 

cutoff_table <- df %>% 
  filter(max >= cutoff_area) %>% 
  mutate(keep = 1)

post_hw_data <- post_hw_data %>% 
  left_join(cutoff_table) %>% 
  filter(keep == 1)

new_area <- sum(post_hw_data$area) 

100 - 100*(new_area/total_area) # Percent of total removed from the dataset

##### Load Functions ###########################################################
# Small function to determine if a value is within range, used in the next function
# to see if the kelp area value is within the baseline range. 
is_in_range <- function(x, lower, upper) {
  x >= lower & x <= upper
}

# This function calculates the "recovery time" index which basically is the 
# number of years after the heat wave that the kelp area returns to within it's normal 
# range. Set my the standard deviation limit (sd_limit).

# Baseline data is a 30 year baseline (1984 to 2013) to determine the mean and variance of kelp area per pixel
# Post_hw_data is the data after the heat wave (2016-2021) where you have a value for kelp area per year per pixel 
# n_repeats is the number of years where kelp area is back within the normal range to count that it is recovered 
Calculate_recovery_time <- function(baseline_data, post_hw_data, sd_limit, n_repeats) {
  
  # Create table to store resutls, number of years to recover per pixel
  recovery_time_table <- data.frame(PixelID = integer(), Recovery_time = integer())
  
  # All unique pixels 
  pixels <- unique(post_hw_data$PixelID)
  
  # All unique years 
  post_hw_years <- unique(post_hw_data$year)
  
  # Loop which goes through each pixel 
  for (ii in 1:length(pixels)) {
    
    # Post heat wave data per pixel 
    data_per_pixel <- post_hw_data %>% 
      filter(PixelID == pixels[ii])
    
    # Baseline data per pixel 
    baseline_data_per_pixel <- baseline_data %>% 
      filter(PixelID == pixels[ii])
    
    # Create a column that will be filled with TRUE or FALSE in the next step 
    data_per_pixel$within_range <- NA
    
    # For loop that calcualtes for each year, whether the kelp area is within range 
    for (i in 1:length(post_hw_years)) {
      kelp_area <- data_per_pixel$area[i]
      mean_area <- baseline_data_per_pixel$mean_area
      sd_area <- baseline_data_per_pixel$sd_area
      
      data_per_pixel$within_range[i] <- is_in_range(kelp_area, mean_area - (sd_limit*sd_area), mean_area + (sd_limit*sd_area))
    }
    
    # The next two lines determine at which position there are multiple TRUE's (the kelp area is within range)
    test_range <- data_per_pixel$within_range
    repeated_ones_pos <- which(stats::filter(test_range, rep(1, n_repeats), sides=1) == n_repeats)
    
    # Stores the results for this pixel 
    recovery_time <- data.frame(PixelID = pixels[[ii]], 
                                Recovery_time = repeated_ones_pos[1])
    recovery_time_table <- rbind(recovery_time_table, recovery_time)
    
  }
  
  return(recovery_time_table)
}

# This function calculates the cummulative percentages per a group within 
# A table, useful for the figures section of this code 
cum_percent_group <- function(df, group_var, x) {
  df <- df %>%
    group_by({{group_var}}) %>%
    mutate(cum_perc = cumsum({{x}})/sum({{x}})) %>%
    ungroup()
  
  return(df)
}
##### Calculate recovery time ##################################################
sd_limit <- 1 # Limit to determine what is considered within normal variance 
n_repeats <- 2 # Number of required years to be considered back within the normal range


# Run and time tecovery time function
start <- Sys.time()
results <- Calculate_recovery_time(baseline_data, 
                                   post_hw_data, 
                                   sd_limit, 
                                   n_repeats)
end <- Sys.time()

start - end

#results <- read.csv("Processed_data/data_tables/recovery_time.csv")
results # Final results

##### Export results ###########################################################
export_file_name <- paste0("Processed_data/data_tables/recovery_time_", cutoff_percent*100, "percent.csv")
write.csv(results, export_file_name, row.names = F)

##### Explore results ##########################################################
results %>% 
  count(Recovery_time) %>% 
  mutate(perent = (n/nrow(results))*100)

n_not_recovered <- sum(is.na(results$Recovery_time)) # NA's mean that they have not recovered yet

# Questions: 
## Do these percentages change with protection?

pixel_properties <- kelp_data_all %>% 
  filter(year == 2017) %>% 
  dplyr::select(lat, long, region, PixelID, mpa_status) %>% 
  distinct() %>% 
  full_join(results, by = "PixelID")

# Investigating protection 
n_mpa_status <- pixel_properties %>% 
  count(mpa_status) %>% 
  rename(n_total = n)

df_protection <- pixel_properties %>% 
  group_by(mpa_status) %>% 
  count(Recovery_time) %>% 
  rename(n_pixels = n) %>% 
  left_join(n_mpa_status, by = "mpa_status") %>% 
  mutate(percentage = n_pixels/n_total) %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("Full", "Partial", "None")))

df_protection <- cum_percent_group(df_protection, mpa_status, percentage)


plot1 <- ggplot(data = df_protection, aes(x = Recovery_time, 
                                          y = cum_perc, color = mpa_status)) +
  geom_point() +
  geom_line() +
  #ylim(0, .15) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  scale_color_manual(values = viridis(3)) +
  labs(x = "Recovery Time (years)", y = "Cummulative Percent of Pixels Recovered")

plot1

exportfile <- paste0("Figures/Recovery_time_", cutoff_percent*100, "percent.png")
png(exportfile, 
    units = "in", 
    height = 4, 
    width = 6, 
    res = 300)
plot1
dev.off()