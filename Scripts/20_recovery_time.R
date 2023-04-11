# Date: April 6th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Try to calculate recovery time 
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(sf)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")

##### Calculate baseline area per pixel ########################################
baseline_data <- kelp_data_all %>% 
  dplyr::select(PixelID, year, area) %>% 
  filter(year < 2014) %>% # 30 year baseline starting in 1984 to 2013
  group_by(PixelID) %>% 
  summarize(mean_area = mean(area, na.rm = T),
            sd_area = sd(area, na.rm = T)) 

# Create a dataset with the post heat wave data (2016-2021)
post_hw_data <- kelp_data_all %>% 
  dplyr::select(PixelID, year, area) %>% 
  filter(year > 2015) 

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

results # Final results

##### Explore results ##########################################################
results %>% 
  count(Recovery_time) %>% 
  mutate(perent = (n/nrow(results))*100)

n_not_recovered <- sum(is.na(results$Recovery_time)) # NA's mean that they have not recovered yet
# 214 pixels

# Questions: 
## Do these percentages change with protection?
## Do these percentages change with region?

pixel_properties <- kelp_data_all %>% 
  filter(year == 2016) %>% 
  dplyr::select(lat, long, region, PixelID, mpa_status) %>% 
  distinct() %>% 
  full_join(results, by = "PixelID") 

# Investigating protection 
n_mpa_status <- pixel_properties %>% 
  count(mpa_status) %>% 
  rename(n_total = n)

pixel_properties %>% 
  group_by(mpa_status) %>% 
  count(Recovery_time) %>% 
  rename(n_pixels = n) %>% 
  left_join(n_mpa_status, by = "mpa_status") %>% 
  mutate(percentage = 100*(n_pixels/n_total))

# Investigating region
n_region <- pixel_properties %>% 
  count(region) %>% 
  rename(n_total = n)

pixel_properties %>% 
  group_by(region) %>% 
  count(Recovery_time) %>% 
  rename(n_pixels = n) %>% 
  left_join(n_region, by = "region") %>% 
  mutate(percentage = 100*(n_pixels/n_total)) %>% 
  view()


# North coast recovered fastest??, followed by southern california, then north central coast, then central coast 

##### Export results ###########################################################
write.csv(results, "Processed_data/data_tables/recovery_time.csv", row.names = F)
