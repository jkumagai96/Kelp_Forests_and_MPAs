# Date: April 20th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Bootstrap - Permutation approach based on difference in percent recovery instead of kelp area 
#          But testing the other side of the distribution
# BIO 202: Ecological Statistics

##### Set Up: Packages and Data ################################################
# Packages
library(tidyverse)
library(doParallel)
library(foreach)

# Data 
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv") 
points_in_mpas <- read.csv("Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv")

# Would only be looking at data from Central and Southern Coast of California
# Would only be looking at data from 2016-2021 due to calculating percent recovery 

##### Declare Functions ########################################################
Calculate_percent_recovery <- function(x, min_v, max_v) {
  value <- 100*(x - min_v)/(max_v - min_v)
  return(value)
}

##### Format Data ##############################################################
data_area <- kelp_data_all %>% 
  filter(region == "South_Coast" | region == "Central_Coast") %>% 
  filter(year < 2014) %>% 
  select(PixelID, year, area) %>% 
  pivot_wider(names_from = year, values_from = area) 

# Calculate Max and Min value 
data_area$min <- apply(data_area, 1, FUN = min, na.rm = T)
data_area$max <- apply(data_area, 1, FUN = max, na.rm = T)

df <- data_area %>% 
  select(PixelID, min, max)

# Calculate Percent Recovery for each row in the dataset from 2016 to 2021
df_percent_recovery <- kelp_data_all %>% 
  filter(region == "South_Coast" | region == "Central_Coast") %>% 
  filter(year >= 2016) %>% 
  select(PixelID, year, area) %>% 
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
  filter(region == "South_Coast" | region == "Central_Coast") %>% 
  filter(year >= 2016) %>% 
  select(PixelID, year, mpa_status) %>% 
  left_join(df_percent_recovery, by = c("PixelID", "year"))

##### Calculate True Values ####################################################
# Change global option to not print out message about group summaries 
options(dplyr.summarise.inform = FALSE)

# Calculate average/median percent kelp recovery per category per year 
# Change is here don't forget 
true_values <- kelp_data %>% 
  group_by(year, mpa_status) %>% 
  summarise(median_pr = median(percent_recovery)) %>% 
  pivot_wider(names_from = mpa_status,
              values_from = median_pr) %>% 
  mutate(N_P = None - Partial,
         N_F = None - Full,
         P_F = Partial - Full) %>% 
  select(-c(None, Partial, Full))

##### Bootstrapping ############################################################
## Set up:
set.seed(20) # So the results are repeatable

# Set up variables outside of the loops
# bootstrap_list <- list()
years <- 2016:2021

# Remove original mpa_status from kelp_data
kelp_data_r <- kelp_data %>% select(-mpa_status)

all_pixels <- kelp_data$PixelID %>% unique()

#create the cluster
n_cores <- 2

my.cluster <- parallel::makeCluster(n_cores)

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParWorkers()

# within the for loop 
start <- Sys.time()
bootstrap_list <- foreach (j = 1:10000) %dopar% {
  require(tidyverse)
  
  # Sample random pixels that will overlap w/ MPAs 
  r_pixels <- sample(all_pixels, size = nrow(points_in_mpas), replace = F) # number of pixels originally overlapping with MPAs 
  
  # Assign the random pixels to the points in the mpas instead of the original values
  points_in_mpas$PixelID <- r_pixels
  
  # Join randomized mpa data with kelp data 
  kelp_w_mpas_r <- dplyr::left_join(kelp_data_r, points_in_mpas, by = "PixelID") 
  
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
  
  # Calculate table, changes here too!
  values <- new_data %>% 
    group_by(year, mpa_status) %>% 
    summarise(median_pr = median(percent_recovery)) %>% 
    pivot_wider(names_from = mpa_status,
                values_from = median_pr) %>% 
    mutate(N_P = None - Partial,
           N_F = None - Full,
           P_F = Partial - Full) %>% 
    select(-c(None, Partial, Full)) %>% 
    ungroup()
  
  values
}
end <- Sys.time()

start - end
parallel::stopCluster(cl = my.cluster)

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
  
  sig_N_F <- 1 - ecdf(df$N_F)(tv$N_F)
  sig_N_P <- 1 - ecdf(df$N_P)(tv$N_P)
  sig_P_F <- 1 - ecdf(df$P_F)(tv$P_F)
  
  t <- data.frame(years[i],
                  sig_N_F,
                  sig_N_P,
                  sig_P_F)
  
  colnames(t) <- c("year", "N_F", "N_P", "P_F")
  
  final_results <- rbind(final_results, t)
  
}

test <- as.data.frame(final_results <= 0.0004385965)
test$year <- final_results$year

##### Export ###################################################################
final_results

write.csv(final_results, 
          "Processed_data/data_tables/bootstrap_percent_recovery_median_reversed.csv", 
          row.names = F)

print('bootstrap results saved')
##### Figures ##################################################################
# Identify years of marine heat waves based on plot
mhw_years <- data.frame(year = c(2014, 2015),
                        mhw = 1)

# Format Data into long format
results_long <- final_results %>% 
  pivot_longer(names_to = "Comparison", values_to = "pvalues", N_F:P_F) %>% 
  left_join(mhw_years, by = "year") %>% # Join the data together 
  mutate(mhw = replace_na(mhw, 0))

results_long$pvalues[results_long$pvalues == 0] <- .000001

print('now attempting to make figure')
# Plot 
plot1 <- results_long %>% 
  ggplot(aes(x = year, y = -log10(pvalues), group = Comparison)) +
  geom_point(aes(color = Comparison, shape = Comparison), size = 2) +
  geom_hline(yintercept = -log10(0.05/(18))) +
  geom_hline(yintercept = -log10(0.05/3), linetype = "dashed") +
  scale_color_manual(values=c('#FF5C00', '#999999','#000EDD')) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(2016, 2017, 2018, 2019, 2020, 2021)) +
  labs(y = "-log10(P values)", x = "Year") 

# Export
print('attempting export')

png("Figures/Bootstrap_percent_recovery_median_reversed.png", width = 5, height = 3, units = "in", res = 600)
plot1
dev.off() 
