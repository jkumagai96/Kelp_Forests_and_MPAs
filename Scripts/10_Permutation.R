# Date: January 11th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Bootstrap - randomization approach 
# BIO 202: Ecological Statistics

###### Research Question #######################################################
# Do marine protected areas increase the resilience of kelp forests (increase in kelp area)
# following a marine heat wave?

# Is there more kelp area in marine protected areas? Partial vs. fully protected?

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(sf)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")
points_in_mpas <- read.csv("Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv")

##### Structure Data ###########################################################
kelp_data <- kelp_data_all %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("None", "Partial", "Full"))) %>% 
  select(PixelID, year, area, mpa_status) %>% 
  drop_na()

##### Calculate True Values ####################################################
# Change global option to not print out message about group summaries 
options(dplyr.summarise.inform = FALSE)

# Calculate average kelp area per category per year 
true_values <- kelp_data %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_area = mean(area)) %>% 
  pivot_wider(names_from = mpa_status,
              values_from = avg_area) %>% 
  mutate(P_N = Partial - None,
         F_N = Full - None,
         F_P = Full - Partial) %>% 
  select(-c(None, Partial, Full))


##### Explore bias #############################################################
hist(kelp_data_all$lat)

p_mpas <- kelp_data_all %>% filter(mpa_status != "None") 
hist(p_mpas$lat)

##### Bootstrapping ############################################################
## Set up:
set.seed(20) # So the results are repeatable

# Set up variables outside of the loops
bootstrap_list <- list()
years <- 1984:2021

# Remove original mpa_status from kelp_data
kelp_data_r <- kelp_data %>% select(-mpa_status)

all_pixels <- kelp_data$PixelID %>% unique()

# within the for loop 
start <- Sys.time()
for (j in 1:10000) {
  
  # Sample random pixels that will overlap w/ MPAs 
  r_pixels <- sample(all_pixels, size = nrow(points_in_mpas), replace = F) # number of pixels originally overlapping with MPAs 
  
  # Assign the random pixels to the points in the mpas instead of the original values
  points_in_mpas$PixelID <- r_pixels
  
  # Join randomized mpa data with kelp data 
  kelp_w_mpas_r <- left_join(kelp_data_r, points_in_mpas, by = "PixelID") 
  
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
    summarise(avg_area = mean(area)) %>% 
    pivot_wider(names_from = mpa_status,
                values_from = avg_area) %>% 
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

final_results
0.05/(38*3)
test <- as.data.frame(final_results <= 0.0004385965)
test$year <- final_results$year

##### Export ###################################################################
final_results

write.csv(final_results, 
          "Processed_data/data_tables/bootstrap_1984_2021.csv", 
          row.names = F)

##### Figures ##################################################################
final_results <- read.csv("Processed_data/data_tables/bootstrap_1984_2021.csv")

# Identify years of marine heat waves based on plot
mhw_years <- data.frame(year = c(1992, 1997, 1998, 2014, 2015),
                        mhw = 1)

# Format Data into long format
results_long <- final_results %>% 
  pivot_longer(names_to = "Comparison", values_to = "pvalues", F_N:F_P) %>% 
  left_join(mhw_years, by = "year") %>% # Join the data together 
  mutate(mhw = replace_na(mhw, 0))

results_long$pvalues[results_long$pvalues == 0] <- .000001

# Plot 
results_long %>% 
  #filter(year >= 2008) %>% 
  ggplot(aes(x = year, y = -log10(pvalues), group = Comparison)) +
  geom_point(aes(color = Comparison, shape = Comparison), size = 2) +
  geom_hline(yintercept = -log10(0.05/(114))) +
  geom_hline(yintercept = -log10(0.05/3), linetype = "dashed") +
  scale_color_manual(values=c('#FF5C00', '#999999','#000EDD')) +
  theme_bw() +
  theme(legend.position = "bottom") +
  #scale_x_continuous(breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020)) +
  scale_x_continuous(breaks = c(1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)) +
  annotate("rect", xmin = 2014, xmax = 2015, ymin = 0, ymax = 6,
           alpha = .2,fill = "red") + 
  annotate("rect", xmin = 1997, xmax = 1998, ymin = 0, ymax = 6,
           alpha = .2,fill = "red") +
  labs(y = "-log10(P values)", x = "Year") 
  
# Export 
ggsave(last_plot(), filename = "Figures/Bootstrap_1984_to_2021.png",
       dpi = 600,
       units = "in", 
       height = 3, 
       width = 5)


#ggsave(last_plot(), filename = "Figures/Bootstrap_2008_to_2021.png",
#       dpi = 600,
#       units = "in", 
#       height = 3, 
#       width = 5)