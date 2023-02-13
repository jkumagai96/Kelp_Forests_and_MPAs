# Date: February 13th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Visualize bootstrap - create figures
# BIO 202: Ecological Statistics

##### Set Up: Packages and Data ################################################
library(tidyverse)

##### Plot Function #############################################################
mhw_years <- data.frame(year = c(1992, 1997, 1998, 2014, 2015),
                        mhw = 1)

Plot_Pvalues <- function(data_final) {
  
  results_long <- data_final %>% 
    pivot_longer(names_to = "Comparison", values_to = "pvalues", F_N:F_P) %>% 
    left_join(mhw_years, by = "year") %>% # Join the data together 
    mutate(mhw = replace_na(mhw, 0))
  
  results_long$pvalues[results_long$pvalues == 0] <- .000001
  
  # Calcualtes the number of comparisions to implement Bonferoniis correction within the plot 
  n_comparisons <- na.omit(results_long) %>% nrow()
  
  plot1 <- results_long %>% 
    ggplot(aes(x = year, y = -log10(pvalues), group = Comparison)) +
    geom_point(aes(color = Comparison, shape = Comparison), size = 2) +
    geom_hline(yintercept = -log10(0.05/(n_comparisons))) +
    geom_hline(yintercept = -log10(0.05/3), linetype = "dashed") +
    scale_color_manual(values=c('#FF5C00', '#999999','#000EDD')) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = c(1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)) +
    annotate("rect", xmin = 2014, xmax = 2015, ymin = 0, ymax = 6,
             alpha = .2,fill = "red") + 
    annotate("rect", xmin = 1997, xmax = 1998, ymin = 0, ymax = 6,
             alpha = .2,fill = "red") +
    labs(y = "-log10(P values)", x = "Year") 
  
  return(plot1)
}

Plot_Pvalues_reversed <- function(data_final) {
 
  results_long <- final_results %>% 
    pivot_longer(names_to = "Comparison", values_to = "pvalues", N_F:P_F) %>% 
    left_join(mhw_years, by = "year") %>% # Join the data together 
    mutate(mhw = replace_na(mhw, 0))
  
  results_long$pvalues[results_long$pvalues == 0] <- .000001
  
  # Calcualtes the number of comparisions to implement Bonferoniis correction within the plot 
  n_comparisons <- na.omit(results_long) %>% nrow()
  
  # Plot 
  plot2 <- results_long %>% 
    ggplot(aes(x = year, y = -log10(pvalues), group = Comparison)) +
    geom_point(aes(color = Comparison, shape = Comparison), size = 2) +
    geom_hline(yintercept = -log10(0.05/(n_comparisons))) +
    geom_hline(yintercept = -log10(0.05/3), linetype = "dashed") +
    scale_color_manual(values=c('#FF5C00','#000EDD', '#999999')) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = c(1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)) +
    annotate("rect", xmin = 2014, xmax = 2015, ymin = 0, ymax = 6,
             alpha = .2,fill = "red") + 
    annotate("rect", xmin = 1997, xmax = 1998, ymin = 0, ymax = 6,
             alpha = .2,fill = "red") +
    labs(y = "-log10(P values)", x = "Year") 
  
  return(plot2)
}

##### Create Figures and Export ################################################
region <- c("south", "central", "north_central", "north")

for (i in 1:length(region)) {
  file_name <- paste0("Processed_data/data_tables/bootstrap_", region[i] ,"_coast.csv")
  final_results <- read.csv(file_name)
  
  figure_file_name <- paste0("Figures/regional_bootstrap/Bootstrap_", region[i], "_coast.png")
  
  Plot_Pvalues(final_results)
  
  ggsave(last_plot(), filename = figure_file_name,
         dpi = 600,
         units = "in", 
         height = 3, 
         width = 5)
}

for (i in 1:length(region)) {
  file_name <- paste0("Processed_data/data_tables/bootstrap_", region[i] ,"_coast_reversed.csv")
  final_results <- read.csv(file_name)
  
  figure_file_name <- paste0("Figures/regional_bootstrap/Bootstrap_", region[i], "_coast_reversed.png")
  
  Plot_Pvalues_reversed(final_results)
  ggsave(last_plot(), filename = figure_file_name,
         dpi = 600,
         units = "in", 
         height = 3, 
         width = 5)
}

##### End Script ###############################################################