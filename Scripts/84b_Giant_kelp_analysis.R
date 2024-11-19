# Date: Nov. 6th 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Conduct statistical analysis if giant kelp abundances are higher during and 
#          after the heatwave in protected areas
# BIO 202: Ecological Statistics

##### Set up ###################################################################
# Load packages
library(tidyverse)
library(emmeans)
library(glmmTMB)
library(DHARMa)

# Load Data
data_all <- read.csv("Processed_data/MLPA_data_summarized_wo_siteblocks.csv") %>% 
  mutate(mpa_status = ifelse(mpa_status == "Reference", "Unprotected", mpa_status))

# Declare Functions
std <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))

##### Format Data ##############################################################
data_all <- data_all %>% 
  filter(!is.na(MACPYRAD_d)) %>% 
  mutate(
    heatwave = factor(heatwave, levels = c("before", "during", "after")), 
    site_name = factor(site), 
    year_fct = factor(year),
    year_std = (year - 2013)/sd(year),
    mpa_status = factor(mpa_status, levels = c("Unprotected", "Partial", "Full"))
  )

central_df <- data_all %>% 
  filter(region == "Central_Coast")

south_df <- data_all %>% 
  filter(region == "South_Coast")

##### Plot data ################################################################

data_all %>% 
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern")) %>% 
  group_by(year, region) %>% 
  summarise(total = mean(kelp_d)/60,
            total_se = std(kelp_d)/60) %>% 
  ggplot(aes(x = year, y = total)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = total - total_se, 
                    ymax = total + total_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  facet_grid(cols = vars(region)) +
  ylab("Number of Kelp") +
  ylab(bquote('Kelp per ' ~ m^2)) +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))

##### Statistical Modeling #####################################################
### Southern California
# Model 1 - random intercepts
m1_kelp_s <- glmmTMB(
  kelp_d ~ heatwave*mpa_status +
    (1 | site_name),
  data = south_df,
  family = tweedie(link = "log")
)
m1_kelp_s

sim_output1 <- simulateResiduals(m1_kelp_s, plot = F)
plot(sim_output1) 

# Model 2 - random slopes and intercepts
m2_kelp_s <- glmmTMB(
  kelp_d ~ heatwave*mpa_status +
    (1 + year_std | site_name),
  data = south_df,
  family = tweedie(link = "log")
)
m2_kelp_s

sim_output2 <- simulateResiduals(m2_kelp_s, plot = F)
plot(sim_output2)  


# Model 3 - random slopes and intercepts and ar1
m3_kelp_s <- glmmTMB(
  kelp_d ~ heatwave*mpa_status +
    (1 + year_std | site_name) +
    ar1(0 + year_fct | site_name),
  data = south_df,
  family = tweedie(link = "log")
)
m3_kelp_s

sim_output3 <- simulateResiduals(m3_kelp_s, plot = F)
plot(sim_output3) # Fit is much worse, removed from consideration

# Model 4 - random intercepts and ar1
m4_kelp_s <- glmmTMB(
  kelp_d ~ heatwave*mpa_status +
    (1 | site_name) +
    ar1(0 + year_fct | site_name),
  data = south_df,
  family = tweedie(link = "log")
)
m4_kelp_s # seems that the random intercepts are now not doing much that the AR1 is (low std.dev. value)
sim_output4 <- simulateResiduals(m4_kelp_s, plot = F)
plot(sim_output4) # Fit is worse than models 2 and 1 

### Model Comparison for Southern California kelp
AIC(m1_kelp_s, m2_kelp_s)

# When we remove models 3 and 4 as an option due to bad fit, 
# Model 2 is chosen as no convergence issues, the residuals look the same and AIC is lower 
car::Anova(m1_kelp_s)
car::Anova(m2_kelp_s)

summary(m2_kelp_s)


### Central California 
# Random interecepts
m1_kelp_c <- glmmTMB(
  kelp_d ~ heatwave*mpa_status +
    (1 | site_name),
  data = central_df,
  family = tweedie(link = "log")
)
m1_kelp_c

sim_output4 <- simulateResiduals(m1_kelp_c, plot = F)
plot(sim_output4)  

# Random slopes and intercepts
m2_kelp_c <- glmmTMB(
  kelp_d ~ heatwave*mpa_status +
    (1 + year_std | site_name),
  data = central_df,
  family = tweedie(link = "log")
)
m2_kelp_c

sim_output5 <- simulateResiduals(m2_kelp_c, plot = F)
plot(sim_output5)  


# Random slopes and intercepts and ar1
m3_kelp_c <- glmmTMB(
  kelp_d ~ heatwave*mpa_status +
    (1 + year_std | site_name) +
    ar1(0 + year_fct | site_name),
  data = central_df,
  family = tweedie(link = "log")
)
m3_kelp_c # No longer Convergence issues 

sim_output6 <- simulateResiduals(m3_kelp_c, plot = F)
plot(sim_output6) # Better residuals actually 

# Random intercepts and ar1
m4_kelp_c <- glmmTMB(
  kelp_d ~ heatwave*mpa_status +
    (1 + year_std | site_name),
  data = central_df,
  family = tweedie(link = "log")
)
m4_kelp_c # No longer Convergence issues 

sim_output7 <- simulateResiduals(m4_kelp_c, plot = F)
plot(sim_output7) # slighlty worse residuals

### Model comparison for central california kelp
AIC(m1_kelp_c, m2_kelp_c, m3_kelp_c, m4_kelp_c) # Model 3, m3_kelp_c is best 

# Best residuals, lowest AIC, 

car::Anova(m1_kelp_c)
car::Anova(m2_kelp_c)
car::Anova(m3_kelp_c)
car::Anova(m4_kelp_c)

summary(m3_kelp_c)

##### Bias correction ##########################################################
# Mean standardized years corresponding to heatwave times
hw_yrs <- aggregate(data_all$year_std, by = list(data_all$heatwave), mean)[,2]

# Extract covariance matrices of random effects
sigma_re_s <- glmmTMB::VarCorr(m2_kelp_s)$cond$site_name
sigma_re_c <- glmmTMB::VarCorr(m3_kelp_c)$cond$site_name

# Compute total random effects standard deviation at each time point
sigma_re_s <- sapply(hw_yrs, \(t) sqrt(sum(c(1, t) * sigma_re_s %*% c(1, t))))
sigma_re_c <- sapply(hw_yrs, \(t) sqrt(sum(c(1, t) * sigma_re_c %*% c(1, t))))

# Compute bias-corrected estimated marginal means for SoCal model
ems <- ref_grid(m2_kelp_s, at = list(heatwave = levels(data_all$heatwave), mpa_status = levels(data_all$mpa_status)))
ems <- as.data.frame(emmip(
  ems, ~ heatwave + mpa_status, type = "response", CIs = TRUE, plot = FALSE, 
  bias.adjust = TRUE, sigma = rep(sigma_re_s, 3)
))

# Compute bias-corrected estimated marginal means for CenCal model
emc <- ref_grid(m3_kelp_c, at = list(heatwave = levels(data_all$heatwave), mpa_status = levels(data_all$mpa_status)))
emc <- as.data.frame(emmip(
  emc, ~ heatwave + mpa_status, type = "response", CIs = TRUE, plot = FALSE, 
  bias.adjust = TRUE, sigma = rep(sigma_re_c, 3)
))

# Join
em <- bind_rows(list(Southern = ems, Central = emc), .id = "region")

### Export of estimated marginal means 
kelpmeans_s <- emmeans::emmeans(m2_kelp_s, pairwise ~ mpa_status | heatwave, type = "response",
                      bias.adj = TRUE, sigma = rep(sigma_re_s, each = 3))
kelpmeans_s

write.csv(broom::tidy(kelpmeans_s$emmeans),
          "Processed_data/data_tables/emmeans_kelp2_south.csv",
          row.names = F)
write.csv(broom::tidy(kelpmeans_s$contrasts),
          "Processed_data/data_tables/contrasts_kelp2_south.csv",
          row.names = F)


# There is less kelp after the heatwave in Central California

kelpmeans_c <- emmeans::emmeans(m3_kelp_c, pairwise ~ mpa_status | heatwave, type = "response",
                                bias.adj = TRUE, sigma = rep(sigma_re_c, each = 3))
kelpmeans_c

emmeans::emmeans(m3_kelp_c, pairwise ~ heatwave, type = "response",
                bias.adj = TRUE, sigma = sigma_re_c)

write.csv(broom::tidy(kelpmeans_c$emmeans), 
          "Processed_data/data_tables/emmeans_kelp2_central_new.csv", row.names = F)
write.csv(broom::tidy(kelpmeans_c$contrasts), 
          "Processed_data/data_tables/contrasts_kelp2_central_new.csv", 
          row.names = F)

##### Figures ##################################################################
group.colors <- c(Full = "#440154", Unprotected = "#FFBA00", Partial ="#21918c")

interaction_plot <- em %>%
  ggplot(aes(heatwave, color = mpa_status)) +
  geom_pointrange(aes(y = yvar/60, ymin = LCL/60, ymax = UCL/60),
                  position = position_dodge(width = 0.2),
                  linewidth = 1) +
  facet_grid(cols = vars(region)) +
  scale_color_manual(values=group.colors, name = "Protection Status") +
  ylab(bquote('Model Response - Adult Kelps per ' ~m^2)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  xlab("Heatwave period")

interaction_plot

kelp_per_region <- data_all %>% 
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern")) %>% 
  group_by(year, mpa_status, region) %>% 
  summarise(total = mean(kelp_d)/60,
            total_se = std(kelp_d)/60) %>% 
  ggplot(aes(x = year, y = total, color = mpa_status)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = total - total_se, 
                    ymax = total + total_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  facet_grid(cols = vars(region)) +
  scale_color_manual(values=group.colors, name = "Protection Status") +
  ylab("Number of Kelp") +
  ylab(bquote('Adult Kelps per ' ~ m^2)) +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = c(.15, .8),
        legend.background = element_rect(linewidth = 0.2, colour = 1),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))

kelp_per_region

##### Export ###################################################################
kelp_and_interactions_plot <- cowplot::plot_grid(kelp_per_region, interaction_plot,
                                           labels = "AUTO",
                                           nrow = 2)
kelp_and_interactions_plot

png(filename = "Figures/Kelp_interactions_plots.png", 
    width = 8, 
    height = 8,
    units = "in", 
    res = 600)
kelp_and_interactions_plot
dev.off()

png(filename = "Figures/Residuals_kelp2_central.png", 
    width = 8, 
    height = 5,
    units = "in", 
    res = 600)
plot(sim_output6) 
dev.off()

png(filename = "Figures/Residuals_kelp2_south.png", 
    width = 8, 
    height = 5,
    units = "in", 
    res = 600)
plot(sim_output2) 
dev.off()

##### Exploration ##############################################################
# Figure out if there are more sites with kelp or higher densities...
data_all %>% 
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern")) %>% 
  mutate(yes_kelp = kelp_d > 0,
         no_kelp = kelp_d == 0) %>% 
  group_by(year, region) %>% 
  summarize(n_yes = sum(yes_kelp),
            n_no = sum(no_kelp)) %>% 
  pivot_longer(n_yes:n_no, values_to = "n", names_to = "group") %>% 
  ggplot(aes(x = year, y = n, fill = group)) +
  geom_bar(stat = "identity") +
  facet_grid(cols = vars(region)) +
  theme_bw()

plot_sites_w_kelp <- data_all %>% 
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern")) %>% 
  mutate(yes_kelp = kelp_d > 0) %>% 
  group_by(year, region) %>% 
  summarize(p_yes = sum(yes_kelp)/n() ) %>% 
  ggplot(aes(x = year, y = p_yes, col = region)) +
  scale_color_manual(values = c("#DC3220", "#005AB5")) +
  geom_line(linewidth = 1.5) +
  theme_bw() +
  labs(y = "Proporation of sites with kelp", x = "Year")

plot_sites_w_kelp

png(filename = "Figures/Sites_w_kelp_over_time.png", 
    width = 6, 
    height = 4,
    units = "in", 
    res = 600)
plot_sites_w_kelp
dev.off()


