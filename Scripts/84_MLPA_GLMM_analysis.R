# Date: Nov. 6th 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Run linear regression statistics on MLPA data without site blocks
# BIO 202: Ecological Statistics

##### Load Packages and Data ###################################################
library("MASS")
library("ggeffects")
library("nlme")
library("emmeans")
library("DHARMa")
library("glmmTMB")
library("cowplot")
library("patchwork")
library("png")
library("tidyverse")

# Add in little black outlines of organisms
sheephead_3 <- readPNG("Figures/Organisms/Sheephead_3.png", native = TRUE)
spiny_lobster <- readPNG("Figures/Organisms/Spiny_Lobster_Figure_t.png", native = TRUE)

# Load Data
data_all <- read.csv("Processed_data/MLPA_data_summarized_wo_siteblocks.csv") 

# Empirical standard error function
se <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))

##### Process Data #############################################################
glm_data <- data_all |> 
  filter(!is.na(urchin_total)) |> 
  mutate(
    heatwave = factor(heatwave, levels = c("before", "during", "after")), 
    site_name = factor(site), 
    year_fct = factor(year),
    year_std = (year - 2013)/sd(year), 
    mpa_status = factor(mpa_status, levels = c("Reference", "Partial", "Full"))
  )

glm_data_south <- glm_data |> filter(region == "South_Coast")

glm_data_central <- glm_data |> filter(region == "Central_Coast")

##### Southern Cali. GLMM Models for interaction of MPAs and heatwave ##########
# Model 1 - Random intercepts
model_mpa1 <- glmmTMB(
  urchin_total ~ heatwave*mpa_status + offset(log(n_transects)) +
    (1 | site_name),
  data = glm_data_south,
  family = nbinom1(link = "log")
)

simulationOutput_mpa_m1 <- simulateResiduals(model_mpa1, plot = F)
plot(simulationOutput_mpa_m1)

# 
# model_mpa1 <- glmmTMB(
#   urchins_d ~ heatwave*mpa_status +
#     (1 | site_name),
#   data = glm_data_south,
#   family = tweedie(link = "log")
# )


# Model 2 - Random intercept and slopes 
model_mpa2 <- glmmTMB(
  urchin_total ~ heatwave*mpa_status + offset(log(n_transects)) +
    (1 + year_std | site_name), 
  data = glm_data_south, 
  family = nbinom1(link = "log")
) 

simulationOutput_mpa_m2 <- simulateResiduals(model_mpa2, plot = F)
plot(simulationOutput_mpa_m2)


# Random intercepts and autocorrelation structure 
# Fit is not good with nbinom1 or nbinom2 or poisson
model_mpa3 <- glmmTMB(
  urchin_total ~ heatwave*mpa_status + offset(log(n_transects)) +
    (1 + year_std | site_name) + 
    ar1(0 + year_fct | site_name), 
  data = glm_data_south, 
  family = poisson(link = "log")
)

simulationOutput_mpa_m3 <- simulateResiduals(model_mpa3, plot = F)
plot(simulationOutput_mpa_m3)


## Model comparison -----------------------------------------------------------

# Model 1 residuals looks best but the results look the same for the two 
car::Anova(model_mpa1)
car::Anova(model_mpa2)
summary(model_mpa1)
summary(model_mpa2)

AIC(model_mpa1, model_mpa2, model_mpa3)

# Model 2 is chosen, AIC is better and then it is consistent with central california 

## Check autocorrelation -------------------------------------------------------

# Obtain residuals, pivot to wide data frame
resid_ts <- glm_data_south |>
  mutate(resid = resid(model_mpa2)) |>
  dplyr::select(site, year, resid) |>
  pivot_wider(names_from = "site", values_from = "resid")

# Compute normal ACF (partial ACF was giving some weird results)
# Not a problem for AR1 because PACF = ACF at lag 1
acf_df <- bind_rows(apply(
  ts(resid_ts[,-1]), MARGIN = 2,
  \(x) data.frame(lag = 0:5, acf = as.numeric(acf(x, na.action = na.pass, lag.max = 5, plot = FALSE)$acf))),
  .id = "site"
)
acf_df |> group_by(lag) |> summarize(acf = mean(acf, na.rm = TRUE))

# Central California GLMMs ----------------------------------------------------
# Random Intercepts - chose Poisson
model_mpa1_central <- glmmTMB(
  urchin_total ~ heatwave*mpa_status + offset(log(n_transects)) +
    (1 | site_name),
  data = glm_data_central,
  family = poisson(link = "log")
)

simulationOutput_mpa_m1_c <- simulateResiduals(model_mpa1_central, plot = F)
plot(simulationOutput_mpa_m1_c)  

# Random intercepts and slopes - chose nbinom2
model_mpa2_central <- glmmTMB(
  urchin_total ~ heatwave*mpa_status + offset(log(n_transects)) +
    (1 + year_std | site_name),
  data = glm_data_central,
  family = nbinom2(link = "log")
)

simulationOutput_mpa_m2_c <- simulateResiduals(model_mpa2_central, plot = F)
plot(simulationOutput_mpa_m2_c) 

# Random intercepts and autoccorelation - nbinom1 had covergence issues
model_mpa3_central <- glmmTMB(
  urchin_total ~ heatwave*mpa_status + offset(log(n_transects)) +
    (1 + year_std | site_name) +
    ar1(0 + year_fct | site_name),
  data = glm_data_central, 
  family  = poisson(link = "log") # 
)

simulationOutput_mpa_m3_c <- simulateResiduals(model_mpa3_central, plot = F)
plot(simulationOutput_mpa_m3_c)

## Model Comparison ------------------------------------------------------------
AIC(model_mpa1_central, model_mpa2_central)

car::Anova(model_mpa1_central)
summary(model_mpa1_central)

car::Anova(model_mpa2_central)
summary(model_mpa2_central)

# Drastically different results between model_mpa1_central and model_mpa2_central

# Model 2 chosen as it has better AIC and the fit is ok and the results make sense 

##### Bias correctiosn for plots and estimated marginal means ------------------

# Mean standardized years corresponding to heatwave times
hw_yrs <- aggregate(glm_data$year_std, by = list(glm_data$heatwave), mean)[,2]

# Extract covariance matrices of random effects
sigma_re_s <- glmmTMB::VarCorr(model_mpa2)$cond$site_name
sigma_re_c <- glmmTMB::VarCorr(model_mpa2_central)$cond$site_name

# Compute total random effects standard deviation at each time point
sigma_re_s <- sapply(hw_yrs, \(t) sqrt(sum(c(1, t) * sigma_re_s %*% c(1, t))))
sigma_re_c <- sapply(hw_yrs, \(t) sqrt(sum(c(1, t) * sigma_re_c %*% c(1, t))))

# Compute bias-corrected estimated marginal means for SoCal model
ems <- ref_grid(model_mpa2, at = list(heatwave = levels(glm_data$heatwave), mpa_status = levels(glm_data$mpa_status)))
ems <- as.data.frame(emmip(
  ems, ~ heatwave + mpa_status, type = "response", CIs = TRUE, plot = FALSE, 
  bias.adjust = TRUE, sigma = rep(sigma_re_s, 3), offset = log(1/60)
))

# Compute bias-corrected estimated marginal means for CenCal model
emc <- ref_grid(model_mpa2_central, at = list(heatwave = levels(glm_data$heatwave), mpa_status = levels(glm_data$mpa_status)))
emc <- as.data.frame(emmip(
  emc, ~ heatwave + mpa_status, type = "response", CIs = TRUE, plot = FALSE, 
  bias.adjust = TRUE, sigma = rep(sigma_re_c, 3), offset = log(1/60)
))

# Join
em <- bind_rows(list(Southern = ems, Central = emc), .id = "region")

##### Export emmeans for both sets of models -----------------------------------
# Southern California 
mpa_means <- emmeans::emmeans(model_mpa2, pairwise ~ mpa_status | heatwave, type = "response", offset = log(1/60),
                              bias.adj = TRUE, sigma = rep(sigma_re_s, each = 3))
mpa_means

write.csv(broom::tidy(mpa_means$emmeans),
          "Processed_data/data_tables/emmeans_model_mpa2_new.csv",
          row.names = F)
write.csv(broom::tidy(mpa_means$contrasts),
          "Processed_data/data_tables/contrasts_model_mpa2_new.csv",
          row.names = F)

# Central California 
emmeans::emmeans(model_mpa2_central, pairwise ~ heatwave, type = "response", offset = log(1/60),
                 bias.adj = TRUE, sigma = sigma_re_c)

mpa_means2 <- emmeans::emmeans(model_mpa2_central, pairwise ~ mpa_status | heatwave, type = "response", offset = log(1/60),
                               bias.adj = TRUE, sigma = rep(sigma_re_c, each = 3))
mpa_means2

write.csv(broom::tidy(mpa_means2$emmeans), 
           "Processed_data/data_tables/emmeans_model_mpa2_central_new.csv", row.names = F)
write.csv(broom::tidy(mpa_means2$contrasts), 
           "Processed_data/data_tables/contrasts_model_mpa2_central_new.csv", 
           row.names = F)

# Plot estimated marginal means
group.colors <- c(Full = "#440154", Reference = "#FFBA00", Partial ="#21918c")

interaction_plot <- em %>%
  ggplot(aes(heatwave, color = mpa_status)) +
  geom_pointrange(aes(y = yvar, ymin = LCL, ymax = UCL),
                  position = position_dodge(width = 0.2),
                  linewidth = 1) +
  facet_grid(cols = vars(region)) +
  scale_color_manual(values=group.colors, name = "Protection Status") +
  ylab(bquote('Model Response - Urchins per ' ~m^2)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank()) +
  xlab("Heatwave period")
interaction_plot

# Plot empirical trends over time
Urchins_per_region <- data_all %>%
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern")) %>%
  group_by(year, mpa_status, region) %>%
  summarise(total = mean(urchins_d, na.rm = T)/60,
            se_total = se(urchins_d)/60) %>%
  ggplot(aes(x = year, y = total, color = mpa_status)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = total - se_total,
                    ymax = total + se_total),
                width = 0.2,
                linewidth = .5,
                alpha = 0.4) +
  facet_grid(cols = vars(region)) +
  scale_color_manual(values=group.colors, name = "Protection Status",
                     labels = c("Full", "Partial", "Unprotected")) +
  ylab(bquote('Urchins per ' ~m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2,
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = c(.1, .75))

urchins_and_interactions_plot <- plot_grid(Urchins_per_region, interaction_plot,
                                           labels = "AUTO",
                                           nrow = 2)
urchins_and_interactions_plot

png(filename = "Figures/GLMM_heatwave_protection_plot_PISCO_new.png", 
    width = 8, 
    height = 8,
    units = "in", 
    res = 600)
urchins_and_interactions_plot
dev.off()

##### Export Residual Plots ####################################################
# Southern California
simulationOutput_mpa_m2 <- simulateResiduals(model_mpa2, plot = F)

png(filename = "Figures/Residuals_mpa2_south.png", 
    width = 8, 
    height = 5,
    units = "in", 
    res = 600)
plot(simulationOutput_mpa_m2)
dev.off()

# Central California
simulationOutput_mpa_m2_c <- simulateResiduals(model_mpa2_central, plot = F)

png(filename = "Figures/Residuals_mpa2_central.png", 
    width = 8, 
    height = 5,
    units = "in", 
    res = 600)
plot(simulationOutput_mpa_m2_c)
dev.off()
