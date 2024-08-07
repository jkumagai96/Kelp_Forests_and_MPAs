# Date:July 24th 2024
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

# Declare Functions
std <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))

##### Offsets ##################################################################
data_all <- data_all %>% 
  mutate(log_transects = log(n_transects)) 

##### Process Data #############################################################
glm_data_south <- data_all %>% 
  filter(region == "South_Coast") %>% 
  filter(!is.na(urchin_total)) %>% 
  mutate(
    heatwave = factor(heatwave, levels = c("before", "during", "after")), 
    site_name = factor(site), 
    year_fct = factor(year),
    mpa_status = factor(mpa_status, levels = c("Reference", "Partial", "Full"))
  )

# Central california data 
glm_data_central <- data_all %>% 
  filter(region == "Central_Coast") %>% 
  filter(!is.na(urchin_total)) %>% 
  mutate(
    heatwave = factor(heatwave, levels = c("before", "during", "after")), 
    site_name = factor(site), 
    year_fct = factor(year),
    mpa_status = factor(mpa_status, levels = c("Reference", "Partial", "Full"))
  )%>% 
  mutate(urchins_d = round(urchins_d))

##### GLMM Models for interaction of MPAs and heatwave #########################
# Model 1 - Random intercepts
model_mpa1 <- glmmTMB(
  urchin_total ~ heatwave*mpa_status + offset(log_transects) +
    (1 | site_name),
  data = glm_data_south,
  family = nbinom1(link = "log")
)
# 
# model_mpa1 <- glmmTMB(
#   urchins_d ~ heatwave*mpa_status +
#     (1 | site_name),
#   data = glm_data_south,
#   family = tweedie(link = "log")
# )


simulationOutput_mpa_m1 <- simulateResiduals(model_mpa1, plot = F)
plot(simulationOutput_mpa_m1)

# Model 2 - Random intercept and slopes 
model_mpa2 <- glmmTMB(
  urchin_total ~ heatwave*mpa_status + offset(log_transects) +
    (1 + year | site_name), 
  data = glm_data_south, 
  family = nbinom1(link = "log")
) # For some reason only nbinom2 is working out of poisson and nbinom1 

# model_mpa2 <- glmmTMB(
#   urchins_d ~ heatwave*mpa_status +
#     (1 + year | site_name),
#   data = glm_data_south,
#   family = tweedie(link = "log")
# )



simulationOutput_mpa_m2 <- simulateResiduals(model_mpa2, plot = F)
plot(simulationOutput_mpa_m2)

# Random intercepts and autocorrelation structure 
model_mpa3 <- glmmTMB(
  urchin_total ~ heatwave*mpa_status + offset(log(n_transects)) +
    (1 + year | site_name) + 
    ar1(0 + year_fct | site_name), 
  data = glm_data_south, 
  family = poisson(link = "log")
) # Fit is not good with nbinom1 or nbinom2 or poisson

simulationOutput_mpa_m3 <- simulateResiduals(model_mpa3, plot = F)
plot(simulationOutput_mpa_m3)

# Model 1 residuals looks best but the results look the same for the two 
car::Anova(model_mpa1)
car::Anova(model_mpa2)
summary(model_mpa1)
summary(model_mpa2)

AIC(model_mpa1, model_mpa2)
# Model 2 is chosen, AIC is better and then it is consistent with central california 

emmeans::emmeans(model_mpa1, pairwise ~ heatwave, type = "response", offset = log(1/60))

z <- emmeans::emmeans(model_mpa1, pairwise ~ mpa_status | heatwave, type = "response", offset = log(1/60))
z

write.csv(broom::tidy(z$emmeans),
          "Processed_data/data_tables/emmeans_model_mpa2_new.csv",
          row.names = F)
write.csv(broom::tidy(z$contrasts),
          "Processed_data/data_tables/contrasts_model_mpa2_new.csv",
          row.names = F)

### Check autocorrelation
## Obtain residuals, pivot to wide data frame
resid_ts <- glm_data_south |>
  mutate(resid = resid(model_mpa2)) |>
  dplyr::select(site, year, resid) |>
  pivot_wider(names_from = "site", values_from = "resid")

## Compute normal ACF (partial ACF was giving some weird results)
## Not a problem for AR1 because PACF = ACF at lag 1
acf_df <- bind_rows(apply(
  ts(resid_ts[,-1]), MARGIN = 2,
  \(x) data.frame(lag = 0:5, acf = as.numeric(acf(x, na.action = na.pass, lag.max = 5, plot = FALSE)$acf))),
  .id = "site"
)
acf_df |> group_by(lag) |> summarize(acf = mean(acf, na.rm = TRUE))

##### Central california model
# Random Intercepts 
model_mpa1_central <- glmmTMB(
  urchin_total ~ heatwave*mpa_status + offset(log_transects) +
    (1 | site_name),
  data = glm_data_central,
  family = poisson(link = "log")
) # Poisson looks nice

# model_mpa1_central <- glmmTMB(
#   urchins_d ~ heatwave*mpa_status +
#     (1 | site_name),
#   data = glm_data_central,
#   family = tweedie(link = "log")
# )

simulationOutput_mpa_m1_c <- simulateResiduals(model_mpa1_central, plot = F)
plot(simulationOutput_mpa_m1_c)  

# Random intercepts and slopes 
model_mpa2_central <- glmmTMB(
  urchin_total ~ heatwave*mpa_status + offset(log_transects) +
    (1 + year | site_name),
  data = glm_data_central,
  family = nbinom2(link = "log")
) #nbinom2 is the best out of poisson and nbinom1


# model_mpa2_central <- glmmTMB(
#   urchins_d ~ heatwave*mpa_status +
#     (1 + year | site_name),
#   data = glm_data_central,
#   family = tweedie(link = "log")
# ) #


simulationOutput_mpa_m2_c <- simulateResiduals(model_mpa2_central, plot = F)
plot(simulationOutput_mpa_m2_c) 

# Random intercepts and autoccorelation
model_mpa3_central <- glmmTMB(
  urchin_total ~ heatwave*mpa_status + offset(log(n_transects)) +
    (1 | site_name) +
    ar1(0 + year_fct | site_name),
  data = glm_data_central, 
  family  = poisson(link = "log") # nbinom1 had covergence issues
)

simulationOutput_mpa_m3_c <- simulateResiduals(model_mpa3_central, plot = F)
plot(simulationOutput_mpa_m3_c)

AIC(model_mpa1_central, model_mpa2_central)

car::Anova(model_mpa1_central)
summary(model_mpa1_central)

car::Anova(model_mpa2_central)
summary(model_mpa2_central)

# Drastically different results between model_mpa1_central and model_mpa2_central

# Model 2 chosen as it has better AIC and the fit is ok and the results make sense 
emmeans::emmeans(model_mpa2_central, pairwise ~ heatwave, type = "response", offset = log(1/60))

z2 <- emmeans::emmeans(model_mpa2_central, pairwise ~ mpa_status | heatwave, type = "response", offset = log(1/60))
z2

write.csv(broom::tidy(z2$emmeans), 
           "Processed_data/data_tables/emmeans_model_mpa2_central_new.csv", row.names = F)
write.csv(broom::tidy(z2$contrasts), 
           "Processed_data/data_tables/contrasts_model_mpa2_central_new.csv", 
           row.names = F)

##### Model effects plots #####################################################
is_data <- as.data.frame(emmeans::emmeans(model_mpa2, "heatwave", by = "mpa_status", type = "response", offset = log(1/60))) %>%
  mutate(region = "Southern")
ic_data <- as.data.frame(emmeans::emmeans(model_mpa2_central, "heatwave", by = "mpa_status", type = "response", offset = log(1/60))) %>%
  mutate(region = "Central") 
interaction_data <- rbind(is_data, ic_data)

group.colors <- c(Full = "#440154", Reference = "#FFBA00", Partial ="#21918c")

interaction_plot <- interaction_data %>%
  ggplot(aes(heatwave, color = mpa_status)) +
  geom_pointrange(aes(y = response, ymin = asymp.LCL, ymax = asymp.UCL),
                  position = position_dodge(width = 0.2),
                  linewidth = 1) +
  facet_grid(cols = vars(region)) +
  scale_color_manual(values=group.colors, name = "Protection Status") +
  ylab(bquote('Model Response - Urchins per ' ~m^2)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  xlab("Heatwave period")
interaction_plot

Urchins_per_region <- data_all %>%
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern")) %>%
  group_by(year, mpa_status, region) %>%
  summarise(total = mean(urchins_d, na.rm = T)/60,
            se_total = std(urchins_d)/60) %>%
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
