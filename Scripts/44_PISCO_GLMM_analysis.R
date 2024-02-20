# Date: February 8th 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Run linear regression statistics on PISCO data 
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
data_all <- read.csv("Processed_data/PISCO_data_summarized.csv") %>% 
  rename(total_urchins = urchin_d)

# Declare Functions
std <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))
##### Process Data #############################################################
# Filter to just southern california 
data_urchins <- data_all %>% 
  filter(region == "South_Coast")

glm_data <- data_urchins %>% 
  mutate(
    heatwave = factor(heatwave, levels = c("before", "during", "after")), 
    site_name = factor(site), 
    year_fct = factor(year),
    mpa_status = factor(mpa_status, levels = c("Reference", "Partial", "Full"))
  )

# Central california data 
glm_data3 <- data_all %>% 
  filter(region == "Central_Coast") %>% 
  mutate(
    heatwave = factor(heatwave, levels = c("before", "during", "after")), 
    site_name = factor(site), 
    year_fct = factor(year),
    mpa_status = factor(mpa_status, levels = c("Reference", "Partial", "Full"))
  )

##### GLMM Models for interaction of MPAs and heatwave #########################
# Model 1 - Random intercepts
model_mpa1 <- glmmTMB(
  total_urchins ~ heatwave*mpa_status + 
    (1 | site_name), 
  data = glm_data, 
  family  = tweedie(link = "log")
)

# Model 2 - Random intercept and slopes 
model_mpa2 <- glmmTMB(
total_urchins ~ heatwave*mpa_status + 
    (1 + year | site_name), 
  data = glm_data, 
  family  = tweedie(link = "log")
 ) 

# Random intercepts and autocorrelation structure 
model_mpa3 <- glmmTMB(
   total_urchins ~ heatwave*mpa_status + 
     (1 | site_name) +
     ar1(0 + year_fct | site_name), 
   data = glm_data, 
   family = tweedie(link = "log")
 )

# DHARMa package 
simulationOutput_mpa_m1 <- simulateResiduals(model_mpa1, plot = F)
plot(simulationOutput_mpa_m1)

simulationOutput_mpa_m2 <- simulateResiduals(model_mpa2, plot = F)
plot(simulationOutput_mpa_m2) 

simulationOutput_mpa_m3 <- simulateResiduals(model_mpa3, plot = F)
plot(simulationOutput_mpa_m3)

AIC(model_mpa1, model_mpa2, model_mpa3)

# Model 2 is chosen 
car::Anova(model_mpa2) # Compares models with vs. without terms to each other 
summary(model_mpa2)
z <- emmeans::emmeans(model_mpa2, pairwise ~ mpa_status | heatwave, type = "response")
z

write.csv(broom::tidy(z$emmeans), 
          "Processed_data/data_tables/emmeans_model_mpa2.csv", 
          row.names = F)
write.csv(broom::tidy(z$contrasts), 
          "Processed_data/data_tables/contrasts_model_mpa2.csv", 
          row.names = F)

##### Central california model
# Random Intercepts 
model_mpa1_central <- glmmTMB(
  total_urchins ~ heatwave*mpa_status + 
    (1 | site_name),
  data = glm_data3, 
  family  = tweedie(link = "log")
)

# Random intercepts and slopes 
model_mpa2_central <- glmmTMB(
  total_urchins ~ heatwave*mpa_status + 
    (1 + year | site_name),
  data = glm_data3, 
  family  = tweedie(link = "log")
)

# Random intercepts and autoccorelation
model_mpa3_central <- glmmTMB(
  total_urchins ~ heatwave*mpa_status + 
    (1 | site_name) +
    ar1(0 + year_fct | site_name),
  data = glm_data3, 
  family  = tweedie(link = "log")
)

# DHARMa package 
simulationOutput_mpa_m1_c <- simulateResiduals(model_mpa1_central, plot = F)
plot(simulationOutput_mpa_m1_c)

simulationOutput_mpa_m2_c <- simulateResiduals(model_mpa2_central, plot = F)
plot(simulationOutput_mpa_m2_c) # Best residuals 

simulationOutput_mpa_m3_c <- simulateResiduals(model_mpa3_central, plot = F)
plot(simulationOutput_mpa_m3_c)

AIC(model_mpa1_central, model_mpa2_central)
# Model 2 has the best residuals and lowest AIC


car::Anova(model_mpa2_central)
summary(model_mpa2_central)
z2 <- emmeans::emmeans(model_mpa2_central, pairwise ~ mpa_status | heatwave, type = "response")
z2

emmeans(model_mpa2_central, pairwise ~ heatwave, type = "response")

write.csv(broom::tidy(z2$emmeans), 
          "Processed_data/data_tables/emmeans_model_mpa2_central.csv", row.names = F)
write.csv(broom::tidy(z2$contrasts), 
          "Processed_data/data_tables/contrasts_model_mpa2_central.csv", 
          row.names = F)

##### GLMM Predator Prey models ################################################
fish_data <- glm_data %>% 
  na.omit(total_urchins) %>% 
  na.omit(SPUL_d)

# Random intercepts with quadratics
model_tc1 <- glmmTMB(
  total_urchins ~ poly(SPUL_d, 2) + poly(PANINT_d, 2) + (1 | site_name), 
  data = fish_data, family = tweedie(link = "log")
) 

# Random intercepts without quadratics
model_tc1_1 <- glmmTMB(
   total_urchins ~ SPUL_d + PANINT_d + (1 | site_name),
   data = fish_data, family = tweedie(link = "log")
)

# random intercepts and sloeps 
#model_tc_rsi <- glmmTMB(
#    total_urchins ~ SPUL_d + PANINT_d + (year_fct | site_name), 
#    data = fish_data, family = tweedie(link = "log")
#  ) #does not converge

# Random intercepts and autocorrelation structure quadratics
model_tc2 <- glmmTMB(
  total_urchins ~ poly(SPUL_d, 2) + poly(PANINT_d, 2) + (1 | site_name) + 
    ar1(0 + year_fct | site_name),
  data = fish_data, family = tweedie(link = "log")
)

# Random intercepts and autocorrelation structure without quadratics 
 model_tc2_1 <- glmmTMB(
   total_urchins ~ SPUL_d + PANINT_d + (1 | site_name) + 
     ar1(0 + year_fct | site_name),
   data = fish_data, family = tweedie(link = "log")
 )
 
simulationOutput_m1 <- simulateResiduals(model_tc1, plot = F)
plot(simulationOutput_m1)

simulationOutput_m1_1 <- simulateResiduals(model_tc1_1, plot = F)
plot(simulationOutput_m1_1) 

simulationOutput_m2 <- simulateResiduals(model_tc2, plot = F)
plot(simulationOutput_m2)

simulationOutput_m2_1 <- simulateResiduals(model_tc2_1, plot = F)
plot(simulationOutput_m2_1)


car::Anova(model_tc1)
summary(model_tc1)

car::Anova(model_tc2)
summary(model_tc2)

# Quadratic models chosen as I assume these effects of predation can level off

##### Model effects plots #####################################################
is_data <- as.data.frame(emmeans::emmeans(model_mpa2, "heatwave", by = "mpa_status", type = "response")) %>% 
  mutate(region = "Southern")
ic_data <- as.data.frame(emmeans::emmeans(model_mpa2_central, "heatwave", by = "mpa_status", type = "response")) %>% 
  mutate(region = "Central")

interaction_data <- rbind(is_data, ic_data)

group.colors <- c(Full = "#440154", Reference = "#FFBA00", Partial ="#21918c")

interaction_plot <- interaction_data %>% 
  ggplot(aes(heatwave, color = mpa_status)) + 
  geom_pointrange(aes(y = response, ymin = asymp.LCL, ymax = asymp.UCL), 
                  position = position_dodge(width = 0.2), 
                  linewidth = 1) +
  facet_grid(cols = vars(region)) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab("Response (NUmber of Urchins)") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  xlab("Heatwave period")
interaction_plot

Urchins_per_region <- data_all %>% 
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern")) %>% 
  group_by(year, mpa_status, region) %>% 
  summarise(total = mean(total_urchins, na.rm = T),
            se_total = std(total_urchins)) %>% 
  ggplot(aes(x = year, y = total, color = mpa_status)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = total - se_total, 
                    ymax = total + se_total), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 0.4) +
  facet_grid(cols = vars(region)) +
  scale_color_manual(values=group.colors, name = "MPA Category", 
                     labels = c("Full", "Partial", "Unprotected")) +
  ylab('Number of Urchins') +
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

##### Autocorrelation #########################################################

## Obtain residuals, pivot to wide data frame
resid_ts <- fish_data |> 
  mutate(resid = resid(model_tc1_1)) |> 
  dplyr::select(site, year, resid) |> 
  pivot_wider(names_from = "site", values_from = "resid")

## Compute normal ACF (partial ACF was giving some weird results)
## Not a problem for AR1 because PACF = ACF at lag 1
acf_df <- bind_rows(apply(
  ts(resid_ts[,-1]), MARGIN = 2, 
  \(x) data.frame(lag = 0:5, acf = as.numeric(acf(x, na.action = na.pass, lag.max = 5, plot = FALSE)$acf))), 
  .id = "site"
)

## Plot ACF
## Suggests very little autocorrelation, on average
violin_plot <- acf_df |> 
  filter(lag > 0) |> 
  mutate(lag = factor(lag)) |> 
  #left_join(data %>% group_by(site = site_name) %>% tally(), by = "site") %>% 
  ggplot(aes(lag, acf)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_violin(scale = "width", alpha = 0.4) +
  geom_violin() +
  geom_jitter(width = 0.2) + 
  geom_smooth() +
  theme_bw()

## Summary of average autocorrelation
acf_df |> group_by(lag) |> summarize(acf = mean(acf, na.rm = TRUE))

# Suggests little autocorrelation

##### Create Publication Figure ################################################
em1 <- ref_grid(model_tc1, at = list(SPUL_d = seq(0, 12, length.out = 100)))
em2 <- ref_grid(model_tc2, at = list(SPUL_d = seq(0, 12, length.out = 100)))

sheephead1 <- as.data.frame(emmip(em1, ~SPUL_d, type = "response", CIs = TRUE, plot = FALSE))
sheephead2 <- as.data.frame(emmip(em2, ~SPUL_d, type = "response", CIs = TRUE, plot = FALSE))

em1 <- ref_grid(model_tc1, at = list(PANINT_d = seq(0, 12, length.out = 100)))
em2 <- ref_grid(model_tc2, at = list(PANINT_d = seq(0, 12, length.out = 100)))

lobster1 <- as.data.frame(emmip(em1, ~PANINT_d, type = "response", CIs = TRUE, plot = FALSE))
lobster2 <- as.data.frame(emmip(em2, ~PANINT_d, type = "response", CIs = TRUE, plot = FALSE))

plotdata <- do.call("rbind", list(
  rename(mutate(sheephead1, model = "random intercepts", species = "Sheephead"), x = SPUL_d),
  rename(mutate(sheephead2, model = "random intercepts + AR(1)", species = "Sheephead"), x = SPUL_d),
  rename(mutate(lobster1, model = "random intercepts", species = "Lobster"), x = PANINT_d),
  rename(mutate(lobster2, model = "random intercepts + AR(1)", species = "Lobster"), x = PANINT_d)))

plot1 <- plotdata |>
  ggplot(aes(x, yvar)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = model), alpha = 0.5) +
  geom_line(aes(color = model)) +
  facet_grid(model~species, scales = "free_x", switch = "x") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()  +
  theme(
    panel.border = element_rect(fill = NA),
    legend.position = "bottom",
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  ) + scale_x_continuous(breaks=c(0,5,10)) +
  ylab("Urchins") 

plot1

dat_text <- data.frame(
  label = c("p < 0.0001; p = 0.0008", "p = 0.0001; p = 0.31", "p < 0.0001; p = 0.035 ", "p = 0.74; p = 0.56"),
  species = c("Sheephead", "Sheephead", "Lobster", "Lobster"),
  model = c("random intercepts", "random intercepts + AR(1)", "random intercepts", "random intercepts + AR(1)")
)

plot_w_labels <- plot1 + geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.1,
  vjust   = -1
)

library(magick)

plot_final <- ggdraw() +
  draw_plot(plot_w_labels) + 
  draw_image(spiny_lobster, x = -0.33, y = 0.44, scale = 0.09) +
  draw_image(spiny_lobster, x = -0.33, y = 0.02, scale = 0.09) +
  draw_image(sheephead_3, x = 0.1, y = 0.45, scale = 0.1) +
  draw_image(sheephead_3, x = 0.1, y = 0.02, scale = 0.1)

##### Export Publication figures ###############################################
# Modeling results plot
png("Figures/GLMM_model_results_PISCO.png", width = 6, height = 6, 
    units = "in", res = 600)
plot_final
dev.off() 

# Modeling interaction between protected areas and heatwave plot
png("Figures/GLMM_heatwave_protection_plot_PISCO.png", width = 8, height = 8,
    units = "in", res = 600)
urchins_and_interactions_plot
dev.off()

# Violin plot
png("Figures/Autocorrelation_PISCO.png", width = 6, height = 6, 
    units = "in", res = 600)
violin_plot
dev.off() 

# Residuals plots for model_tc1 and model_tc2
png("Figures/Residuals_m1_PISCO.png", width = 8, height = 5, 
    units = "in", res = 600)
plot(simulationOutput_m1)
dev.off()

png("Figures/Residuals_m2_ar1_PISCO.png", width = 8, height = 5, 
    units = "in", res = 600)
plot(simulationOutput_m2)
dev.off()

png("Figures/Residuals_mpa_m2_PISCO.png", width = 8, height = 5,
    units = "in", res = 600)
plot(simulationOutput_mpa_m2)
dev.off()

png("Figures/Residuals_mpa_m2_central_PISCO.png", width = 8, height = 5,
    units = "in", res = 600)
plot(simulationOutput_mpa_m2_c)
dev.off()
