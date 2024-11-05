# Date: July 25th 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Explore correlations between kelp and urchins and their predators
# BIO 202: Ecological Statistics

##### Set up ###################################################################
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
# Load Data

data_all <- read.csv("Processed_data/MLPA_data_summarized_wo_siteblocks.csv") %>% 
  mutate(site = factor(site))

# Add in little black outlines of organisms
urchin_sheephead <- readPNG("Figures/Organisms/urchin_sheephead.png", native = TRUE)
urchin_lobster <- readPNG("Figures/Organisms/urchin_lobster.png", native = TRUE)
kelp_urchin <- readPNG("Figures/Organisms/kelp_urchin.png", native = TRUE)
kelp_sheephead <- readPNG("Figures/Organisms/kelp_sheephead.png", native = TRUE)
kelp_lobster <- readPNG("Figures/Organisms/kelp_lobster.png", native = TRUE)

##### Exploring ################################################################

data_all %>% 
  filter(region == "South_Coast") %>% 
  ggplot(aes(x = urchins_d, y = MACPYRAD_d)) +
  geom_point()

data_all %>% 
  filter(region == "Central_Coast") %>% 
  ggplot(aes(x = urchins_d, y = MACPYRAD_d)) +
  geom_point()

data_all %>% 
  ggplot(aes(x = urchins_d, y = MACPYRAD_d)) +
  geom_point()

data_all %>% 
  filter(region == "South_Coast") %>% 
  ggplot(aes(x = SPUL_d, y = MACPYRAD_d)) +
  geom_point()

data_all %>% 
  filter(region == "South_Coast") %>% 
  ggplot(aes(x = PANINT_d, y = MACPYRAD_d)) +
  geom_point()

##### Statistics ###############################################################
glm_data_all <- data_all %>% filter(!is.na(urchin_total)) %>% 
  mutate(year_fct = factor(year))
glm_data_south <- data_all %>% filter(!is.na(urchin_total),
                                      !is.na(SPUL_d),
                                      region == "South_Coast") %>% 
  mutate(year_fct = factor(year))

glm_data_central <- data_all %>% filter(!is.na(urchin_total),
                                      !is.na(SPUL_d),
                                      region == "Central_Coast") %>% 
  mutate(year_fct = factor(year))


## Giant kelp and urchin correlation
# Southern California 
model_giant1_s <- glmmTMB(
  MACPYRAD_d ~ poly(urchins_d, 2) +
    (1 | site), 
  data = glm_data_south, family = tweedie(link = "log")
) 

# model_giant1_s <- glmmTMB(
#   MACPYRAD_d ~ poly(urchins_d, 2) +
#     (1 + year_fct | site), 
#   data = glm_data_south, family = tweedie(link = "log")
# )
# Convergence issues 

model_giant1_s_ar1 <- glmmTMB(
  MACPYRAD_d ~ poly(urchins_d, 2) +
    (1 | site) +
    ar1(0 + year_fct | site),
  data = glm_data_south, family = tweedie(link = "log")
) 


simulationOutput_giant1_s <- simulateResiduals(model_giant1_s, plot = F)
plot(simulationOutput_giant1_s) # Residuals are better 

simulationOutput_giant1_s_ar1 <- simulateResiduals(model_giant1_s_ar1, plot = F)
plot(simulationOutput_giant1_s_ar1)


summary(model_giant1_s)

sigma_re_giant1_s <- sqrt(as.numeric(glmmTMB::VarCorr(model_giant1_s)$cond$site))
em1 <- ref_grid(model_giant1_s, at = list(urchins_d = seq(0, 1000, length.out = 100)))

giantkelp1_s <- as.data.frame(emmip(em1, ~urchins_d, type = "response", CIs = TRUE, plot = FALSE, bias.adjust = TRUE, sigma = sigma_re_giant1_s))

plotdata <- do.call("rbind", list(
  rename(mutate(giantkelp1_s, x = urchins_d))))

plot1_s <- plotdata |>
  ggplot(aes(x, yvar)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "darkorchid4", alpha = 0.5) +
  geom_line(color = "darkorchid4") +
  theme_minimal() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.position = "bottom",
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  ylab(bquote('Giant kelp per 60 ' ~m^2)) +
  xlab(bquote('Urchins per 60 ' ~m^2)) 

plot1_s

# Central California
# model_giant1_c <- glmmTMB(
#   MACPYRAD_d ~ poly(urchins_d, 2) +
#     (1 | site), 
#   data = glm_data_central, family = tweedie(link = "log")
# ) 
# 
# simulationOutput_giant1_c <- simulateResiduals(model_giant1_c, plot = F)
# plot(simulationOutput_giant1_c)
# 
# summary(model_giant1_c)
# 
# sigma_re_giant1_c <- sqrt(as.numeric(glmmTMB::VarCorr(model_giant1_c)$cond$site))
# em1 <- ref_grid(model_giant1_c, at = list(urchins_d = seq(0, 1000, length.out = 100)))
# 
# giantkelp1_c <- as.data.frame(emmip(em1, ~urchins_d, type = "response", CIs = TRUE, plot = FALSE, adjust.bias = TRUE, sigma = sigma_re_giant1_c))
# 
# plotdata <- do.call("rbind", list(
#   rename(mutate(giantkelp1_c, x = urchins_d))))
# 
# plot1_c <- plotdata |>
#   ggplot(aes(x, yvar)) +
#   geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "darkorchid4", alpha = 0.5) +
#   geom_line(color = "darkorchid4") +
#   theme_minimal()  +
#   theme(
#     panel.border = element_rect(fill = NA),
#     legend.position = "bottom",
#     strip.placement = "outside",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()) +
#   ylab("Giant Kelp") +
#   xlab("Urchins")
# 
# plot1_c


# Indirect trophic cascade: 
## Giant kelp and predators correlation
model_giant2 <- glmmTMB(
  MACPYRAD_d ~ SPUL_d + PANINT_d +
    (1 | site), 
  data = glm_data_south, family = tweedie(link = "log")
) # Between nbinom1, nbinom2, and poisson. nbinom1 gives the best residuals

model_giant2_ar1 <- glmmTMB(
  MACPYRAD_d ~ SPUL_d + PANINT_d +
    (1 | site) +
    ar1(0 + year_fct | site), 
  data = glm_data_south, family = tweedie(link = "log")
) 

# model_giant2 <- glmmTMB(
#   MACPYRAD_d ~ SPUL_d + PANINT_d +
#     (1 + year_fct | site), 
#   data = glm_data_south, family = tweedie(link = "log")
# ) # Convergence issues 

model_giant2b <- glmmTMB(
  MACPYRAD_d ~ poly(SPUL_d, 2) + poly(PANINT_d, 2) +
    (1 | site), 
  data = glm_data_south, family = tweedie(link = "log")
) 

AIC(model_giant2, model_giant2b) # 2b has the lower AIC, but I don't think it makes sense
# to have a quadratic here... I think the simplest hypothesis is just a + relationship 

simulationOutput_giant2 <- simulateResiduals(model_giant2, plot = F)
plot(simulationOutput_giant2)

simulationOutput_giant2_ar1 <- simulateResiduals(model_giant2_ar1, plot = F)
plot(simulationOutput_giant2_ar1)


simulationOutput_giant2b <- simulateResiduals(model_giant2b, plot = F)
plot(simulationOutput_giant2b)

summary(model_giant2)
summary(model_giant2b)
summary(model_giant2_ar1) # Does not change the interpretation of these relationships 

sigma_re_giant2 <- sqrt(as.numeric(glmmTMB::VarCorr(model_giant2)$cond$site))
em1 <- ref_grid(model_giant2, at = list(SPUL_d = seq(0, 12, length.out = 100)))
sheephead1 <- as.data.frame(emmip(em1, ~SPUL_d, type = "response", CIs = TRUE, plot = FALSE, 
                                  adjust.bias = T, sigma = sigma_re_giant2))

em1 <- ref_grid(model_giant2, at = list(PANINT_d = seq(0, 12, length.out = 100)))
lobster1 <- as.data.frame(emmip(em1, ~PANINT_d, type = "response", CIs = TRUE, plot = FALSE,
                                adjust.bias = T, sigma = sigma_re_giant2))

plotdata2 <- do.call("rbind", list(
  rename(mutate(sheephead1, model = "random intercepts", species = "Sheephead"), x = SPUL_d),
  rename(mutate(lobster1, model = "random intercepts", species = "Lobsters"), x = PANINT_d)))


plotdata_x <- transform(plotdata2, species = factor(species, levels=c("Lobsters", "Sheephead"),
                                    labels=c("Lobster~per~60~m^2", "Sheephead~per~60~m^2")))
plot2 <- plotdata_x |>
  ggplot(aes(x, yvar)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = species), alpha = 0.5) +
  geom_line(aes(color = species)) +
  facet_grid(~species, 
             scales = "free_x", switch = "x", 
             labeller = label_parsed) +
  scale_color_manual(values = c("#D55E00", "#D55E00")) +
  scale_fill_manual(values = c("#D55E00", "#D55E00")) +
  theme_minimal()  +
  theme(
    panel.border = element_rect(fill = NA),
    legend.position = "none",
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  ) + scale_x_continuous(breaks=c(0,5,10)) +
  ylab(bquote('Giant kelp per 60 ' ~m^2)) 

plot2

## Urchins as predicted by predators
model_pp <- glmmTMB(
  urchin_total ~ poly(SPUL_d, 2) + poly(PANINT_d, 2) + 
    offset(log(n_transects)) + 
    (1 | site), 
  data = glm_data_south, family = nbinom1(link = "log")
) 

simulationOutput_pp <- simulateResiduals(model_pp, plot = F)
plot(simulationOutput_pp) 

summary(model_pp)

# Test fit with tweeide family 
model_pp2 <- glmmTMB(
  urchins_d ~ poly(SPUL_d, 2) + poly(PANINT_d, 2) + 
    (1 | site), 
  data = glm_data_south, family = tweedie(link = "log")
) 

simulationOutput_pp2 <- simulateResiduals(model_pp2, plot = F)
plot(simulationOutput_pp2) 

# Test random slopes 
# model_pp3 <- glmmTMB(
#   urchin_total ~ poly(SPUL_d, 2) + poly(PANINT_d, 2) + 
#     offset(log(n_transects)) + 
#     (1 + year_fct | site), 
#   data = glm_data_south, family = nbinom1(link = "log")) 
# convergence issues 


# Extract standard deviation of random effect
sigma_re <- sqrt(as.numeric(glmmTMB::VarCorr(model_pp)$cond$site))

# Calculate estimated marginal means 
em1 <- ref_grid(model_pp, at = list(SPUL_d = seq(0, 12, length.out = 100)), offset = log(1))

sheephead1 <- as.data.frame(emmip(em1, ~SPUL_d, type = "response", CIs = TRUE, plot = FALSE, bias.adjust = TRUE, sigma = sigma_re))

em1 <- ref_grid(model_pp, at = list(PANINT_d = seq(0, 12, length.out = 100)), offset = log(1))

lobster1 <- as.data.frame(emmip(em1, ~PANINT_d, type = "response", CIs = TRUE, plot = FALSE, bias.adjust = TRUE, sigma = sigma_re))

plotdata3 <- do.call("rbind", list(
  rename(mutate(sheephead1, model = "random intercepts", species = "Sheephead"), x = SPUL_d),
  rename(mutate(lobster1, model = "random intercepts", species = "Lobsters"), x = PANINT_d)))

plotdata_3 <- transform(plotdata3, species = factor(species, levels=c("Lobsters", "Sheephead"),
                                                    labels=c("Lobster~per~60~m^2", "Sheephead~per~60~m^2")))


plot3 <- plotdata_3 |>
  ggplot(aes(x, yvar)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = species), alpha = 0.5) +
  geom_line(aes(color = species)) +
  facet_grid(~species, scales = "free_x", switch = "x", 
             labeller = label_parsed) +
  scale_color_manual(values = c("#0072B2", "#0072B2")) +
  scale_fill_manual(values = c("#0072B2", "#0072B2")) +
  theme_minimal()  +
  theme(
    panel.border = element_rect(fill = NA),
    legend.position = "none",
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  ) + scale_x_continuous(breaks=c(0,5,10)) +
  ylab(bquote('Urchins per 60 ' ~m^2)) 

plot3

### Reporting relationships

# Number of giant kelp when there are no urchins
plotdata %>% filter(x == 0)
plotdata %>% filter(abs(x) == x[which.min(abs(x - median(glm_data_south$urchins_d)))])

# Number of urchins when there are 5 lobsters or sheephead when other predator is held constant 
plotdata3 %>% filter(species == "Lobsters", x == 0)
plotdata3 %>% filter(species == "Sheephead", x == 0)

plotdata3 %>% filter(species == "Lobsters", abs(x) == x[which.min(abs(x - 5))])
plotdata3 %>% filter(species == "Sheephead", abs(x) == x[which.min(abs(x - 5))])

# Number of giant kelp when lobsters are at 0 and 5, and sheephead is held constant at mean
plotdata2 %>% filter(species == "Lobsters", x == 0)
plotdata2 %>% filter(species == "Lobsters", abs(x) == x[which.min(abs(x - 5))])

### Double checking urchin model
mean((glm_data_south$urchin_total/glm_data_south$n_transects)[glm_data_south$PANINT_d == 0])

glm_data_south$urchin_fit <- predict(model_pp, type = "response")
mean((glm_data_south$urchin_fit/glm_data_south$n_transects)[glm_data_south$PANINT_d == 0])

##### Plottting #####
library(cowplot)
plot_raw <- plot_grid(plot1_s, plot3, plot2, nrow = 3, labels = "AUTO")

plot_final <- ggdraw() +
  draw_plot(plot_raw) + 
  draw_image(urchin_lobster, x = -0.34, y = -0.06, scale = 0.11) +
  draw_image(urchin_sheephead, x = 0.12, y = -0.06, scale = 0.11) +
  draw_image(kelp_urchin, x = -0.34, y = 0.27, scale = 0.11) + # giant kelp vs. urchins
  draw_image(kelp_lobster, x = -0.34, y = -0.33, scale = 0.11)  + # giant kelp vs. lobsters
  draw_image(kelp_sheephead, x = 0.12, y = -0.35, scale = 0.11) + # giant kelp vs. sheephead 
  annotate("text", x = 0.86, y = 0.96, label = "p < 0.0001; p = 0.1", fontface = "bold") +
  annotate("text", x = 0.39, y = 0.63, label = "p < 0.0001; p < 0.0001", fontface = "bold") +
  annotate("text", x = 0.85, y = 0.63, label = "p < 0.0001; p = 0.059", fontface = "bold") +
  annotate("text", x = 0.43, y = 0.3, label = "p < 0.0001", fontface = "bold") +
  annotate("text", x = 0.92, y = 0.3, label = "p = 0.8")

png(plot_final, filename = "Figures/Trophic_cascade.png", units = "in", 
    width = 6, height = 7.5, res = 600)
plot_final
dev.off()

##### Export Residual Plots ####################################################
# Southern California urchins vs. giant kelp 
simulationOutput_giant1_s <- simulateResiduals(model_giant1_s, plot = F)

png(filename = "Figures/Residuals_kelp_v_urchins.png", 
    width = 8, 
    height = 5,
    units = "in", 
    res = 600)
plot(simulationOutput_giant1_s)
dev.off()

# Urchins vs. predators 
simulationOutput_pp <- simulateResiduals(model_pp, plot = F)

png(filename = "Figures/Residuals_predators_v_urchins.png", 
    width = 8, 
    height = 5,
    units = "in", 
    res = 600)
plot(simulationOutput_pp) 
dev.off()

# Predators vs. kelp

simulationOutput_giant2 <- simulateResiduals(model_giant2, plot = F)

png(filename = "Figures/Residuals_predators_v_kelp.png", 
    width = 8, 
    height = 5,
    units = "in", 
    res = 600)
plot(simulationOutput_giant2)
dev.off()

