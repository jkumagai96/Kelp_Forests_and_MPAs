# Date: October 26th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Run stats on PISCO data to compare with reefcheck data
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(png)
library(patchwork)
library(cowplot)

# Load Data
PISCO_data <- read.csv("Data/PISCO subtidal data clearinghouse/PISCO_kelpforest_combined_data_bysite.csv")
sites_for_joining <- read.csv("Processed_data/data_tables/PISCO_sites_with_MPAs.csv") 
PISCO_biomass_data <- read.csv("Processed_data/Data_tables/PISCO_biomass_data.csv")

# Declare Functions
std <- function(x) sd(x)/sqrt(length(x))

# Add in little black outlines of organisms
urchins <- readPNG("Figures/Organisms/Urchin_Figure_t.png",native = TRUE)
sheephead <- readPNG("Figures/Organisms/California_Sheephead.png", native = TRUE)
sheephead_3 <- readPNG("Figures/Organisms/Sheephead_3.png", native = TRUE)
spiny_lobster <- readPNG("Figures/Organisms/Spiny_Lobster_Figure_t.png", native = TRUE)

##### Filter Data ##############################################################
### Clean Data
data <- PISCO_data %>% 
  left_join(sites_for_joining, by = "SITE") %>% 
  dplyr::select(SITE,
                LATITUDE,
                LONGITUDE,
                region,
                mpa_status,
                Estab_Yr_1,
                year,
                fish_SPUL,
                swath_PANINT,
                swath_CENCOR,
                swath_MESFRAAD,
                swath_STRPURAD
                ) %>% 
  filter(year >= 2007,
         region == "South_Coast") %>% 
  mutate(total_urchins = rowSums(across(c(swath_CENCOR,
                                 swath_MESFRAAD,
                                 swath_STRPURAD)))) %>% 
  drop_na(fish_SPUL, swath_PANINT, swath_CENCOR, swath_MESFRAAD, swath_STRPURAD) %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("Reference", "Partial", "Full"))) %>% 
  filter(total_urchins >= 0)

##### Explore data #############################################################
group.colors <- c(Full = "#440154", Reference = "#FFBA00", Partial ="#21918c")

# Urchins plot 
urchin_data <- data %>% 
  group_by(year, mpa_status) %>%
  summarize(avg_urchins = mean(total_urchins),
            se_urchins = std(total_urchins))

urchins_plot <- ggplot(aes(x = year, y = avg_urchins, color = mpa_status),
                       data = urchin_data) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  geom_errorbar(aes(ymin = avg_urchins - se_urchins, 
                    ymax = avg_urchins + se_urchins), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 0.4) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab(bquote('Urchins per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  theme_bw() +
  theme(legend.position = c(.85, .75)) +
  inset_element(p = urchins,
                left = 0.1,
                bottom = 0.7,
                right = 0.24,
                top = 0.95) 
urchins_plot

# Sheephead plot 
sheephead_plot <- data %>% 
  group_by(year, mpa_status) %>%
  summarize(avg_sheephead = mean(fish_SPUL),
            se_sheephead = std(fish_SPUL)) %>% 
  ggplot(aes(x = year, y = avg_sheephead, color = mpa_status)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  geom_errorbar(aes(ymin = avg_sheephead - se_sheephead, 
                    ymax = avg_sheephead + se_sheephead), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 0.4) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab(bquote('Sheephead per 120' ~ m^3)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  theme_bw() +
  theme(legend.position = "none") +
  inset_element(p = sheephead_3,
                left = 0.1,
                bottom = 0.7,
                right = 0.24,
                top = 0.95)
sheephead_plot

# Lobster Plot 
lobster_data <- data %>% 
  group_by(year, mpa_status) %>%
  summarize(avg_lobster = mean(swath_PANINT),
            se_lobster = std(swath_PANINT))

lobster_plot <- ggplot(aes(x = year, y = avg_lobster, color = mpa_status),
                       data = lobster_data) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  geom_errorbar(aes(ymin = avg_lobster - se_lobster, 
                    ymax = avg_lobster + se_lobster), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 0.4) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab(bquote('Lobsters per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  theme_bw() +
  theme(legend.position = "none") +
  inset_element(p = spiny_lobster,
                left = 0.07,
                bottom = 0.7,
                right = 0.27,
                top = 0.95)
lobster_plot

##### Create biomass figure ####################################################
biomass_plot <- ggplot(data = PISCO_biomass_data, 
                       aes(x = year,
                           y = avg_biomass,
                           color = mpa_status)) +
  geom_line(aes(color = mpa_status), linewidth = 1) +
  geom_errorbar(aes(ymin = avg_biomass-se_biomass, ymax = avg_biomass + se_biomass),
                width = 0.2, 
                linewidth = .5, 
                alpha = 0.4)+
  theme_bw() +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  xlab("Year") +
  ylab(bquote('Sheephead Biomass per 120' ~ m^3)) +
  theme(legend.position = "none") +
  annotate("rect", fill = "red", alpha = 0.2,
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = .7) +
  inset_element(p = sheephead,
                left = 0.1,
                bottom = 0.7,
                right = 0.24,
                top = 0.95)
biomass_plot

combo_plot <- plot_grid(urchins_plot, sheephead_plot, lobster_plot, biomass_plot,
                        labels = "AUTO", nrow = 2)
combo_plot
##### Export ###################################################################
png("Figures/Urchins_sheephead_lobsters_PISCO.png", width = 12, height = 8, 
    units = "in", res = 600)
combo_plot
dev.off() 

##### Statistics ###############################################################
# Are the differences seen in the urchin and lobster plots significant?
# Simple ANOVA on lobster
lobster_data <- lobster_data %>% 
  filter(year >= 2012) # after all the MPAs have been implemented

l_aov <- aov(avg_lobster ~ mpa_status,
               data = lobster_data
)
# Are the residuals normal? 
hist(l_aov$residuals)

# QQ-plot
library(car)
qqPlot(l_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)

summary(l_aov)

# Tukey HSD test:
library(multcomp)
post_test <- glht(l_aov,
                  linfct = mcp(mpa_status = "Tukey")
)

summary(post_test)


# Simple ANOVA on urchins
urchin_data <- urchin_data %>% 
  filter(year >= 2012) # after all the MPAs have been implemented

u_aov <- aov(avg_urchins ~ mpa_status,
             data = urchin_data
)
# Are the residuals normal? 
hist(u_aov$residuals)

# QQ-plot
library(car)
qqPlot(u_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)

summary(u_aov)

# Tukey HSD test:
library(multcomp)
post_test <- glht(u_aov,
                  linfct = mcp(mpa_status = "Tukey")
)

summary(post_test)

##### Linear Regression Statistics #############################################
library("MASS")
library("ggeffects")
library("nlme")
library("emmeans")

glm_data <- data %>% 
  mutate(
    heatwave = case_when(year < 2014 ~ "before", year > 2016 ~ "after", .default = "during"),
    heatwave = factor(heatwave, levels = c("before", "during", "after")), 
    site_name = factor(SITE), 
    year_fct = factor(year)
  )

glm_data$total_urchins[glm_data$total_urchins == 0] <- min(glm_data$total_urchins[glm_data$total_urchins != 0])

glm_data2 %>% 
  left_join(PISCO_biomass_data, by = )


# Model 1 - Random intercepts
model_mpa1 <- glmmPQL(
  total_urchins ~ heatwave*mpa_status, 
  random = ~ 1 | site_name, 
  data = glm_data, family = Gamma(link = "log")
)

model_mpa2 <- glmmPQL(
  total_urchins ~ heatwave*mpa_status, 
  random = ~ 1 + year| site_name, 
  data = glm_data, family = Gamma(link = "log")
) # Does not converge

# Random intercepts and autocorrelation structure 
model_mpa3 <- glmmPQL(
  total_urchins ~ heatwave*mpa_status, 
  random = ~ 1 | site_name, 
  correlation = corAR1(form = ~ year | site_name),
  data = glm_data, family = Gamma(link = "log")
)

plot(model_mpa1) # Plot residuals 
car::Anova(model_mpa1) # Compares models with vs. without terms to each other 
summary(model_mpa1)

plot(model_mpa3)
car::Anova(model_mpa3)
summary(model_mpa3)

# Perhaps we can conclude that the interaction between the heatwave and protection status is highly significant

# Random intercepts
model_tc1 <- glmmPQL(
  total_urchins ~ fish_SPUL + swath_PANINT,
  random = ~ 1 | site_name, 
  data = glm_data, family = Gamma(link = "log")
)

# Random intercepts and autocorrelation structure
model_tc2 <- glmmPQL(
  total_urchins ~ fish_SPUL + swath_PANINT,
  random = ~ 1 | site_name, 
  correlation = corAR1(form = ~ year | site_name),
  data = glm_data, family = Gamma(link = "log")
)

# Random intercepts and slopes (does not converge)

plot(model_tc1)
car::Anova(model_tc1)
summary(model_tc1)


plot(model_tc2)
car::Anova(model_tc2)
summary(model_tc2)

##### Model effects plots #####################################################
interaction_data <- as.data.frame(emmeans::emmeans(model_mpa1, "heatwave", by = "mpa_status", type = "response"))

interaction_plot <- interaction_data %>% 
  ggplot(aes(heatwave, color = mpa_status)) + 
  geom_pointrange(aes(y = response, ymin = lower.CL, ymax = upper.CL), 
                  position = position_dodge(width = 0.2), 
                  linewidth = 1) +
  scale_color_manual(values=group.colors, name = "MPA Category") +
  ylab("Response (Urchins)") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank()) 
interaction_plot

##### Autocorrelation #########################################################

## Obtain residuals, pivot to wide data frame
resid_ts <- data |> 
  mutate(resid = resid(model_tc1)) |> 
  dplyr::select(SITE, year, resid) |> 
  pivot_wider(names_from = "SITE", values_from = "resid")

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

# Suggests little autocorrelation to me 

##### Create Publication Figure ################################################
em1 <- ref_grid(model_tc1, at = list(fish_SPUL = seq(0, max(glm_data$fish_SPUL), length.out = 100)))
em2 <- ref_grid(model_tc2, at = list(fish_SPUL = seq(0, max(glm_data$fish_SPUL), length.out = 100)))

sheephead1 <- as.data.frame(emmip(em1, ~fish_SPUL, type = "response", CIs = TRUE, plot = FALSE))
sheephead2 <- as.data.frame(emmip(em2, ~fish_SPUL, type = "response", CIs = TRUE, plot = FALSE))

em1 <- ref_grid(model_tc1, at = list(swath_PANINT = seq(0, max(glm_data$swath_PANINT), length.out = 100)))
em2 <- ref_grid(model_tc2, at = list(swath_PANINT = seq(0, max(glm_data$swath_PANINT), length.out = 100)))

lobster1 <- as.data.frame(emmip(em1, ~swath_PANINT, type = "response", CIs = TRUE, plot = FALSE))
lobster2 <- as.data.frame(emmip(em2, ~swath_PANINT, type = "response", CIs = TRUE, plot = FALSE))


plotdata <- do.call("rbind", list(
  rename(mutate(sheephead1, model = "random intercepts", species = "Sheephead"), x = fish_SPUL),
  rename(mutate(sheephead2, model = "random intercepts + AR(1)", species = "Sheephead"), x = fish_SPUL),
  rename(mutate(lobster1, model = "random intercepts", species = "Lobster"), x = swath_PANINT),
  rename(mutate(lobster2, model = "random intercepts + AR(1)", species = "Lobster"), x = swath_PANINT)
))

plot1 <- plotdata |>
  ggplot(aes(x, yvar)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = model), alpha = 0.5) +
  geom_line(aes(color = model)) +
  facet_grid(model~species, scales = "free_x", switch = "x") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.position = "bottom",
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  ) +
  ylab("Urchins")

dat_text <- data.frame(
  label = c("p < .0005", "p = 0.0028", "p < 0.0005", "p = 0.6605"),
  species = c("Sheephead", "Sheephead", "Lobster", "Lobster"),
  model = c("random intercepts", "random intercepts + AR(1)", "random intercepts", "random intercepts + AR(1)")
)

plot_final <- plot1 + geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.1,
  vjust   = -1
)
plot_final

# Modeling results plot
png("Figures/GLMM_model_results_PISCO.png", width = 6, height = 6, 
    units = "in", res = 600)
plot_final
dev.off() 

# Modeling iteraction between protected areas and heatwave plot
png("Figures/GLMM_heatwave_protection_plot_PISCO.png", width = 5, height = 4,
    units = "in", res = 600)
interaction_plot
dev.off()

# Violin plot
png("Figures/Autocorrelation_PISCO.png", width = 6, height = 6, 
    units = "in", res = 600)
violin_plot
dev.off() 

