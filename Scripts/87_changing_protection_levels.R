# Date: Nov. 11th 2024
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Visualizing MLPA data (abundance graphs)
# BIO 202: Ecological Statistics

###### Set up ##################################################################
# Load packages
library(tidyverse)
library(cowplot)
library(png)
library(patchwork)
library(readxl)

# Load Data
data <- read.csv("Processed_data/MLPA_data_summarized_wo_siteblocks.csv") %>% 
  mutate(mpa_status = ifelse(mpa_status == "Reference", "Unprotected", mpa_status))

MLPA_sites_w_mpas <- read.csv("Processed_data/data_tables/MLPA_sites_with_MPAs.csv")
MPAs_protection_levels <- read_excel("Drafts/Second Submission/Marine_protected_areas_in_manuscript.xlsx")

# Declare Functions 
std <- function(x) sd(x)/sqrt(length(x))

# Declare group colors
group.colors <- c(Full = "#440154", Unprotected = "#FFBA00", Partial ="#21918c")

##### Load organism photos #####################################################
urchins <- readPNG("Figures/Organisms/Urchin_Figure_t.png",native = TRUE)
sheephead <- readPNG("Figures/Organisms/California_Sheephead.png", native = TRUE)
sheephead_3 <- readPNG("Figures/Organisms/Sheephead_3.png", native = TRUE)
spiny_lobster <- readPNG("Figures/Organisms/Spiny_Lobster_Figure_t.png", native = TRUE)

##### Format Data ##############################################################
# Adjust protection levels 
species_protection <- MPAs_protection_levels %>% 
  filter(MLPA == 1) %>% 
  mutate(Sheephead_protected = case_when(Sheephead_protected == 1 ~ "Full",
                                         .default = "Partial"),
         Lobster_protected = case_when(Lobster_protected == 1 ~ "Full",
                                         .default = "Partial"),
         Urchins_protected = case_when(Urchins_protected == 1 ~ "Full",
                                         .default = "Partial"),
         Kelp_protected = case_when(Kelp_protected == 1 ~ "Full",
                                         .default = "Partial")) %>% 
  select(Mpa_ID, Sheephead_protected, Lobster_protected, Urchins_protected, Kelp_protected) %>% 
  rename("Site_ID_12" = "Mpa_ID")

sites_w_protection_status <- MLPA_sites_w_mpas %>% 
  left_join(species_protection, by = "Site_ID_12") %>% 
  replace_na(list(Sheephead_protected = "Unprotected",
                  Lobster_protected = "Unprotected",
                  Urchins_protected = "Unprotected", 
                  Kelp_protected = "Unprotected")) 

data <- data %>% left_join(sites_w_protection_status, by = join_by(site, region))

# Remove NA values
inverts_data <- data %>% filter(!is.na(urchin_total))
fish_data <- data %>% filter(!is.na(biomass_d))

##### Urchins through time #####################################################
Urchins_per_region <- inverts_data %>% 
  mutate(region = ifelse(region == "Central_Coast", "Central", "Southern")) %>% 
  group_by(year, Urchins_protected, region) %>% 
  summarise(total = mean(urchins_d),
            total_se = std(urchins_d)) %>% 
  ggplot(aes(x = year, y = total, color = Urchins_protected)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = total - total_se, 
                    ymax = total + total_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  facet_grid(cols = vars(region)) +
  scale_color_manual(values=group.colors, name = "Protection Status") +
  ylab("Number of Urchnis") +
  ylab(bquote('Urchins per 60' ~ m^2)) +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = c(.1, .75),
        legend.background = element_rect(linewidth = 0.2, colour = 1),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Urchins_per_region



##### Lobster through time #####################################################
lobster_plot <- inverts_data %>% 
  filter(region == "South_Coast") %>% 
  group_by(year, mpa_status.x) %>% 
  summarise(lobster = mean(PANINT_d), 
            lobster_se = std(PANINT_d)) %>% 
  ggplot(aes(x = year, y = lobster, color = mpa_status.x)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = lobster - lobster_se, 
                    ymax = lobster + lobster_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  scale_color_manual(values=group.colors, name = "Protection Status") +
  ylab(bquote('Spiny lobsters per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = c(0.13, 0.78)) +
  inset_element(p = spiny_lobster,
                left = 0.24,
                bottom = 0.75,
                right = 0.41,
                top = 0.95) 

lobster_plot_new <- inverts_data %>% 
  filter(region == "South_Coast") %>% 
  group_by(year, Lobster_protected) %>% 
  summarise(lobster = mean(PANINT_d), 
            lobster_se = std(PANINT_d)) %>% 
  ggplot(aes(x = year, y = lobster, color = Lobster_protected)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = lobster - lobster_se, 
                    ymax = lobster + lobster_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  scale_color_manual(values=group.colors, name = "Protection Status") +
  ylab(bquote('Spiny lobsters per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = c(0.13, 0.78)) +
  inset_element(p = spiny_lobster,
                left = 0.24,
                bottom = 0.75,
                right = 0.41,
                top = 0.95) 
lobster_plot
lobster_plot_new

##### Sheephead through time ###################################################
sheephead_plot <- fish_data %>% 
  filter(region == "South_Coast") %>% 
  group_by(year, mpa_status.x, region) %>% 
  summarise(sheepheads = mean(SPUL_d), 
            sheepheads_se = std(SPUL_d)) %>% 
  ggplot(aes(x = year, y = sheepheads, color = mpa_status.x)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = sheepheads - sheepheads_se, 
                    ymax = sheepheads + sheepheads_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  scale_color_manual(values=group.colors, name = "Protection Status") +
  ylab(bquote('California sheephead per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  inset_element(p = sheephead_3,
                left = 0.04,
                bottom = 0.7,
                right = 0.24,
                top = 0.95)

sheephead_plot_new <- fish_data %>% 
  filter(region == "South_Coast") %>% 
  group_by(year, Sheephead_protected, region) %>% 
  summarise(sheepheads = mean(SPUL_d), 
            sheepheads_se = std(SPUL_d)) %>% 
  ggplot(aes(x = year, y = sheepheads, color = Sheephead_protected)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = sheepheads - sheepheads_se, 
                    ymax = sheepheads + sheepheads_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  scale_color_manual(values=group.colors, name = "Protection Status") +
  ylab(bquote('California sheephead per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  inset_element(p = sheephead_3,
                left = 0.04,
                bottom = 0.7,
                right = 0.24,
                top = 0.95)

sheephead_plot
sheephead_plot_new


biomass_plot <- fish_data %>% 
  filter(region == "South_Coast") %>% 
  group_by(year, mpa_status.x, region) %>% 
  summarise(biomass = mean(biomass_d/1000), 
            biomass_se = std(biomass_d/1000)) %>% 
  ggplot(aes(x = year, y = biomass, color = mpa_status.x)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = biomass - biomass_se, 
                    ymax = biomass + biomass_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  scale_color_manual(values=group.colors, name = "Protection Status") +
  ylab(bquote('Sheephead biomass (kg) per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  inset_element(p = sheephead,
                left = 0.04,
                bottom = 0.7,
                right = 0.24,
                top = 0.95) 

biomass_plot_new <- fish_data %>% 
  filter(region == "South_Coast") %>% 
  group_by(year, Sheephead_protected, region) %>% 
  summarise(biomass = mean(biomass_d/1000), 
            biomass_se = std(biomass_d/1000)) %>% 
  ggplot(aes(x = year, y = biomass, color = Sheephead_protected)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = biomass - biomass_se, 
                    ymax = biomass + biomass_se), 
                width = 0.2, 
                linewidth = .5, 
                alpha = 1) +
  scale_color_manual(values=group.colors, name = "Protection Status") +
  ylab(bquote('Sheephead biomass (kg) per 60' ~ m^2)) +
  xlab("Year") +
  annotate("rect", fill = "red", alpha = 0.2, 
           xmin = 2014, xmax = 2016,
           ymin = -Inf, ymax = Inf) +
  geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 0.7) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  inset_element(p = sheephead,
                left = 0.04,
                bottom = 0.7,
                right = 0.24,
                top = 0.95) 

biomass_plot
biomass_plot_new

combo_plot <- plot_grid(lobster_plot, 
                        sheephead_plot,
                        biomass_plot, 
                        labels = c("A", "B", "C"), ncol = 1, hjust = 0.1)

combo_plot_new <- plot_grid(lobster_plot_new, 
                        sheephead_plot_new,
                        biomass_plot_new, 
                        labels = c("D", "E", "F"), ncol = 1, hjust = 0.1)

mega_plot <- plot_grid(combo_plot, combo_plot_new, nrow = 1)


##### Export ###################################################################
png("Figures/Sheephead_lobsters_MLPA_mega.png", width = 12, height = 10, 
    units = "in", res = 600)
mega_plot
dev.off()

##### Statistical Analyses for abundances through time and protection ##########
# load packages
library(glmmTMB)
library(emmeans)
library(DHARMa)


# format data, start from data because it will be treated the same as script 86
fish_post2012 <- data %>% 
  filter(region == "South_Coast") %>% 
  filter(year >= 2012) %>% 
  mutate(year_fct = as.factor(year),
         year_std = (year - 2017)/sd(year)) %>% 
  mutate(Sheephead_protected = factor(Sheephead_protected, levels = c("Unprotected", "Partial", "Full"))) 

inverts_post2012 <- data %>% 
  filter(region == "South_Coast") %>% 
  filter(year >= 2012) %>% 
  mutate(year_fct = as.factor(year),
         year_std = (year - 2017)/sd(year)) %>% 
  mutate(Lobster_protected = factor(Lobster_protected, levels = c("Unprotected", "Partial", "Full")))


# sheephead
model_sheephead <- glmmTMB(
  SPUL_d ~ Sheephead_protected + 
    (1 | site) + # random intercept 
    ar1(0 + year_fct | site), #  autoregressive function 
  data = fish_post2012, family = tweedie(link = "log")
)

model_sheephead
simulationOutput_sheephead <- simulateResiduals(model_sheephead, plot = F)
plot(simulationOutput_sheephead)

car::Anova(model_sheephead)
plot(emmeans(model_sheephead, specs = pairwise ~ Sheephead_protected, type = "response"))
summary(model_sheephead)

# sheephead biomass
model_biomass <- glmmTMB(
  biomass_d ~ Sheephead_protected + 
    (1 + year_std | site) + 
    ar1(0 + year_fct | site),
  data = fish_post2012, family = tweedie(link = "log")
)

model_biomass
simulationOutput_biomass <- simulateResiduals(model_biomass, plot = F)
plot(simulationOutput_biomass)

car::Anova(model_biomass)
plot(emmeans(model_biomass, specs = pairwise ~ Sheephead_protected, type = "response"))
summary(model_biomass)

# lobsters
model_lobsters <- glmmTMB(
  PANINT_d ~ Lobster_protected + 
    (1 + year_std | site) + # random intercept 
    ar1(0 + year_fct | site), # autoregressive function 
  data = inverts_post2012, family = tweedie(link = "log")
)

model_lobsters
simulationOutput_lobsters <- simulateResiduals(model_lobsters, plot = F)
plot(simulationOutput_lobsters)

car::Anova(model_lobsters)
plot(emmeans(model_lobsters, specs = pairwise ~ Lobster_protected, type = "response"))
summary(model_lobsters)

