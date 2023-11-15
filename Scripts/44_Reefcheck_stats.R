# Date: August 14th 2023
# Author: Maruice Goodman and Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Run GLMs on clean reef check data
# BIO 202: Ecological Statistics

##### Load Data and Packages ###################################################
library("MASS")
library("tidyverse")
library("ggeffects")
library("nlme")
library("emmeans")

data <- read.csv("Processed_data/data_tables/subtidal_surveys.csv") 

##### Format Data ##############################################################

data <- data |>  
  filter(region == "South_Coast") |>
  dplyr::select(site_name, year, mpa_status, avg_urchins, avg_sheephead, avg_lobster) |>  
  group_by(site_name, year, mpa_status) |> 
  summarize_all(mean) |> ungroup() |> 
  mutate(
    mpa_status = factor(mpa_status, levels = c("None", "Partial", "Full")),
    heatwave = case_when(year < 2014 ~ "before", year > 2016 ~ "after", .default = "during"),
    heatwave = factor(heatwave, levels = c("before", "during", "after")), 
    site_name = factor(site_name), 
    year_fct = factor(year)
  )

data$avg_urchins[data$avg_urchins == 0] <- min(data$avg_urchins[data$avg_urchins != 0])

##### Fit Models ###############################################################
# Random intercepts
model_mpa1 <- glmmPQL(
  avg_urchins ~ heatwave*mpa_status, 
  random = ~ 1 | site_name, 
  data = data, family = Gamma(link = "log")
)

model_mpa2 <- glmmPQL(
  avg_urchins ~ heatwave*mpa_status, 
  random = ~ 1 + year| site_name, 
  data = data, family = Gamma(link = "log")
) # Does not converge

# Random intercepts and autocorrelation structure 
model_mpa3 <- glmmPQL(
  avg_urchins ~ heatwave*mpa_status, 
  random = ~ 1 | site_name, 
  correlation = corAR1(form = ~ year | site_name),
  data = data, family = Gamma(link = "log")
) # Same significance levels as model_mpa1 

# Now model trophic interactions, don't include MPAs because we don't want to isolate 
# the effect of MPAs on these interactions, but rather see if more lobsters = less urchins, and same with sheephead 

# Random intercepts
model_tc1 <- glmmPQL(
  avg_urchins ~ avg_sheephead + avg_lobster,
  random = ~ 1 | site_name, 
  data = data, family = Gamma(link = "log")
)

# Random intercepts and autocorrelation structure
model_tc2 <- glmmPQL(
  avg_urchins ~ avg_sheephead + avg_lobster,
  random = ~ 1 | site_name, 
  correlation = corAR1(form = ~ year | site_name),
  data = data, family = Gamma(link = "log")
)

##### Autocorrelation #########################################################

## Obtain residuals, pivot to wide data frame
resid_ts <- data |> 
  mutate(resid = resid(model_tc1)) |> 
  dplyr::select(site_name, year, resid) |> 
  pivot_wider(names_from = "site_name", values_from = "resid")

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

## Sites with more observations suggest mild autocorrelation 
# Maurice How do I interpret this? These numbers have all changed
arima(resid_ts$`120 Reef`, c(1, 0, 0)) # 0.2
arima(resid_ts$`Casino Point`, c(1, 0, 0)) # 0.7
arima(resid_ts$`Malaga Cove`, c(1, 0, 0)) # 0.48
arima(resid_ts$`La Jolla Cove`, c(1, 0, 0)) # 0.52
arima(resid_ts$`Christmas Tree Cove`, c(1, 0, 0)) # 0.1

##### Model effects plots #####################################################
cowplot::plot_grid(
  plot(emmeans::emmeans(model_mpa1, "mpa_status"), horizontal = FALSE), 
  plot(emmeans::emmeans(model_mpa1, "heatwave"), horizontal = FALSE)
)

interaction_data <- as.data.frame(emmeans::emmeans(model_mpa1, "heatwave", by = "mpa_status", type = "response"))

interaction_data %>% 
  ggplot(aes(heatwave, color = mpa_status)) + 
  geom_pointrange(aes(y = response, ymin = lower.CL, ymax = upper.CL), position = position_dodge(width = 0.2), linewidth = 1)

##### Main effects tests ######################################################
car::Anova(model_mpa1)
car::Anova(model_tc1)

summary(model_mpa1)
summary(model_tc1)
summary(model_tc2)

##### Create Publication Figure ################################################
em1 <- ref_grid(model_tc1, at = list(avg_sheephead = seq(0, max(data$avg_sheephead), length.out = 100)))
em2 <- ref_grid(model_tc2, at = list(avg_sheephead = seq(0, max(data$avg_sheephead), length.out = 100)))

sheephead1 <- as.data.frame(emmip(em1, ~avg_sheephead, type = "response", CIs = TRUE, plot = FALSE))
sheephead2 <- as.data.frame(emmip(em2, ~avg_sheephead, type = "response", CIs = TRUE, plot = FALSE))

em1 <- ref_grid(model_tc1, at = list(avg_lobster = seq(0, max(data$avg_lobster), length.out = 100)))
em2 <- ref_grid(model_tc2, at = list(avg_lobster = seq(0, max(data$avg_lobster), length.out = 100)))

lobster1 <- as.data.frame(emmip(em1, ~avg_lobster, type = "response", CIs = TRUE, plot = FALSE))
lobster2 <- as.data.frame(emmip(em2, ~avg_lobster, type = "response", CIs = TRUE, plot = FALSE))


plotdata <- do.call("rbind", list(
  rename(mutate(sheephead1, model = "random intercepts", species = "sheephead"), x = avg_sheephead),
  rename(mutate(sheephead2, model = "random intercepts + AR(1)", species = "sheephead"), x = avg_sheephead),
  rename(mutate(lobster1, model = "random intercepts", species = "lobster"), x = avg_lobster),
  rename(mutate(lobster2, model = "random intercepts + AR(1)", species = "lobster"), x = avg_lobster)
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
    axis.title.x = element_blank()
  ) +
  ylab("Mean urchin density")
  

dat_text <- data.frame(
  label = c("p = 0.0001", "p = 0.0212", "p = 0.1502", "p = 0.9081"),
  species = c("sheephead", "sheephead", "lobster", "lobster"),
  model = c("random intercepts", "random intercepts + AR(1)", "random intercepts", "random intercepts + AR(1)")
)

plot_final <- plot1 + geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.1,
  vjust   = -1
)

##### Export Plots #############################################################
png("Figures/GLMM_model_results.png", width = 6, height = 6, 
    units = "in", res = 600)
plot_final
dev.off() 

# Violin plot
png("Figures/Autocorrelation.png", width = 6, height = 6, 
    units = "in", res = 600)
violin_plot
dev.off() 

