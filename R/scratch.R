# scratch.R
#
# document for code that doesn't fit anywhere. For example exploring options for how to analyse data, that doesn't become the final analysis
#
# Created: 13/11/26
# Created by: Davina Gifford


# explorting the new variables

ggpairs(polyp_abund)

m1 <- lm(relative_abundance ~ salinity_level, data = polyp_abund)
summary(m1)
intercept <- coef(m1)[1]  # Intercept (b0)
slope <- coef(m1)[2]      # Slope (b1)

# Plot with custom regression line
ggplot(data = polyp_abund, aes(x = salinity_level, y = relative_abundance)) +
  geom_point() + 
  geom_abline(intercept = intercept, slope = slope) +
  labs(x = "Salinity (ppt)", y = "relative abundance of polyps")

plot(m1)
gglm(m1)

m2 <- lm(log1p(relative_abundance) ~ salinity_level, data = polyp_abund) # A log10 transformation of relative abundance

summary(m2)
plot(m2)
gglm(m2)

polyp_abund <- polyp_abund %>% 
  mutate(across(tidyselect::where(is.character),as.factor))

str(polyp_abund)

m3 <- lm(relative_abundance ~ salinity_level, data = polyp_abund)
summary(m3)

log_abund <- polyp_abund$relative_abundance
log_abund <- log1p(log_abund)

sqrt_abund <- cbrt(polyp_abund$relative_abundance)

hist(sqrt_abund)


library(dplyr)
library(tidyr)
library(ggplot2)

transformed_data <- polyp_abund %>%
  mutate(
    raw = relative_abundance,
    sqrt = sqrt(relative_abundance),
    fourth_root = relative_abundance^(1/4),
    log1p = log1p(relative_abundance)
  ) %>%
  select(raw, sqrt, fourth_root, log1p) %>%
  pivot_longer(
    everything(),
    names_to = "transformation",
    values_to = "value"
  )

ggplot(transformed_data, aes(x = value)) +
  geom_histogram(bins = 20, fill = "grey70", colour = "black") +
  facet_wrap(~ transformation, scales = "free_x") +
  theme_classic()


str(polyp_abund)



