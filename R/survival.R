# survival.R
#
# Calculating survival rates of polyps. 
# calculate relative abundance of polyps and proportion of budding
# 
# Created: 13/11/2026
# Created by: Davina Gifford

# load libraries
library(tidyverse)
library(GGally)
library(gglm)
library(glmmTMB)
library(marginaleffects)
library(DHARMa)



# load data
polyp_data <- readRDS("var/polyp_tidy.rds")
treatment_data <- readRDS("var/treatment_tidy.rds")
medusa_data <- readRDS("var/medusae_tidy.rds")

treatment_data <- treatment_data %>%
  rename(salinity_level = Treatment)

# combine into a single var for analysis

plot_data <- polyp_data %>%
  left_join(
    treatment_data,
    by = c(
      "salinity_level" = "salinity_level",
      "Day" = "Day"
    )
  )

# check for missing data
polyp_abund %>%
  group_by(salinity_level, replicate) %>%
  summarise(
    has_day1 = any(Day == 15),
    initial_polyps = first(initial_polyps),
    .groups = "drop"
  ) %>%
  filter(!has_day1)

# check for duplicates
polyp_abund %>%
  count(Day, salinity_level, replicate) %>%
  filter(n > 1)

# create survival data

survival_polyps <- plot_data %>%
  group_by(salinity_level, replicate) %>%
  mutate(
    initial_polyps = polyp_number[Day == 1][1],
    survival = polyp_number / initial_polyps
  ) %>%
  ungroup()

# net population change

survival_polyps <- plot_data %>%
  group_by(salinity_level, replicate) %>%
  #arrange(Day, .by_group = TRUE) %>%
  mutate(
    dead = replace_na(dead, 0),
    initial_polyps = polyp_number[Day == 1][1],
    previous_polyps = lag(polyp_number), # polyps present at previous observation
    pop_abundance = polyp_number - previous_polyps, # population change relative to day 1
    pop_change = if_else(
      initial_polyps > 0,
      polyp_number / initial_polyps,
      NA_real_
    ),
    recruitment = pop_abundance + dead, #  inferred recruitment
    unexpected_loss = pmax(   # unexpected losses
      -(pop_abundance + dead), 0
    ),
    mortality_rate = if_else(
      previous_polyps > 0,
      dead / previous_polyps,
      NA_real_
    ), # mortality rate relative to previous observation
    
  ) %>%
  ungroup()



# different approach due to limitations of data collection

polyp_abund <- plot_data %>%
  group_by(salinity_level, replicate) %>%
  #arrange(Day, .by_group = TRUE) %>%
  mutate(
    dead = replace_na(dead, 0),
    
    # Initial population
    initial_polyps = polyp_number[Day == 1][1],
    
    # Relative abundance since Day 1
    relative_polyp_abundance = if_else(
      initial_polyps > 0,
      polyp_number / initial_polyps,
      NA_real_
    ),
    
    # Previous observation
    previous_polyps = lag(polyp_number),
    days_elapsed = Day - lag(Day),
    
    # Net population change
    population_delta = polyp_number - previous_polyps,
    
    # Daily net growth rate
    growth_rate = if_else(
      previous_polyps > 0 &
        polyp_number > 0 &
        days_elapsed > 0,
      log(polyp_number / previous_polyps) / days_elapsed,
      NA_real_
    ),
    
    # Observed mortality (minimum recorded deaths)
    observed_deaths = dead
  ) %>%
  # budding ratio
  mutate(
    budding_ratio = if_else(
      polyp_number > 0,
      replace_na(buds, 0) / polyp_number,
      NA_real_
    )
  ) %>% 
  ungroup()

# set factors
polyp_abund <- polyp_abund %>% 
  mutate(across(where(is.character), as.factor))

polyp_abund <- polyp_abund %>% 
  mutate(plot_day = Day - 2)


# create unique ids for each replicate
polyp_abund <- polyp_abund %>%
  mutate(
    replicate_id = interaction(salinity_level, replicate)
  )



# Plotting budding ratio --------------------------------------------------


# plot the budding ratio for the individual replicates 

replicates <- polyp_abund %>%
  distinct(salinity_level, replicate)

for(i in seq_len(nrow(replicates))) {
  
  treatment_i <- replicates$salinity_level[i]
  replicate_i <- replicates$replicate[i]
  
  dat <- polyp_abund %>%
    filter(
      salinity_level == treatment_i,
      replicate == replicate_i
    )
  
  
  scale_factor <- max(dat$budding_ratio, na.rm = TRUE) / 40
  
  print(
    dat %>%
      select(plot_day, salinity_level, replicate, budding_ratio)
  )
  
  p <- ggplot(dat, aes(x = plot_day)) +
    
    geom_line(
      aes(y = budding_ratio),
      colour = "#0072B2",
      linewidth = 1
    ) +
    
    geom_point(
      aes(y = budding_ratio),
      colour = "#0072B2",
      size = 2
    ) +
    
    geom_step(
      aes(y = Salinity * scale_factor),
      colour = "#D55E00",
      linewidth = 1
    ) +
    
    scale_y_continuous(
      name = "Ratio of polyps with buds",
      
      sec.axis = sec_axis(
        ~ . / scale_factor,
        name = "Salinity (ppt)"
      )
    ) +
    
    scale_x_continuous(
      breaks = -1:15
    ) +
    
    theme_bw() +
    
    labs(
      title = paste(
        treatment_i,
        "- Replicate",
        replicate_i
      ),
      x = "Experiment day"
    )
  
  ggsave(
    filename = paste0(
      "output/bud_ratio/",
      treatment_i,
      "_bud_ratio_rep",
      replicate_i,
      ".png"
    ),
    plot = p,
    width = 7,
    height = 4.5,
    dpi = 300
  )
}

# Plot by treatment


for(treatment_i in unique(polyp_abund$salinity_level)) {
  
  dat <- polyp_abund %>%
    filter(salinity_level == treatment_i) %>%
    arrange(replicate, plot_day)
  
  scale_factor <- max(dat$budding_ratio, na.rm = TRUE) / 40
  
  phase_data <- tibble(
    xmin = c(-0.5, 3.5),
    xmax = c(3.5, 8.5),
    Phase = c("Exposure", "Recovery")
  )
  
  p <- ggplot(dat, aes(x = plot_day)) +
    
    # Treatment phases
    geom_rect(
      data = phase_data,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = Phase
      ),
      inherit.aes = FALSE,
      alpha = 0.1
    ) +
    
    # Replicate trajectories
    geom_line(
      aes(
        y = budding_ratio,
        colour = factor(replicate),
        group = replicate
      ),
      linewidth = 1
    ) +
    
    geom_point(
      aes(
        y = budding_ratio,
        colour = factor(replicate)
      ),
      size = 2
    ) +
    
    # Salinity profile
    geom_step(
      aes(
        y = Salinity * scale_factor,
        linetype = "Salinity"
      ),
      colour = "black",
      linewidth = 1.2
    ) +
    
    scale_fill_manual(
      name = "Experimental phase",
      values = c(
        Exposure = "firebrick",
        Recovery = "steelblue"
      )
    ) +
    
    scale_linetype_manual(
      name = "",
      values = c(
        Salinity = "solid"
      )
    ) +
    
    scale_colour_brewer(
      palette = "Dark2",
      name = "Replicate"
    ) +
    
    scale_y_continuous(
      name = "Ratio of polyps with buds",
      sec.axis = sec_axis(
        ~ . / scale_factor,
        name = "Salinity (ppt)"
      )
    ) +
    
    theme_bw()
  
  
  
  ggsave(
    filename = paste0(
      "output/bud_ratio/",
      treatment_i,
      "_bud_ratio_polyps.png"
    ),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# facet wrap


scale_factor <- max(polyp_abund$budding_ratio, na.rm = TRUE) / 40

# mean polyp values per treatment
mean_data_budding <- polyp_abund %>%
  group_by(salinity_level, plot_day) %>%
  summarise(
    mean_budding = mean(budding_ratio, na.rm = TRUE),
    sd = sd(budding_ratio, na.rm = TRUE),
    .groups = "drop"
  )




p <- ggplot(
  polyp_abund,
  aes(
    x = plot_day,
    y = budding_ratio,
    colour = factor(replicate),
    group = replicate
  )
) +
  
  # Exposure period
  annotate(
    "rect",
    xmin = -0.5,
    xmax = 3.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "firebrick",
    alpha = 0.08
  ) +
  
  # Recovery period
  annotate(
    "rect",
    xmin = 3.5,
    xmax = 8.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "steelblue",
    alpha = 0.08
  ) +
  
  # SE ribbon
  geom_ribbon(
    data = mean_data_budding,
    aes(
      x = plot_day,
      ymin = mean_budding - sd,
      ymax = mean_budding + sd,
      group = 1
    ),
    inherit.aes = FALSE,
    fill = "grey40",
    alpha = 0.2
  ) +
  
  # Individual replicate trajectories
  geom_line(
    linewidth = 1,
    alpha = 0.6
  ) +
  
  geom_point(
    size = 2,
    alpha = 0.6
  ) +
  
  # Treatment mean
  geom_line(
    data = mean_data_budding,
    aes(
      x = plot_day,
      y = mean_budding,
      group = 1
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.5
  ) +
  
  geom_point(
    data = mean_data_budding,
    aes(
      x = plot_day,
      y = mean_budding
    ),
    inherit.aes = FALSE,
    colour = "black",
    size = 2.5
  ) +
  
  # Phase labels
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    label = "Exposure",
    vjust = 1.5,
    fontface = "bold",
    colour = "firebrick"
  ) +
  
  annotate(
    "text",
    x = 6,
    y = Inf,
    label = "Recovery",
    vjust = 1.5,
    fontface = "bold",
    colour = "steelblue"
  ) +
  
  annotate(
    "text",
    x = 11,
    y = Inf,
    label = "Post-recovery",
    vjust = 1.5,
    fontface = "bold"
  ) +
  
  facet_wrap(
    ~ salinity_level,
    scales = "fixed"
  ) +
  
  scale_x_continuous(
    breaks = -1:15
  ) +
  
  scale_colour_brewer(
    palette = "Dark2",
    name = "Replicate"
  ) +
  
  labs(
    x = "Experiment Day",
    y = "Ratio of polyps budding"
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) 

ggsave(
  filename = 
    "output/budding_ratio_results_by_treatment.png",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)



# plot relative abundance -------------------------------------------------

# plot the relative_abundances for the individual replicates 

replicates <- polyp_abund %>%
  distinct(salinity_level, replicate)

for(i in seq_len(nrow(replicates))) {
  
  treatment_i <- replicates$salinity_level[i]
  replicate_i <- replicates$replicate[i]
  
  dat <- polyp_abund %>%
    filter(
      salinity_level == treatment_i,
      replicate == replicate_i
    )
  
  
  scale_factor <- max(dat$relative_polyp_abundance, na.rm = TRUE) / 40
  
  print(
    dat %>%
      select(plot_day, salinity_level, replicate, relative_polyp_abundance)
  )
  
  p <- ggplot(dat, aes(x = plot_day)) +
    
    geom_line(
      aes(y = relative_polyp_abundance),
      colour = "#0072B2",
      linewidth = 1
    ) +
    
    geom_point(
      aes(y = relative_polyp_abundance),
      colour = "#0072B2",
      size = 2
    ) +
    
    geom_step(
      aes(y = Salinity * scale_factor),
      colour = "#D55E00",
      linewidth = 1
    ) +
    
    scale_y_continuous(
      name = "Relative abundance of polyps",
      
      sec.axis = sec_axis(
        ~ . / scale_factor,
        name = "Salinity (ppt)"
      )
    ) +
    
    scale_x_continuous(
      breaks = -1:15
    ) +
    
    theme_bw() +
    
    labs(
      title = paste(
        treatment_i,
        "- Replicate",
        replicate_i
      ),
      x = "Experiment day"
    )
  
  ggsave(
    filename = paste0(
      "output/rel_abund/",
      treatment_i,
      "_rel_abund_rep",
      replicate_i,
      ".png"
    ),
    plot = p,
    width = 7,
    height = 4.5,
    dpi = 300
  )
}

# Plot by treatment


for(treatment_i in unique(polyp_abund$salinity_level)) {
  
  dat <- polyp_abund %>%
    filter(salinity_level == treatment_i) %>%
    arrange(replicate, plot_day)
  
  scale_factor <- max(dat$relative_polyp_abundance, na.rm = TRUE) / 40
  
  phase_data <- tibble(
    xmin = c(0.5, 3.5),
    xmax = c(3.5, 8.5),
    Phase = c("Exposure", "Recovery")
  )
  
  p <- ggplot(dat, aes(x = plot_day)) +
    
    # Treatment phases
    geom_rect(
      data = phase_data,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = Phase
      ),
      inherit.aes = FALSE,
      alpha = 0.1
    ) +
    
    # Replicate trajectories
    geom_line(
      aes(
        y = relative_polyp_abundance,
        colour = factor(replicate),
        group = replicate
      ),
      linewidth = 1
    ) +
    
    geom_point(
      aes(
        y = relative_polyp_abundance,
        colour = factor(replicate)
      ),
      size = 2
    ) +
    
    # Salinity profile
    geom_step(
      aes(
        y = Salinity * scale_factor,
        linetype = "Salinity"
      ),
      colour = "black",
      linewidth = 1.2
    ) +
    
    scale_fill_manual(
      name = "Experimental phase",
      values = c(
        Exposure = "firebrick",
        Recovery = "steelblue"
      )
    ) +
    
    scale_linetype_manual(
      name = "",
      values = c(
        Salinity = "solid"
      )
    ) +
    
    scale_colour_brewer(
      palette = "Dark2",
      name = "Replicate"
    ) +
    
    scale_y_continuous(
      name = "Relative polyp abundance",
      sec.axis = sec_axis(
        ~ . / scale_factor,
        name = "Salinity (ppt)"
      )
    ) +
    
    theme_bw()
  
  
  
  ggsave(
    filename = paste0(
      "output/rel_abund/",
      treatment_i,
      "_rel_abund_polyps.png"
    ),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# facet wrap

treatment_data <- treatment_data %>%
  rename(salinity_level = Treatment)

scale_factor <- max(polyp_abund$relative_polyp_abundance, na.rm = TRUE) / 40

# mean polyp values per treatment
mean_data_abund <- polyp_abund %>%
  group_by(salinity_level, plot_day) %>%
  summarise(
    mean_abund = mean(relative_polyp_abundance, na.rm = TRUE),
    sd = sd(relative_polyp_abundance, na.rm = TRUE),
    .groups = "drop"
  )






p <- ggplot(
  polyp_abund,
  aes(
    x = plot_day,
    y = relative_polyp_abundance,
    colour = factor(replicate),
    group = replicate
  )
) +
  
  # Exposure period
  annotate(
    "rect",
    xmin = -0.5,
    xmax = 3.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "firebrick",
    alpha = 0.08
  ) +
  
  # Recovery period
  annotate(
    "rect",
    xmin = 3.5,
    xmax = 8.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "steelblue",
    alpha = 0.08
  ) +
  
  # SE ribbon
  geom_ribbon(
    data = mean_data_abund,
    aes(
      x = plot_day,
      ymin = mean_abund - sd,
      ymax = mean_abund + sd,
      group = 1
    ),
    inherit.aes = FALSE,
    fill = "grey40",
    alpha = 0.2
  ) +
  
  # Individual replicate trajectories
  geom_line(
    linewidth = 1,
    alpha = 0.6
  ) +
  
  geom_point(
    size = 2,
    alpha = 0.6
  ) +
  
  # Treatment mean
  geom_line(
    data = mean_data_abund,
    aes(
      x = plot_day,
      y = mean_abund,
      group = 1
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.5
  ) +
  
  geom_point(
    data = mean_data_abund,
    aes(
      x = plot_day,
      y = mean_abund
    ),
    inherit.aes = FALSE,
    colour = "black",
    size = 2.5
  ) +
  
  # Phase labels
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    label = "Exposure",
    vjust = 1.5,
    fontface = "bold",
    colour = "firebrick"
  ) +
  
  annotate(
    "text",
    x = 6,
    y = Inf,
    label = "Recovery",
    vjust = 1.5,
    fontface = "bold",
    colour = "steelblue"
  ) +
  
  annotate(
    "text",
    x = 11,
    y = Inf,
    label = "Post-recovery",
    vjust = 1.5,
    fontface = "bold"
  ) +
  
  facet_wrap(
    ~ salinity_level,
    scales = "fixed"
  ) +
  
  scale_x_continuous(
    breaks = -1:15
  ) +
  
  scale_colour_brewer(
    palette = "Dark2",
    name = "Replicate"
  ) +
  
  labs(
    x = "Experiment Day",
    y = "Relative polyp abundance"
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) 

ggsave(
  filename = 
    "output/polyp_rel_abund_results_by_treatment.png",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)



# analysis ----------------------------------------------------------------



# exploring GLMMs
m1 <- glm(relative_polyp_abundance ~ salinity_level, data = polyp_abund)
summary(m1)
gglm(m1)


m2 <- lm(relative_abundance ~ salinity_level + phase + salinity_level:phase, data = polyp_abund)
summary(m2)

m3 <- lm(relative_abundance ~ salinity_level + phase, data = polyp_abund)
summary(m3)


m4 <- lm(relative_abundance ~ salinity_level + salinity_level:phase, data = polyp_abund)
summary(m4)

AIC(m1, m2, m3, m4)
AIC(m1, m2, m3, m4, k=log(nrow(dat)))

m5 <- glmmTMB(relative_abundance ~ salinity_level + (1|replicate_id), data = polyp_abund)
summary(m5)




m6 <- glmmTMB(relative_abundance ~ salinity_level + (1 + salinity_level|replicate_id), data = polyp_abund, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(m6)
diagnose(m7)

AIC(m5, m6)
anova(m5,m6)

simulationOutput <- simulateResiduals(fittedModel = m6, plot = FALSE)
plotQQunif(simulationOutput, testUniformity = FALSE, testOutliers = FALSE, testDispersion = FALSE)
plotResiduals(simulationOutput, quantreg = FALSE)



m7 <- glmmTMB(polyp_number ~ salinity_level * Day + (1 | replicate_id), family = nbinom2(), data = polyp_abund)
summary(m7)

m8 <- glmmTMB(relative_abundance ~ salinity_level * Day + (1 | replicate_id), family = gaussian(), data = polyp_abund)
summary(m8)

anova(m5,m8)

simulationOutput <- simulateResiduals(fittedModel = m8, plot = FALSE)
plotQQunif(simulationOutput, testUniformity = FALSE, testOutliers = FALSE, testDispersion = FALSE)

m9 <- glmmTMB(
  polyp_number ~ salinity_level * Day +
    offset(log(initial_polyps)) +
    (1 | replicate_id),
  family = nbinom2(),
  data = polyp_abund
)
summary(m9)

anova(m9,m8)




# exploring anovas
# want to do Anovas at 3 points in time. 1) end of exposure, 2) end of recovery, 3) end of experiment

exposure_data <- polyp_abund %>% 
  filter(Day == 6)

summary(exposure_data)
ggplot(exposure_data) +
  aes(x = salinity_level, y = relative_polyp_abundance, color = salinity_level) +
  geom_jitter() +
  theme(legend.position = "none")

m10 <- lm(relative_polyp_abundance ~ salinity_level * replicate_id, data = exposure_data)
anova(m10)


# anova

aov1 <- aov(relative_polyp_abundance ~ salinity_level, data = exposure_data)

hist(aov1$residuals)
library(car)
qqPlot(aov1$residuals,
       id = FALSE # id = FALSE to remove point identification
)
shapiro.test(aov1$residuals)

boxplot(relative_polyp_abundance ~ salinity_level, data = exposure_data)

summary(aov1)

kruskal.test(relative_polyp_abundance ~ salinity_level,
             data = exposure_data)



# transform data 
polyp_abund <- polyp_abund %>%
  mutate(
    log_relative_abundance = logit(relative_polyp_abundance)
  )

# extract data for the 3 days for each phase
key_times <- polyp_abund %>%
  filter(Day %in% c(6, 8, 15))


library(rstatix)

key_times %>%
  group_by(Day) %>%
  anova_test(
    log_relative_abundance ~ salinity_level
  )



exposure_data %>%
  select(salinity_level,
         replicate,
         polyp_number,
         relative_polyp_abundance)

summary(polyp_abund$relative_polyp_abundance)

hist(polyp_abund$relative_polyp_abundance)

summary(exposure_data$relative_polyp_abundance)
ggplot(exposure_data,
       aes(sample = relative_polyp_abundance)) +
  stat_qq() +
  stat_qq_line()

# filter data for each time point

day6 <- polyp_abund %>%
  filter(Day == 6)

day8 <- polyp_abund %>%
  filter(Day == 8)

day15 <- polyp_abund %>%
  filter(Day == 15)

 
day6_aov <- 
  aov(relative_polyp_abundance ~ salinity_level, data = day6)
summary(day6_aov)

ggplot(day15,
       aes(sample = relative_polyp_abundance)) +
  stat_qq() +
  stat_qq_line()

boxplot(relative_polyp_abundance ~ salinity_level, data = day6)

shapiro.test(day6_aov$residuals)
qqPlot(day6_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)

day8_aov <- 
  aov(relative_polyp_abundance ~ salinity_level, data = day8)
summary(day8_aov)

boxplot(relative_polyp_abundance ~ salinity_level, data = day8)
qqPlot(day8_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)

day15_aov <- 
  aov(relative_polyp_abundance ~ salinity_level, data = day15)
summary(day15_aov)

boxplot(relative_polyp_abundance ~ salinity_level, data = day15)


qqPlot(day15_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)


library(xlsx)
library(rstatix)
library(reshape)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(plyr)
library(datarium)


# repeated measures anova at the three time points

# prepare the data

rep_measures_data <- day6 %>% 
  bind_rows(day8)

rep_measures_data <- rep_measures_data %>% 
  bind_rows(day15)

rep_measures_data <- rep_measures_data %>% 
  select(replicate_id, phase, relative_polyp_abundance)

# boxplot of each time point
ggplot(rep_measures_data) +
  aes(x = phase, y = relative_polyp_abundance) +
  geom_boxplot() +
  labs(
    x = "Phase",
    y = "Relative polyp abundance"
  )

# spaghetti plot

ggplot(rep_measures_data) +
  aes(x = phase, y = relative_polyp_abundance, group = replicate_id) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  stat_summary(aes(group = 1),
               fun = mean, geom = "line",
               linewidth = 1.2, color = "steelblue"
  ) +
  stat_summary(aes(group = 1),
               fun = mean, geom = "point",
               size = 3, color = "steelblue"
  ) +
  labs(
    x = "Experiment phase",
    y = "Relative polyp abundance"
  )

res_lm <- lm(relative_polyp_abundance ~ phase + replicate_id, data = rep_measures_data)

par(mfrow = c(1, 2)) # combine plots
hist(residuals(res_lm),
     main = "Histogram of the residuals",
     xlab = "Residuals"
)
qqPlot(residuals(res_lm),
       id = FALSE # remove point identification
)
shapiro.test(residuals(res_lm))

# normality is not being met, even with a log transform
# repeated measures anova might not be best option


# sphericity
# from the long format to the wide format
dat_wide <- rep_measures_data %>%
  pivot_wider(names_from = phase, values_from = relative_polyp_abundance)

c(
  "recovery - exposure" = var(dat_wide$Recovery - dat_wide$Exposure),
  "post-recovery - exposure" = var(dat_wide$`Post-Recovery` - dat_wide$Exposure),
  "post-recovery - recovery" = var(dat_wide$`Post-Recovery` - dat_wide$Recovery)
)


res_aov <- anova_test(
  data = rep_measures_data,
  dv = relative_polyp_abundance,
  wid = replicate_id,
  within = phase
)

res_aov

get_anova_table(res_aov, correction = "GG")


summary(aov(relative_polyp_abundance ~ phase + Error(patient / time),
            data = dat
))


# repeated measures anova for salinity level

# prepare the data

sal_data_1 <- polyp_abund %>% 
  filter(Day <= 1)

sal_data_end <- polyp_abund %>% 
  filter(Day == 15)

rep_measures_data_salinity <- sal_data_1 %>% 
  bind_rows(sal_data_end)


rep_measures_data_salinity <- rep_measures_data_salinity %>% 
  select(replicate_id, salinity_level, relative_polyp_abundance)

# boxplot of each time point
ggplot(rep_measures_data_salinity) +
  aes(x = salinity_level, y = relative_polyp_abundance) +
  geom_boxplot() +
  labs(
    x = "salinity level",
    y = "Relative polyp abundance"
  )

# spaghetti plot

ggplot(rep_measures_data_salinity) +
  aes(x = salinity_level, y = relative_polyp_abundance, group = replicate_id) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  stat_summary(aes(group = 1),
               fun = mean, geom = "line",
               linewidth = 1.2, color = "steelblue"
  ) +
  stat_summary(aes(group = 1),
               fun = mean, geom = "point",
               size = 3, color = "steelblue"
  ) +
  labs(
    x = "Salinity level",
    y = "Relative polyp abundance"
  )

res_lm_sal <- lm(relative_polyp_abundance ~ salinity_level + replicate_id, data = rep_measures_data_salinity)

par(mfrow = c(1, 2)) # combine plots
hist(residuals(res_lm_sal),
     main = "Histogram of the residuals",
     xlab = "Residuals"
)
qqPlot(residuals(res_lm_sal),
       id = FALSE # remove point identification
)
shapiro.test(residuals(res_lm_sal))

# normality is not being met, even with a log transform
# repeated measures anova might not be best option


# sphericity
# from the long format to the wide format
dat_wide <- rep_measures_data_salinity %>%
  pivot_wider(names_from = salinity_level, values_from = relative_polyp_abundance)

c(
  "recovery - exposure" = var(dat_wide$Recovery - dat_wide$Exposure),
  "post-recovery - exposure" = var(dat_wide$`Post-Recovery` - dat_wide$Exposure),
  "post-recovery - recovery" = var(dat_wide$`Post-Recovery` - dat_wide$Recovery)
)


res_aov <- anova_test(
  data = rep_measures_data,
  dv = relative_polyp_abundance,
  wid = replicate_id,
  within = phase
)

res_aov

get_anova_table(res_aov, correction = "GG")


summary(aov(relative_polyp_abundance ~ phase + Error(patient / time),
            data = dat
))



# another attempt

library(dplyr)
library(ggplot2)
library(afex)
library(emmeans)

new_rep_measures <- polyp_abund %>% 
  select(salinity_level, Day, phase, replicate_id, polyp_number, relative_polyp_abundance)

model <- aov_ez(
  id = "replicate_id",
  dv = "polyp_number",
  data = new_rep_measures,
  between = "salinity_level",
  within = "phase",
  type = 3
)

summary(model)
nice(model)


str(rep_measures_data_salinity)


phase_df <- new_rep_measures %>%
  group_by(replicate_id, salinity_level, phase) %>%
  summarise(
    mean_polyps = mean(polyp_number, na.rm = TRUE),
    .groups = "drop"
  )

phase_model <- aov_ez(
  id = "replicate_id",
  dv = "mean_polyps",
  data = phase_df,
  between = "salinity_level",
  within = "phase",
  type = 3
)

summary(phase_model)
nice(phase_model)


# 19 November ---------------------------------
# restart! 
library(tidyverse)
library(glmmTMB)
library(car)      # For overall ANOVA p-values
library(emmeans)  # For pairwise comparisons

# remove day 1 from data to be analysed
dat_new <- polyp_abund %>% 
  filter(phase != "Acclimation")



model_new <- glmmTMB(
  polyp_number ~ salinity_level * phase + offset(log(initial_polyps)) + (1 | replicate_id),
  family = nbinom2,
  ziformula = ~1,
  data = dat_new
)

Anova(model_new, type = "III")

library(DHARMa)

# Simulate residuals and check for overdispersion, zero-inflation, and outliers
sim_res <- simulateResiduals(model_new, plot = TRUE)
testZeroInflation(sim_res)

# Option A: Allow zero-inflation to vary by salinity and phase (additive)
model_zi_add <- glmmTMB(
  polyp_number ~ salinity_level * phase + offset(log(initial_polyps)) + (1 | replicate_id),
  family = nbinom2,
  ziformula = ~ salinity_level,
  data = dat_new
)

# Option B: Allow zero-inflation to vary by the full interaction (if Option A still shows ZI)
model_zi_int <- glmmTMB(
  polyp_number ~ salinity_level * phase + offset(log(initial_polyps)) + (1 | replicate_id),
  family = nbinom2,
  ziformula = ~ salinity_level * phase,
  data = dat_new
)

# 1. Zero-inflation driven ONLY by salinity
model_zi_sal <- glmmTMB(
  polyp_number ~ salinity_level * phase + offset(log(initial_polyps)) + (1 | replicate_id),
  family = nbinom2,
  ziformula = ~ salinity_level,
  data = dat_new
)

# 2. Zero-inflation driven ONLY by phase
model_zi_phase <- glmmTMB(
  polyp_number ~ salinity_level * phase + offset(log(initial_polyps)) + (1 | replicate_id),
  family = nbinom2,
  ziformula = ~ phase,
  data = dat_new
)

# 3. Tweedie distribution (handles excess zeros natively without ziformula)
model_tweedie <- glmmTMB(
  polyp_number ~ salinity_level * phase + offset(log(initial_polyps)) + (1 | replicate_id),
  family = tweedie,
  data = dat_new
)

# Compare AIC across all converged models
AIC(model_new, model_zi_sal, model_zi_phase, model_tweedie)

summary(model_zi_sal)



# try something new
fullModel <- with(dat_new, lm(relative_abundance ~ salinity_level * phase * initial_polyps))
qqnorm(fullModel$residuals)                                         
# View qq plot to check assumption 1
qqline(fullModel$residuals) 
plot(fullModel$fitted.values, fullModel$residuals)  

logFullModel <- with(dat_new, lm(log1p(relative_abundance) ~ salinity_level * phase *  replicate_id))
sqrtFullModel <- with(dat_new, lm(sqrt(relative_abundance) ~ salinity_level * phase * initial_polyps))

par(mfrow = c(1, 3))  # set the plot area so that we can plot three plots side by side (1 row, 3 columns)

## first panel in the plot
qqnorm(fullModel$residuals, main = "Untransformed Y")                                        
qqline(fullModel$residuals) 

## second panel in the plot
qqnorm(sqrtFullModel$residuals, main = "Square Root Transformed")                           
qqline(sqrtFullModel$residuals) 
## third panel in the plot
qqnorm(logFullModel$residuals, main = "Log Transformed")
qqline(logFullModel$residuals) 

par(mfrow = c(3, 1))  # Graphic parameters: 3 rows, 1 column
#having all plots in one panel can help to compare, but if they get too squished to see, revert to par(mfrow = c(1, 1)), or try another combination of rows and columns that does work for you. 
plot(fullModel$fitted.values, fullModel$residuals, main="Untransformed Y")  
plot(sqrtFullModel$fitted.values, sqrtFullModel$residuals, main="Square Root Transformed")
plot(logFullModel$fitted.values, logFullModel$residuals, main="Log Transformed")
par(mfrow = c(1, 1)) # this resets the plot area to show one plot at a time

dat_new.subset <- dat_new %>% 
  select(salinity_level, phase, initial_polyps)

cor(dat_new.subset)# use the `cor` function from base R to generate the matrix of correlation coefficients
round(cor(dat_new.subset), 3)# use the `round` function to make it easier to inspect by truncating the number of decimal places to 3

library(stats)# Load the library that contains the "step" program
stepModel <- step(logFullModel, direction = "both")

anova(stepModel)# Look at the significance tests for the best model from `step`
qqnorm(stepModel$residuals, main = "Step Model")# Check the diagnostics
qqline(stepModel$residuals)
plot(stepModel$fitted.values, stepModel$residuals, main = "Step Model")

Anova(logFullModel, Type = "II")
summary(logFullModel)



# another approach



library(glmmTMB)

model_3 <- glmmTMB(
  polyp_number ~ salinity_level * Day +
    initial_polyps +
    (1 | replicate_id),
  family = nbinom2,
  data = dat_new
)

summary(model_3)

model_4 <- glmmTMB(
  polyp_number ~ salinity_level * Day +
    offset(log(initial_polyps)) +
    (1 | replicate_id),
  family = nbinom2,
  data = dat_new
)

summary(model_4)


summary(dat_new)


# OTAY. lets try again.
# negative binomial, piecewise model 

# Load packages
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(performance)

polyp_abund <- polyp_abund %>% 
  mutate(
    # piecewise time variables
    t_exp = pmin(pmax(Day - 1, 0), 4),
    t_rec = pmin(pmax(Day - 5, 0), 5),
    t_post = pmax(Day - 10, 0)
  )

# Inspect the time variables
polyp_abund %>%
  distinct(Day, phase, t_exp, t_rec, t_post) %>%
  arrange(Day) %>%
  print(n = 20)

# visualise the data

ggplot(
  polyp_abund,
  aes(
    x = Day,
    y = polyp_number,
    group = replicate_id
  )
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1) +
  facet_wrap(~ salinity_level) +
  geom_vline(
    xintercept = c(5, 10),
    linetype = "dashed",
    colour = "grey50"
  ) +
  labs(
    x = "Experimental day",
    y = "Number of polyps",
    title = "Polyp population trajectories"
  ) +
  theme_classic()



polyp_abund %>%
  group_by(salinity_level) %>%
  summarise(
    Mean = mean(polyp_number, na.rm = TRUE),
    Variance = var(polyp_number, na.rm = TRUE),
    Zeros = sum(polyp_number == 0, na.rm = TRUE),
    Observations = n(),
    .groups = "drop"
  )


ggplot(polyp_abund, aes(x = polyp_number)) +
  geom_histogram(
    binwidth = 2,
    boundary = 0,
    colour = "white"
  ) +
  facet_wrap(~ salinity_level) +
  labs(
    x = "Polyp count",
    y = "Frequency",
    title = "Distribution of polyp abundance"
  ) +
  theme_classic()


# Exclude baseline observations
model_data <- polyp_abund %>%
  filter(Day > 1) %>%
  droplevels()


model_data2 <- model_data %>%
  mutate(
    post_10 = t_post * (salinity_level == "10 ppt"),
    post_15 = t_post * (salinity_level == "15 ppt"),
    post_20 = t_post * (salinity_level == "20 ppt"),
    post_control =
      t_post * (salinity_level == "Control"),
    post_procedural =
      t_post *
      (salinity_level == "Procedural control")
  )

# Fit negative binomial GLMM
model_full <- glmmTMB(
  polyp_number ~
    salinity_level * (t_exp + t_rec + t_post) +
    log(initial_polyps) +
    (1 | replicate_id),
  family = nbinom2(link = "log"),
  data = model_data
)

summary(model_full)

model_reduced <- glmmTMB(
  polyp_number ~
    salinity_level +
    t_exp + t_rec + t_post +
    log(initial_polyps) +
    (1 | replicate_id),
  family = nbinom2(link = "log"),
  data = model_data
)

# Likelihood ratio test
anova(model_reduced, model_full)

# check assumptions

# Simulate residuals
res <- simulateResiduals(
  fittedModel = model_full,
  n = 1000
)

# Diagnostic plots
plot(res)

# Formal diagnostic tests
testDispersion(res)
testZeroInflation(res)
testUniformity(res)


# Example: inspect one replicate
one_id <- levels(model_data$replicate_id)[1]

one_rep <- model_data %>%
  filter(replicate_id == one_id) %>%
  arrange(Day)

# Plot raw time series
plot(
  one_rep$Day,
  one_rep$polyp_number,
  type = "b",
  xlab = "Day",
  ylab = "Polyp count"
)



# Exposure slopes
exp_slopes <- emtrends(
  model_full,
  ~ salinity_level,
  var = "t_exp"
)

summary(exp_slopes)

# Recovery slopes
rec_slopes <- emtrends(
  model_full,
  ~ salinity_level,
  var = "t_rec"
)

summary(rec_slopes)

# Post-recovery slopes
post_slopes <- emtrends(
  model_full,
  ~ salinity_level,
  var = "t_post"
)

summary(post_slopes)


# Example: exposure slopes
exp_results <- as.data.frame(exp_slopes)

exp_results <- exp_results %>%
  mutate(
    Percent_Change_Per_Day =
      100 * (exp(t_exp.trend) - 1)
  )

print(exp_results)

# Pairwise comparisons during exposure
pairs(exp_slopes, adjust = "tukey")

# Pairwise comparisons during recovery
pairs(rec_slopes, adjust = "tukey")

# Pairwise comparisons during post-recovery
pairs(post_slopes, adjust = "tukey")


model_full$fit$convergence
model_full$sdr$pdHess

levels(model_data$salinity_level)



model_data %>%
  filter(
    salinity_level == "5 ppt",
    phase == "Post-recovery"
  ) %>%
  group_by(replicate_id) %>%
  summarise(
    Minimum = min(polyp_number),
    Maximum = max(polyp_number),
    Mean = mean(polyp_number),
    Zero_Days = sum(polyp_number == 0),
    Total_Days = n(),
    .groups = "drop"
  )


# deal with the 5ppt
model_data <- model_data %>%
  mutate(
    is_5ppt = as.integer(salinity_level == "5 ppt")
  )

model_data <- model_data %>%
  mutate(
    post_10 = t_post * (salinity_level == "10 ppt"),
    post_15 = t_post * (salinity_level == "15 ppt"),
    post_20 = t_post * (salinity_level == "20 ppt"),
    post_control =
      t_post * (salinity_level == "Control"),
    post_procedural =
      t_post *
      (salinity_level == "Procedural control")
  )


model_modified <- glmmTMB(
  polyp_number ~
    salinity_level * (t_exp + t_rec) +
    post_10 + post_15 + post_20 +
    post_control + post_procedural +
    log(initial_polyps) +
    (1 | replicate_id),
  family = nbinom2(link = "log"),
  data = model_data
)

summary(model_modified)


# Convergence
model_modified$fit$convergence
model_modified$sdr$pdHess

# Compare information criteria
AIC(model_full, model_modified)

# Simulated residual diagnostics
library(DHARMa)

res_modified <- simulateResiduals(
  fittedModel = model_modified,
  n = 1000
)

plot(res_modified)

testDispersion(res_modified)
testZeroInflation(res_modified)
testUniformity(res_modified)

model_data %>%
  group_by(replicate_id) %>%
  summarise(
    n_days = n_distinct(Day),
    first_day = min(Day),
    last_day = max(Day),
    .groups = "drop"
  )

# Extract conditional Pearson residuals
model_data$residual_pearson <- residuals(
  model_modified,
  type = "pearson"
)

# Select one replicate
rep_example <- model_data %>%
  filter(
    replicate_id == levels(replicate_id)[7]
  ) %>%
  arrange(Day)

# Plot autocorrelation
acf(
  rep_example$residual_pearson,
  main = "Residual autocorrelation: example replicate"
)


five_ppt_predictions <- model_data %>%
  filter(
    salinity_level == "5 ppt",
    phase == "Post-Recovery"
  ) %>%
  mutate(
    Predicted = predict(
      model_modified,
      newdata = .,
      type = "response"
    )
  )

five_ppt_predictions %>%
  group_by(Day) %>%
  summarise(
    Observed_Mean = mean(polyp_number),
    Predicted_Mean = mean(Predicted),
    .groups = "drop"
  ) %>%
  print(n = Inf)



library(dplyr)
library(ggplot2)

# Extract residuals
model_data$residual_pearson <- residuals(
  model_modified,
  type = "pearson"
)

# Calculate lag-1 correlations
acf_summary <- model_data %>%
  arrange(replicate_id, Day) %>%
  group_by(replicate_id, salinity_level) %>%
  summarise(
    Lag1_ACF = {
      x <- residual_pearson
      
      if (
        length(x) >= 3 &&
        sd(x[-1]) > 0 &&
        sd(x[-length(x)]) > 0
      ) {
        cor(
          x[-1],
          x[-length(x)]
        )
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

# View results
print(acf_summary, n = Inf)

# Plot
ggplot(
  acf_summary,
  aes(
    x = salinity_level,
    y = Lag1_ACF
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_jitter(
    width = 0.1,
    height = 0,
    size = 2
  ) +
  labs(
    x = "Salinity treatment",
    y = "Lag-1 residual correlation",
    title = "Residual autocorrelation across replicates"
  ) +
  theme_classic()



# Simulate residuals from the modified model
res_modified <- DHARMa::simulateResiduals(
  model_modified,
  n = 1000
)

# Plot residuals against experimental day
DHARMa::plotResiduals(
  res_modified,
  form = model_data$Day
)

# Plot residuals against salinity treatment
DHARMa::plotResiduals(
  res_modified,
  form = model_data$salinity_level
)



library(dplyr)

polyp_abund <- polyp_abund %>%
  mutate(
    Phase = case_when(
      Day == 1 ~ "Baseline",
      Day <= 6 ~ "Exposure",
      Day <= 10 ~ "Recovery",
      Day >= 11 ~ "Post-recovery"
    ),
    
    Phase = factor(
      Phase,
      levels = c(
        "Baseline",
        "Exposure",
        "Recovery",
        "Post-recovery"
      )
    ),
    
    # Updated piecewise time variables
    t_exp = pmin(pmax(Day - 1, 0), 5),
    t_rec = pmin(pmax(Day - 6, 0), 4),
    t_post = pmax(Day - 10, 0)
  )

# Inspect the corrected structure
polyp_abund %>%
  distinct(Day, Phase, t_exp, t_rec, t_post) %>%
  arrange(Day)


model_data2 <- polyp_abund %>%
  filter(Day > 1) %>%
  droplevels()


model_data2 <- model_data2 %>%
  mutate(
    post_10 = t_post * (salinity_level == "10 ppt"),
    post_15 = t_post * (salinity_level == "15 ppt"),
    post_20 = t_post * (salinity_level == "20 ppt"),
    post_control =
      t_post * (salinity_level == "Control"),
    post_procedural =
      t_post *
      (salinity_level == "Procedural control")
  )


library(glmmTMB)

model_corrected <- glmmTMB(
  polyp_number ~
    salinity_level * (t_exp + t_rec) +
    post_10 + post_15 + post_20 +
    post_control + post_procedural +
    log(initial_polyps) +
    (1 | replicate_id),
  family = nbinom2(link = "log"),
  data = model_data2
)

summary(model_corrected)


# diagnostics

library(DHARMa)

# Check convergence
model_corrected$fit$convergence
model_corrected$sdr$pdHess

# Simulated residuals
res_corrected <- simulateResiduals(
  fittedModel = model_corrected,
  n = 1000
)

plot(res_corrected)

testDispersion(res_corrected)
testZeroInflation(res_corrected)
testUniformity(res_corrected)

# Residual patterns across time and treatments
plotResiduals(
  res_corrected,
  form = model_data$Day
)

plotResiduals(
  res_corrected,
  form = model_data$salinity_level
)




library(DHARMa)
library(dplyr)

# Test temporal autocorrelation separately
# within each replicate

replicate_ids <- unique(model_data2$replicate_id)

autocorrelation_results <- lapply(
  replicate_ids,
  function(id) {
    
    idx <- which(
      model_data2$replicate_id == id
    )
    
    # Extract residuals for this replicate
    res_subset <- recalculateResiduals(
      res_corrected,
      sel = idx
    )
    
    # Test temporal autocorrelation
    testTemporalAutocorrelation(
      res_subset,
      time = model_data2$Day[idx],
      plot = FALSE
    )
  }
)

# Extract p-values
acf_pvalues <- sapply(
  autocorrelation_results,
  function(x) x$p.value
)

summary(acf_pvalues)



# Number of significant tests
sum(acf_pvalues < 0.05, na.rm = TRUE)

# Proportion of significant tests
mean(acf_pvalues < 0.05, na.rm = TRUE)

# Identify affected replicates
data.frame(
  replicate_id = replicate_ids,
  p_value = acf_pvalues
) %>%
  filter(p_value < 0.05)



# Adjust p-values for multiple testing
acf_adjusted <- p.adjust(
  acf_pvalues,
  method = "BH"
)

# Create summary table
acf_results <- data.frame(
  replicate_id = replicate_ids,
  p_value = acf_pvalues,
  adjusted_p = acf_adjusted
)

# Display smallest p-values
acf_results %>%
  arrange(p_value)

# Number remaining significant
sum(acf_adjusted < 0.05, na.rm = TRUE)



library(dplyr)
library(ggplot2)

# Extract observed and fitted values
procedural_check <- model_data2 %>%
  filter(
    as.character(replicate_id) %in% c(
      "Procedural control.2",
      "Procedural control.3"
    )
  ) %>%
  mutate(
    Predicted = predict(
      model_corrected,
      newdata = .,
      type = "response"
    )
  )

# Plot observed and predicted trajectories
ggplot(procedural_check, aes(x = Day)) +
  geom_line(
    aes(y = polyp_number, colour = "Observed"),
    linewidth = 0.9
  ) +
  geom_point(
    aes(y = polyp_number, colour = "Observed"),
    size = 2
  ) +
  geom_line(
    aes(y = Predicted, colour = "Predicted"),
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  facet_wrap(~ replicate_id) +
  scale_colour_manual(
    values = c(
      "Observed" = "black",
      "Predicted" = "blue"
    )
  ) +
  labs(
    x = "Experimental day",
    y = "Number of polyps",
    colour = NULL,
    title = "Procedural control: observed vs predicted"
  ) +
  theme_classic()



library(emmeans)

# Estimate exposure slopes
exposure_slopes <- emtrends(
  model_corrected,
  ~ salinity_level,
  var = "t_exp"
)

# Display slopes and confidence intervals
summary(
  exposure_slopes,
  infer = c(TRUE, TRUE)
)


exposure_results <- as.data.frame(
  confint(exposure_slopes)
)

exposure_results <- exposure_results %>%
  mutate(
    Percent_Change = 100 * (
      exp(t_exp.trend) - 1
    ),
    Lower_Percent = 100 * (
      exp(asymp.LCL) - 1
    ),
    Upper_Percent = 100 * (
      exp(asymp.UCL) - 1
    )
  )

print(exposure_results)



# Compare each treatment against Control
exposure_contrasts <- contrast(
  exposure_slopes,
  method = "trt.vs.ctrl",
  ref = which(
    levels(model_data2$salinity_level) == "Control"
  ),
  adjust = "dunnettx"
)

summary(exposure_contrasts)


# did polyp polulations begin recovering while salinity was being restored


library(emmeans)
library(dplyr)

# Estimate recovery-phase slopes
recovery_slopes <- emtrends(
  model_corrected,
  ~ salinity_level,
  var = "t_rec"
)

# Display estimates and confidence intervals
summary(
  recovery_slopes,
  infer = c(TRUE, TRUE)
)

# Convert to percentage change per day
recovery_results <- as.data.frame(
  confint(recovery_slopes)
) %>%
  mutate(
    Percent_Change =
      100 * (exp(t_rec.trend) - 1),
    
    Lower_Percent =
      100 * (exp(asymp.LCL) - 1),
    
    Upper_Percent =
      100 * (exp(asymp.UCL) - 1)
  )

print(recovery_results)



recovery_contrasts <- contrast(
  recovery_slopes,
  method = "trt.vs.ctrl",
  ref = which(
    levels(model_data2$salinity_level) == "Control"
  ),
  adjust = "dunnettx"
)

summary(recovery_contrasts)


library(emmeans)

# Compare each recovery slope with the control
recovery_contrasts <- contrast(
  recovery_slopes,
  method = "trt.vs.ctrl",
  ref = which(
    levels(model_data2$salinity_level) == "Control"
  ),
  adjust = "dunnettx"
)

# View treatment comparisons
summary(
  recovery_contrasts,
  infer = TRUE
)


# post recovery


library(dplyr)

post_recovery_results <- data.frame(
  salinity_level = c(
    "10 ppt",
    "15 ppt",
    "20 ppt",
    "Control",
    "Procedural control"
  ),
  
  term = c(
    "post_10",
    "post_15",
    "post_20",
    "post_control",
    "post_procedural"
  )
)

# Extract coefficients and confidence intervals
coefs <- summary(model_corrected)$coefficients$cond

post_recovery_results <- post_recovery_results %>%
  mutate(
    Estimate = coefs[term, "Estimate"],
    SE = coefs[term, "Std. Error"],
    
    Lower = Estimate - 1.96 * SE,
    Upper = Estimate + 1.96 * SE,
    
    Percent_Change =
      100 * (exp(Estimate) - 1),
    
    Lower_Percent =
      100 * (exp(Lower) - 1),
    
    Upper_Percent =
      100 * (exp(Upper) - 1)
  )

print(post_recovery_results)


library(emmeans)

# Obtain the model coefficient estimates and covariance
beta <- fixef(model_corrected)$cond
V <- as.matrix(vcov(model_corrected)$cond)

# Post-recovery slope terms
post_terms <- c(
  "post_10",
  "post_15",
  "post_20",
  "post_control",
  "post_procedural"
)

# Construct contrasts against control
contrast_matrix <- sapply(
  setdiff(post_terms, "post_control"),
  function(term) {
    x <- setNames(rep(0, length(beta)), names(beta))
    x[term] <- 1
    x["post_control"] <- -1
    x
  }
)

# Contrast estimates
est <- as.vector(t(contrast_matrix) %*% beta)

# Standard errors
se <- sqrt(diag(
  t(contrast_matrix) %*%
    V %*%
    contrast_matrix
))

# Wald statistics and unadjusted p-values
z <- est / se
p <- 2 * pnorm(-abs(z))

post_contrasts <- data.frame(
  Comparison = paste(
    c("10 ppt", "15 ppt", "20 ppt",
      "Procedural control"),
    "- Control"
  ),
  Estimate = est,
  SE = se,
  z = z,
  p_unadjusted = p,
  p_Holm = p.adjust(p, method = "holm")
)

print(post_contrasts)

# look at final day

library(dplyr)

observed_day15 <- model_data2 %>%
  filter(Day == 15) %>%
  group_by(salinity_level) %>%
  summarise(
    n = n(),
    Mean_Count = mean(polyp_number),
    SD_Count = sd(polyp_number),
    Median_Count = median(polyp_number),
    Min_Count = min(polyp_number),
    Max_Count = max(polyp_number),
    Zero_Replicates = sum(polyp_number == 0),
    .groups = "drop"
  )

print(observed_day15)


final_recovery <- model_data2 %>%
  filter(Day == 15) %>%
  mutate(
    Final_Ratio = polyp_number / initial_polyps,
    Final_Percent = 100 * Final_Ratio
  )

recovery_summary <- final_recovery %>%
  group_by(salinity_level) %>%
  summarise(
    Mean_Percent = mean(Final_Percent),
    SD_Percent = sd(Final_Percent),
    Median_Percent = median(Final_Percent),
    Min_Percent = min(Final_Percent),
    Max_Percent = max(Final_Percent),
    .groups = "drop"
  )

print(recovery_summary)

# day 15 analysis


library(glmmTMB)
library(dplyr)

# Extract final observations
day15_model_data <- model_data2 %>%
  filter(Day == 15) %>%
  mutate(
    salinity_level = relevel(
      salinity_level,
      ref = "Control"
    )
  )

# Fit a negative binomial model
day15_model <- glmmTMB(
  polyp_number ~ salinity_level +
    log(initial_polyps),
  family = nbinom2(link = "log"),
  data = day15_model_data
)

summary(day15_model)

# remove 5ppt

library(glmmTMB)
library(dplyr)
library(emmeans)

day15_nonzero <- day15_model_data %>%
  filter(salinity_level != "5 ppt") %>%
  mutate(
    salinity_level = droplevels(salinity_level)
  )

day15_model2 <- glmmTMB(
  polyp_number ~ salinity_level +
    log(initial_polyps),
  family = nbinom2(link = "log"),
  data = day15_nonzero
)

summary(day15_model2)

# Check model convergence
day15_model2$fit$convergence
day15_model2$sdr$pdHess


# Estimate adjusted treatment means
day15_emmeans <- emmeans(
  day15_model2,
  ~ salinity_level
)

# Compare treatments against Control
day15_contrasts <- contrast(
  day15_emmeans,
  method = "trt.vs.ctrl",
  ref = which(
    levels(day15_nonzero$salinity_level) ==
      "Control"
  ),
  adjust = "dunnettx"
)

# Dunnett-adjusted tests and confidence intervals
summary(
  day15_contrasts,
  infer = c(TRUE, TRUE)
)

# Express treatment differences as abundance ratios
summary(
  day15_contrasts,
  type = "response",
  infer = c(TRUE, TRUE)
)


library(DHARMa)

# Check convergence
day15_model2$fit$convergence
day15_model2$sdr$pdHess

# Simulated residual diagnostics
set.seed(123)
day15_residuals <- simulateResiduals(
  fittedModel = day15_model2,
  n = 1000
)

plot(day15_residuals)

testDispersion(day15_residuals)
testZeroInflation(day15_residuals)
testUniformity(day15_residuals)


# sensitivity analysis

library(dplyr)
library(purrr)
library(glmmTMB)

# Identify replicate cultures
replicates <- unique(day15_nonzero$replicate_id)

# Fit model repeatedly, leaving out one replicate
loo_results <- map_dfr(replicates, function(rep_id) {
  
  reduced_data <- day15_nonzero %>%
    filter(replicate_id != rep_id)
  
  fit <- tryCatch(
    glmmTMB(
      polyp_number ~ salinity_level +
        log(initial_polyps),
      family = nbinom2(link = "log"),
      data = reduced_data
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble(
      Omitted_Replicate = as.character(rep_id),
      Estimate = NA_real_,
      Abundance_Ratio = NA_real_,
      Converged = FALSE
    ))
  }
  
  beta <- fixef(fit)$cond["salinity_level10 ppt"]
  
  tibble(
    Omitted_Replicate = as.character(rep_id),
    Estimate = unname(beta),
    Abundance_Ratio = exp(unname(beta)),
    Converged = (
      fit$fit$convergence == 0 &&
        isTRUE(fit$sdr$pdHess)
    )
  )
})

print(loo_results)


loo_summary <- loo_results %>%
  summarise(
    Number_of_Fits = n(),
    Number_Converged = sum(Converged),
    Min_Ratio = min(
      Abundance_Ratio[Converged],
      na.rm = TRUE
    ),
    Max_Ratio = max(
      Abundance_Ratio[Converged],
      na.rm = TRUE
    ),
    Median_Ratio = median(
      Abundance_Ratio[Converged],
      na.rm = TRUE
    )
  )

print(loo_summary)


loo_results %>%
  arrange(Abundance_Ratio) %>%
  print(n = 25)



library(dplyr)

day15_nonzero %>%
  filter(salinity_level == "10 ppt") %>%
  select(
    replicate_id,
    initial_polyps,
    polyp_number
  ) %>%
  mutate(
    Final_Percent =
      100 * polyp_number / initial_polyps
  ) %>%
  arrange(desc(Final_Percent))


# starting with corrected data ---------------

library(dplyr)

# Check number of observations and replicates
model_data2 %>%
  summarise(
    observations = n(),
    replicates = n_distinct(replicate_id),
    treatments = n_distinct(salinity_level)
  )

# Check observations per treatment and day
model_data2 %>%
  count(salinity_level, Day) %>%
  tidyr::pivot_wider(
    names_from = Day,
    values_from = n
  )

# Check missing values and count ranges
model_data2 %>%
  summarise(
    missing_counts = sum(is.na(polyp_number)),
    missing_initial = sum(is.na(initial_polyps)),
    minimum_count = min(polyp_number, na.rm = TRUE),
    maximum_count = max(polyp_number, na.rm = TRUE),
    minimum_initial = min(initial_polyps, na.rm = TRUE)
  )


# negative binomial glmm


library(dplyr)
library(glmmTMB)

model_data2 <- model_data2 %>%
  mutate(
    t_exp = pmin(pmax(Day - 1, 0), 5),
    t_rec = pmin(pmax(Day - 6, 0), 4),
    t_post = pmax(Day - 10, 0),
    
    post_10 = t_post *
      (salinity_level == "10 ppt"),
    
    post_15 = t_post *
      (salinity_level == "15 ppt"),
    
    post_20 = t_post *
      (salinity_level == "20 ppt"),
    
    post_control = t_post *
      (salinity_level == "Control"),
    
    post_procedural = t_post *
      (salinity_level == "Procedural control")
  )



model_data2 %>%
  filter(Day >= 11) %>%
  group_by(salinity_level) %>%
  summarise(
    Total_Polyps = sum(polyp_number),
    Positive_Observations = sum(polyp_number > 0),
    Zero_Observations = sum(polyp_number == 0),
    .groups = "drop"
  )


library(glmmTMB)

model_corrected <- glmmTMB(
  polyp_number ~
    salinity_level * (t_exp + t_rec) +
    post_10 +
    post_15 +
    post_20 +
    post_control +
    post_procedural +
    log(initial_polyps) +
    (1 | replicate_id),
  
  family = nbinom2(link = "log"),
  data = model_data2
)

summary(model_corrected)

# Check convergence and Hessian
model_corrected$fit$convergence
model_corrected$sdr$pdHess

# Model fit statistics
AIC(model_corrected)


# try a poisson


library(glmmTMB)

model_poisson <- glmmTMB(
  polyp_number ~
    salinity_level * (t_exp + t_rec) +
    post_10 +
    post_15 +
    post_20 +
    post_control +
    post_procedural +
    log(initial_polyps) +
    (1 | replicate_id),
  
  family = poisson(link = "log"),
  data = model_data2
)

summary(model_poisson)

# Check numerical stability
model_poisson$fit$convergence
model_poisson$sdr$pdHess

# Model fit statistics
AIC(model_poisson)



library(DHARMa)

set.seed(123)

res_poisson <- simulateResiduals(
  fittedModel = model_poisson,
  n = 1000
)

plot(res_poisson)

testDispersion(res_poisson)
testZeroInflation(res_poisson)
testUniformity(res_poisson)



library(DHARMa)

# Residual patterns across treatments
plotResiduals(
  res_poisson,
  form = model_data2$salinity_level
)

# Residual patterns across sampling days
plotResiduals(
  res_poisson,
  form = model_data2$Day,
  asFactor = TRUE
)


testCategorical(
  res_poisson,
  catPred = model_data2$salinity_level
)


testCategorical(
  res_poisson,
  catPred = factor(model_data2$Day)
)

# check temporal autocorrelation


library(dplyr)
library(DHARMa)
library(purrr)

# Identify observation positions for each culture
replicate_indices <- split(
  seq_len(nrow(model_data2)),
  model_data2$replicate_id
)

# Test temporal autocorrelation separately for each replicate
autocorrelation_results <- map_dfr(
  names(replicate_indices),
  function(rep_id) {
    
    idx <- replicate_indices[[rep_id]]
    
    idx <- idx[order(model_data2$Day[idx])]
    
    residual_subset <- recalculateResiduals(
      res_poisson,
      sel = idx
    )
    
    test_result <- testTemporalAutocorrelation(
      residual_subset,
      time = model_data2$Day[idx],
      plot = FALSE
    )
    
    tibble(
      replicate_id = rep_id,
      DW_statistic = unname(test_result$statistic),
      p_value = test_result$p.value
    )
  }
)

# Adjust for testing 30 cultures
autocorrelation_results <- autocorrelation_results %>%
  mutate(
    adjusted_p = p.adjust(p_value, method = "BH")
  )

print(autocorrelation_results, n = 30)


autocorrelation_results %>%
  summarise(
    Total_Replicates = n(),
    Significant_Unadjusted =
      sum(p_value < 0.05, na.rm = TRUE),
    Significant_Adjusted =
      sum(adjusted_p < 0.05, na.rm = TRUE),
    Minimum_p = min(p_value, na.rm = TRUE),
    Minimum_Adjusted_p =
      min(adjusted_p, na.rm = TRUE)
  )

# which rep is causing autocorr issues


autocorrelation_results %>%
  arrange(adjusted_p) %>%
  print(n = 10)


# Identify the replicate with the smallest adjusted p-value
rep_id <- autocorrelation_results %>%
  arrange(adjusted_p) %>%
  slice(1) %>%
  pull(replicate_id)

# Display its observed counts over time
model_data2 %>%
  filter(as.character(replicate_id) == rep_id) %>%
  arrange(Day) %>%
  select(
    replicate_id,
    salinity_level,
    Day,
    polyp_number
  )


library(dplyr)
library(ggplot2)

trajectory_check <- model_data2 %>%
  mutate(
    fitted_count = predict(
      model_poisson,
      type = "response",
      re.form = NULL
    )
  ) %>%
  filter(
    as.character(replicate_id) ==
      "Procedural control.2"
  ) %>%
  arrange(Day)

ggplot(
  trajectory_check,
  aes(x = Day)
) +
  geom_point(
    aes(y = polyp_number),
    size = 3
  ) +
  geom_line(
    aes(y = polyp_number),
    linewidth = 1
  ) +
  geom_line(
    aes(y = fitted_count),
    colour = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    title = "Procedural control.2",
    subtitle = "Observed versus fitted polyp abundance",
    x = "Experimental day",
    y = "Number of polyps"
  ) +
  theme_classic()


model_without_pc2 <- update(
  model_poisson,
  data = model_data2 %>%
    filter(
      as.character(replicate_id) !=
        "Procedural control.2"
    ) %>%
    droplevels()
)

# Check convergence
model_without_pc2$fit$convergence
model_without_pc2$sdr$pdHess

# Repeat DHARMa diagnostics
set.seed(123)

res_without_pc2 <- simulateResiduals(
  model_without_pc2,
  n = 1000
)

testUniformity(res_without_pc2)
testDispersion(res_without_pc2)
testZeroInflation(res_without_pc2)


# random slope model


library(glmmTMB)
library(dplyr)

# Centre and scale experimental day
model_data2 <- model_data2 %>%
  mutate(
    Day_scaled = (Day - mean(Day)) / sd(Day)
  )

# Fit Poisson GLMM with replicate-specific time slopes
model_random_slope <- glmmTMB(
  polyp_number ~
    salinity_level * (t_exp + t_rec) +
    post_10 +
    post_15 +
    post_20 +
    post_control +
    post_procedural +
    log(initial_polyps) +
    (1 + Day_scaled | replicate_id),
  
  family = poisson(link = "log"),
  data = model_data2
)

summary(model_random_slope)

# Check convergence
model_random_slope$fit$convergence
model_random_slope$sdr$pdHess

# Compare AIC values
AIC(model_poisson, model_random_slope)



library(DHARMa)

set.seed(123)

res_random_slope <- simulateResiduals(
  fittedModel = model_random_slope,
  n = 1000
)

plot(res_random_slope)

testUniformity(res_random_slope)
testDispersion(res_random_slope)
testZeroInflation(res_random_slope)

testCategorical(
  res_random_slope,
  catPred = model_data2$salinity_level
)

# uncorrelated random slope model


library(glmmTMB)

model_random_slope_uncorrelated <- glmmTMB(
  polyp_number ~
    salinity_level * (t_exp + t_rec) +
    post_10 +
    post_15 +
    post_20 +
    post_control +
    post_procedural +
    log(initial_polyps) +
    (1 | replicate_id) +
    (0 + Day_scaled | replicate_id),
  
  family = poisson(link = "log"),
  data = model_data2
)

summary(model_random_slope_uncorrelated)

# Check numerical stability
model_random_slope_uncorrelated$fit$convergence
model_random_slope_uncorrelated$sdr$pdHess

# Compare all three Poisson models
AIC(
  model_poisson,
  model_random_slope,
  model_random_slope_uncorrelated
)


library(DHARMa)

set.seed(123)

res_random_slope_uncorrelated <- simulateResiduals(
  fittedModel = model_random_slope_uncorrelated,
  n = 1000
)

plot(res_random_slope_uncorrelated)

testUniformity(res_random_slope_uncorrelated)
testDispersion(res_random_slope_uncorrelated)
testZeroInflation(res_random_slope_uncorrelated)

testCategorical(
  res_random_slope_uncorrelated,
  catPred = model_data2$salinity_level
)


library(glmmTMB)

model_nb_random_slope <- glmmTMB(
  polyp_number ~
    salinity_level * (t_exp + t_rec) +
    post_10 +
    post_15 +
    post_20 +
    post_control +
    post_procedural +
    log(initial_polyps) +
    (1 + Day_scaled | replicate_id),
  
  family = nbinom2(link = "log"),
  data = model_data2
)

summary(model_nb_random_slope)

# Check numerical stability
model_nb_random_slope$fit$convergence
model_nb_random_slope$sdr$pdHess

# Check dispersion
sigma(model_nb_random_slope)

# Compare AIC, if available
AIC(
  model_poisson,
  model_random_slope,
  model_nb_random_slope
)


# sensitivity check

library(emmeans)

exposure_intercept <- emtrends(
  model_poisson,
  ~ salinity_level,
  var = "t_exp"
)

exposure_slope <- emtrends(
  model_random_slope,
  ~ salinity_level,
  var = "t_exp"
)

# Dunnett-adjusted comparisons against Control
contrast(
  exposure_intercept,
  method = "trt.vs.ctrl",
  ref = "Control",
  adjust = "dunnett"
)

contrast(
  exposure_slope,
  method = "trt.vs.ctrl",
  ref = "Control",
  adjust = "dunnett"
)



recovery_intercept <- emtrends(
  model_poisson,
  ~ salinity_level,
  var = "t_rec"
)

recovery_slope <- emtrends(
  model_random_slope,
  ~ salinity_level,
  var = "t_rec"
)

contrast(
  recovery_intercept,
  method = "trt.vs.ctrl",
  ref = "Control",
  adjust = "dunnett"
)

contrast(
  recovery_slope,
  method = "trt.vs.ctrl",
  ref = "Control",
  adjust = "dunnett"
)


post_contrasts <- function(model) {
  
  b <- fixef(model)$cond
  V <- vcov(model)$cond
  
  treatments <- c(
    "10 ppt",
    "15 ppt",
    "20 ppt",
    "Procedural control"
  )
  
  terms <- c(
    "post_10",
    "post_15",
    "post_20",
    "post_procedural"
  )
  
  results <- lapply(seq_along(terms), function(i) {
    
    L <- setNames(rep(0, length(b)), names(b))
    
    L[terms[i]] <- 1
    L["post_control"] <- -1
    
    estimate <- sum(L * b)
    se <- sqrt(
      as.numeric(t(L) %*% V %*% L)
    )
    
    z <- estimate / se
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    
    data.frame(
      Treatment = treatments[i],
      Estimate = estimate,
      SE = se,
      z = z,
      p = p
    )
  })
  
  results <- dplyr::bind_rows(results)
  
  results$Holm_p <- p.adjust(
    results$p,
    method = "holm"
  )
  
  results
}

# Random-intercept model
post_contrasts(model_poisson)

# Random-slope model
post_contrasts(model_random_slope)


# leave-one-replicate-out (LOO) sensitivity analysis

library(dplyr)
library(purrr)
library(glmmTMB)
library(emmeans)

replicates <- unique(
  as.character(model_data2$replicate_id)
)

loo_recovery <- map_dfr(replicates, function(rep_id) {
  
  # Remove all observations from one culture
  reduced_data <- model_data2 %>%
    filter(
      as.character(replicate_id) != rep_id
    ) %>%
    droplevels()
  
  # Refit the original random-intercept Poisson GLMM
  fit <- tryCatch(
    update(
      model_poisson,
      data = reduced_data
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble(
      Omitted_Replicate = rep_id,
      Estimate = NA_real_,
      SE = NA_real_,
      Adjusted_p = NA_real_,
      Converged = FALSE
    ))
  }
  
  converged <- (
    fit$fit$convergence == 0 &&
      isTRUE(fit$sdr$pdHess)
  )
  
  if (!converged) {
    return(tibble(
      Omitted_Replicate = rep_id,
      Estimate = NA_real_,
      SE = NA_real_,
      Adjusted_p = NA_real_,
      Converged = FALSE
    ))
  }
  
  # Estimate recovery-phase slopes
  trends <- emtrends(
    fit,
    ~ salinity_level,
    var = "t_rec"
  )
  
  # Dunnett-adjusted comparisons with Control
  comparisons <- as.data.frame(
    summary(
      contrast(
        trends,
        method = "trt.vs.ctrl",
        ref = "Control",
        adjust = "dunnett"
      )
    )
  )
  
  # Extract 10 ppt vs Control
  result <- comparisons %>%
    filter(contrast == "10 ppt - Control")
  
  tibble(
    Omitted_Replicate = rep_id,
    Estimate = result$estimate,
    SE = result$SE,
    Adjusted_p = result$p.value,
    Converged = TRUE
  )
})

print(loo_recovery, n = 30)


loo_recovery_summary <- loo_recovery %>%
  summarise(
    Total_Fits = n(),
    
    Converged_Fits = sum(Converged),
    
    Minimum_Estimate = min(
      Estimate,
      na.rm = TRUE
    ),
    
    Maximum_Estimate = max(
      Estimate,
      na.rm = TRUE
    ),
    
    Median_Estimate = median(
      Estimate,
      na.rm = TRUE
    ),
    
    Significant_Fits = sum(
      Adjusted_p < 0.05,
      na.rm = TRUE
    ),
    
    Minimum_Adjusted_p = min(
      Adjusted_p,
      na.rm = TRUE
    ),
    
    Maximum_Adjusted_p = max(
      Adjusted_p,
      na.rm = TRUE
    )
  )

print(loo_recovery_summary)



loo_recovery %>%
  mutate(
    Change_From_Original = Estimate - (-0.406932)
  ) %>%
  arrange(
    desc(abs(Change_From_Original))
  ) %>%
  print(n = 30)



library(dplyr)
library(tidyr)

model_data2 %>%
  filter(
    salinity_level %in% c("10 ppt", "Control"),
    Day %in% c(6, 7, 8, 11)
  ) %>%
  select(
    replicate_id,
    salinity_level,
    Day,
    polyp_number
  ) %>%
  pivot_wider(
    names_from = Day,
    values_from = polyp_number,
    names_prefix = "Day_"
  ) %>%
  arrange(
    salinity_level,
    replicate_id
  )

# testing the recovery / post-recovery boundary


library(dplyr)
library(purrr)
library(glmmTMB)
library(emmeans)

# Create a dataset with a specified boundary
make_boundary_data <- function(B) {
  
  model_data2 %>%
    mutate(
      t_rec = pmin(
        pmax(Day - 6, 0),
        B - 6
      ),
      
      t_post = pmax(Day - B, 0),
      
      post_10 = t_post *
        (salinity_level == "10 ppt"),
      
      post_15 = t_post *
        (salinity_level == "15 ppt"),
      
      post_20 = t_post *
        (salinity_level == "20 ppt"),
      
      post_control = t_post *
        (salinity_level == "Control"),
      
      post_procedural = t_post *
        (salinity_level == "Procedural control")
    )
}

# Fit the same fixed/random effects structure
fit_boundary <- function(B) {
  
  glmmTMB(
    polyp_number ~
      salinity_level * (t_exp + t_rec) +
      post_10 + post_15 + post_20 +
      post_control + post_procedural +
      log(initial_polyps) +
      (1 | replicate_id),
    
    family = poisson(link = "log"),
    
    data = make_boundary_data(B)
  )
}

boundary_models <- list(
  Day8 = fit_boundary(8),
  Day10 = fit_boundary(10),
  Day11 = fit_boundary(11)
)



boundary_fit_summary <- imap_dfr(
  boundary_models,
  function(model, name) {
    
    tibble(
      Boundary = name,
      AIC = AIC(model),
      LogLik = as.numeric(logLik(model)),
      Convergence = model$fit$convergence,
      Positive_Hessian = isTRUE(model$sdr$pdHess)
    )
  }
)

print(boundary_fit_summary)



get_recovery_contrast <- function(model, name) {
  
  trends <- emtrends(
    model,
    ~ salinity_level,
    var = "t_rec"
  )
  
  result <- as.data.frame(
    summary(
      contrast(
        trends,
        method = "trt.vs.ctrl",
        ref = "Control",
        adjust = "dunnett"
      )
    )
  )
  
  result %>%
    filter(contrast == "10 ppt - Control") %>%
    transmute(
      Boundary = name,
      Estimate = estimate,
      SE = SE,
      z = z.ratio,
      Adjusted_p = p.value
    )
}

boundary_recovery <- imap_dfr(
  boundary_models,
  get_recovery_contrast
)

print(boundary_recovery)



boundary_post <- imap_dfr(
  boundary_models,
  function(model, name) {
    
    post_contrasts(model) %>%
      mutate(Boundary = name, .before = 1)
  }
)

boundary_post %>%
  filter(Treatment == "10 ppt") %>%
  select(
    Boundary,
    Treatment,
    Estimate,
    SE,
    p,
    Holm_p
  ) %>%
  print()

# day 15 analysis------



library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(DHARMa)

day15_data <- model_data2 %>%
  filter(Day == 15) %>%
  select(
    replicate_id,
    salinity_level,
    initial_polyps,
    polyp_number
  ) %>%
  droplevels()

# Verify the dataset
nrow(day15_data)

day15_data %>%
  count(salinity_level)

# Examine individual culture counts
day15_data %>%
  arrange(salinity_level, replicate_id)


day15_summary <- day15_data %>%
  group_by(salinity_level) %>%
  summarise(
    n = n(),
    Mean = mean(polyp_number),
    SD = sd(polyp_number),
    Median = median(polyp_number),
    Minimum = min(polyp_number),
    Maximum = max(polyp_number),
    Zero_Counts = sum(polyp_number == 0),
    Mean_Initial = mean(initial_polyps),
    .groups = "drop"
  )

print(day15_summary)



ggplot(
  day15_data,
  aes(
    x = salinity_level,
    y = polyp_number
  )
) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    fill = "grey85"
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2.5
  ) +
  labs(
    x = "Salinity treatment",
    y = "Polyp abundance at Day 15",
    title = "Final polyp abundance after salinity restoration"
  ) +
  theme_classic()



# Exclude 5 ppt from the count regression
day15_model_data <- day15_data %>%
  filter(salinity_level != "5 ppt") %>%
  droplevels()

# Poisson GLM
day15_poisson <- glm(
  polyp_number ~
    salinity_level +
    log(initial_polyps),
  family = poisson(link = "log"),
  data = day15_model_data
)

# Negative binomial GLM
day15_nb <- glmmTMB(
  polyp_number ~
    salinity_level +
    log(initial_polyps),
  family = nbinom2(link = "log"),
  data = day15_model_data
)

# Numerical checks
summary(day15_poisson)
summary(day15_nb)

day15_nb$fit$convergence
day15_nb$sdr$pdHess
sigma(day15_nb)

# Compare AIC
AIC(day15_poisson, day15_nb)


set.seed(123)

res_day15_poisson <- simulateResiduals(
  fittedModel = day15_poisson,
  n = 1000
)

res_day15_nb <- simulateResiduals(
  fittedModel = day15_nb,
  n = 1000
)

# Poisson diagnostics
plot(res_day15_poisson)
testUniformity(res_day15_poisson)
testDispersion(res_day15_poisson)
testZeroInflation(res_day15_poisson)

# Negative binomial diagnostics
plot(res_day15_nb)
testUniformity(res_day15_nb)
testDispersion(res_day15_nb)
testZeroInflation(res_day15_nb)



library(DHARMa)

set.seed(123)

res_day15_nb <- simulateResiduals(
  fittedModel = day15_nb,
  n = 1000
)

plot(res_day15_nb)

testUniformity(res_day15_nb)
testDispersion(res_day15_nb)
testZeroInflation(res_day15_nb)

# Check for treatment-specific residual problems
testCategorical(
  res_day15_nb,
  catPred = day15_model_data$salinity_level
)



library(emmeans)

# Estimated marginal means
day15_emm <- emmeans(
  day15_nb,
  ~ salinity_level
)

# Dunnett-adjusted comparisons against Control
day15_dunnett <- contrast(
  day15_emm,
  method = "trt.vs.ctrl",
  ref = "Control",
  adjust = "dunnett"
)

# Log-scale treatment comparisons
summary(day15_dunnett)

# Abundance ratios with confidence intervals
summary(
  day15_dunnett,
  type = "response",
  infer = c(TRUE, TRUE)
)

# Model-adjusted abundance estimates
summary(
  day15_emm,
  type = "response"
)


# Negative binomial model with initial abundance as an offset
day15_nb_offset <- glmmTMB(
  polyp_number ~
    salinity_level +
    offset(log(initial_polyps)),
  family = nbinom2(link = "log"),
  data = day15_model_data
)

# Numerical checks
summary(day15_nb_offset)
day15_nb_offset$fit$convergence
day15_nb_offset$sdr$pdHess

# Compare model fit
AIC(day15_nb, day15_nb_offset)

# Dunnett-adjusted comparisons
offset_emm <- emmeans(
  day15_nb_offset,
  ~ salinity_level
)

offset_dunnett <- contrast(
  offset_emm,
  method = "trt.vs.ctrl",
  ref = "Control",
  adjust = "dunnett"
)

summary(
  offset_dunnett,
  type = "response",
  infer = c(TRUE, TRUE)
)


library(DHARMa)

# Numerical stability
day15_nb_offset$fit$convergence
day15_nb_offset$sdr$pdHess
sigma(day15_nb_offset)

# Simulation-based diagnostics
set.seed(123)

res_day15_offset <- simulateResiduals(
  day15_nb_offset,
  n = 1000
)

plot(res_day15_offset)

testUniformity(res_day15_offset)
testDispersion(res_day15_offset)
testZeroInflation(res_day15_offset)

testCategorical(
  res_day15_offset,
  catPred = day15_model_data$salinity_level
)
