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
    has_day1 = any(Day == 1),
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
    relative_abundance = if_else(
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
      select(Day, salinity_level, replicate, budding_ratio)
  )
  
  p <- ggplot(dat, aes(x = Day)) +
    
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
      breaks = 1:15
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
      "output/",
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
    arrange(replicate, Day)
  
  scale_factor <- max(dat$budding_ratio, na.rm = TRUE) / 40
  
  phase_data <- tibble(
    xmin = c(1.5, 5.5),
    xmax = c(5.5, 10.5),
    Phase = c("Exposure", "Recovery")
  )
  
  p <- ggplot(dat, aes(x = Day)) +
    
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
      "output/",
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

treatment_data <- treatment_data %>%
  rename(salinity_level = Treatment)

scale_factor <- max(polyp_abund$budding_ratio, na.rm = TRUE) / 40

# mean polyp values per treatment
mean_data_budding <- polyp_abund %>%
  group_by(salinity_level, Day) %>%
  summarise(
    mean_budding = mean(budding_ratio, na.rm = TRUE),
    se = sd(budding_ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot(
  plot_data,
  aes(
    Day,
    budding_ratio,
    colour = factor(replicate),
    group = replicate
  )
) +
  
  # replicate trajectories
  geom_line(
    linewidth = 1,
    alpha = 0.5
  ) +
  
  geom_point(
    size = 2,
    alpha = 0.5
  ) +
  
  # salinity profile
  geom_step(
    data = treatment_data,
    aes(
      x = Day,
      y = Salinity * scale_factor,
      group = salinity_level
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  
  # treatment mean
  stat_summary(
    aes(
      group = 1
    ),
    fun = mean,
    geom = "line",
    colour = "black",
    linewidth = 1.5
  ) +
  
  stat_summary(
    aes(
      group = 1
    ),
    fun = mean,
    geom = "point",
    colour = "black",
    size = 3
  ) +
  
  facet_wrap(
    ~ salinity_level,
    scales = "free_y"
  ) +
  
  scale_y_continuous(
    name = "Ratio of polyps with buds",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Salinity (ppt)"
    )
  ) +
  
  scale_x_continuous(
    breaks = 1:15
  ) +
  
  
  labs(
    x = "Experiment Day",
    y = "Ratio of polyps with buds",
    colour = "Replicate"
  ) +
  
  theme_bw()




p <- ggplot(
  polyp_abund,
  aes(
    x = Day,
    y = budding_ratio,
    colour = factor(replicate),
    group = replicate
  )
) +
  
  # Exposure period
  annotate(
    "rect",
    xmin = 1.5,
    xmax = 5.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "firebrick",
    alpha = 0.08
  ) +
  
  # Recovery period
  annotate(
    "rect",
    xmin = 5.5,
    xmax = 10.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "steelblue",
    alpha = 0.08
  ) +
  
  # SE ribbon
  geom_ribbon(
    data = mean_data_budding,
    aes(
      x = Day,
      ymin = mean_budding - se,
      ymax = mean_budding + se,
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
      x = Day,
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
      x = Day,
      y = mean_budding
    ),
    inherit.aes = FALSE,
    colour = "black",
    size = 2.5
  ) +
  
  # Phase labels
  annotate(
    "text",
    x = 3.5,
    y = Inf,
    label = "Exposure",
    vjust = 1.5,
    fontface = "bold",
    colour = "firebrick"
  ) +
  
  annotate(
    "text",
    x = 8,
    y = Inf,
    label = "Recovery",
    vjust = 1.5,
    fontface = "bold",
    colour = "steelblue"
  ) +
  
  annotate(
    "text",
    x = 13,
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
    breaks = 1:15
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
  
  
  scale_factor <- max(dat$relative_abundance, na.rm = TRUE) / 40
  
  print(
    dat %>%
      select(Day, salinity_level, replicate, relative_abundance)
  )
  
  p <- ggplot(dat, aes(x = Day)) +
    
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
      name = "Relative abundance of polyps",
      
      sec.axis = sec_axis(
        ~ . / scale_factor,
        name = "Salinity (ppt)"
      )
    ) +
    
    scale_x_continuous(
      breaks = 1:15
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
    arrange(replicate, Day)
  
  scale_factor <- max(dat$relative_abundance, na.rm = TRUE) / 40
  
  phase_data <- tibble(
    xmin = c(1.5, 5.5),
    xmax = c(5.5, 10.5),
    Phase = c("Exposure", "Recovery")
  )
  
  p <- ggplot(dat, aes(x = Day)) +
    
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
        y = relative_abundance,
        colour = factor(replicate),
        group = replicate
      ),
      linewidth = 1
    ) +
    
    geom_point(
      aes(
        y = relative_abundance,
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

scale_factor <- max(polyp_abund$relative_abundance, na.rm = TRUE) / 40

# mean polyp values per treatment
mean_data_abund <- polyp_abund %>%
  group_by(salinity_level, Day) %>%
  summarise(
    mean_abund = mean(relative_abundance, na.rm = TRUE),
    se = sd(relative_abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot(
  polyp_abund,
  aes(
    Day,
    relative_abundance,
    colour = factor(replicate),
    group = replicate
  )
) +
  
  # replicate trajectories
  geom_line(
    linewidth = 1,
    alpha = 0.5
  ) +
  
  geom_point(
    size = 2,
    alpha = 0.5
  ) +
  
  # salinity profile
  geom_step(
    data = treatment_data,
    aes(
      x = Day,
      y = Salinity * scale_factor,
      group = salinity_level
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  
  # treatment mean
  stat_summary(
    aes(
      group = 1
    ),
    fun = mean,
    geom = "line",
    colour = "black",
    linewidth = 1.5
  ) +
  
  stat_summary(
    aes(
      group = 1
    ),
    fun = mean,
    geom = "point",
    colour = "black",
    size = 3
  ) +
  
  facet_wrap(
    ~ salinity_level,
    scales = "free_y"
  ) +
  
  scale_y_continuous(
    name = "Relative polyp abundance",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Salinity (ppt)"
    )
  ) +
  
  scale_x_continuous(
    breaks = 1:15
  ) +
  
  
  labs(
    x = "Experiment Day",
    y = "Relative polyp abundance",
    colour = "Replicate"
  ) +
  
  theme_bw()




p <- ggplot(
  polyp_abund,
  aes(
    x = Day,
    y = relative_abundance,
    colour = factor(replicate),
    group = replicate
  )
) +
  
  # Exposure period
  annotate(
    "rect",
    xmin = 1.5,
    xmax = 5.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "firebrick",
    alpha = 0.08
  ) +
  
  # Recovery period
  annotate(
    "rect",
    xmin = 5.5,
    xmax = 10.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "steelblue",
    alpha = 0.08
  ) +
  
  # SE ribbon
  geom_ribbon(
    data = mean_data_abund,
    aes(
      x = Day,
      ymin = mean_abund - se,
      ymax = mean_abund + se,
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
      x = Day,
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
      x = Day,
      y = mean_abund
    ),
    inherit.aes = FALSE,
    colour = "black",
    size = 2.5
  ) +
  
  # Phase labels
  annotate(
    "text",
    x = 3.5,
    y = Inf,
    label = "Exposure",
    vjust = 1.5,
    fontface = "bold",
    colour = "firebrick"
  ) +
  
  annotate(
    "text",
    x = 8,
    y = Inf,
    label = "Recovery",
    vjust = 1.5,
    fontface = "bold",
    colour = "steelblue"
  ) +
  
  annotate(
    "text",
    x = 13,
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
    breaks = 1:15
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
# create unique ids for each replicate
polyp_abund <- polyp_abund %>%
  mutate(
    replicate_id = interaction(salinity_level, replicate)
  )


m1 <- glm(relative_abundance ~ salinity_level, data = polyp_abund)
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
