# results_plotting.R
#
# Plotting raw results of experiment
#
# Created: 10/11/2026
# Created by: Davina Gifford

# load libraries
library(tidyverse)
library(cowplot)
library(png)

# load data
polyp_data <- readRDS("var/polyp_tidy.rds")
treatment_data <- readRDS("var/treatment_tidy.rds")
medusa_data <- readRDS("var/medusae_tidy.rds")

treatment_data <- treatment_data %>%
  rename(salinity_level = Treatment)


# Polyps ------------------------------------------------------------------



# combine into a single var for plotting

plot_data <- polyp_data %>%
  left_join(
    treatment_data,
    by = c(
      "salinity_level" = "Treatment",
      "Day" = "Day"
    )
  )

# plot the polyp numbers for the individual replicates 

replicates <- plot_data %>%
  distinct(salinity_level, replicate)

for(i in seq_len(nrow(replicates))) {
  
  treatment_i <- replicates$salinity_level[i]
  replicate_i <- replicates$replicate[i]
  
  dat <- plot_data %>%
    filter(
      salinity_level == treatment_i,
      replicate == replicate_i
    )
  
  
  scale_factor <- max(dat$polyp_number, na.rm = TRUE) / 40
  
  print(
    dat %>%
      select(Day, salinity_level, replicate, polyp_number)
  )
  
  p <- ggplot(dat, aes(x = Day)) +
    
    geom_line(
      aes(y = polyp_number),
      colour = "#0072B2",
      linewidth = 1
    ) +
    
    geom_point(
      aes(y = polyp_number),
      colour = "#0072B2",
      size = 2
    ) +
    
    geom_step(
      aes(y = Salinity * scale_factor),
      colour = "#D55E00",
      linewidth = 1
    ) +
    
    scale_y_continuous(
      name = "Number of polyps",
      
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
      "_rep",
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


for(treatment_i in unique(plot_data$salinity_level)) {
  
  dat <- plot_data %>%
    filter(salinity_level == treatment_i) %>%
    arrange(replicate, Day)
  
  scale_factor <- max(dat$polyp_number, na.rm = TRUE) / 40
  
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
        y = polyp_number,
        colour = factor(replicate),
        group = replicate
      ),
      linewidth = 1
    ) +
    
    geom_point(
      aes(
        y = polyp_number,
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
      name = "Number of Polyps",
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
      "_polyps.png"
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

scale_factor <- max(plot_data$polyp_number, na.rm = TRUE) / 40

# mean polyp values per treatment
mean_data <- plot_data %>%
  group_by(salinity_level, Day) %>%
  summarise(
    mean_polyps = mean(polyp_number, na.rm = TRUE),
    se = sd(polyp_number, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot(
  plot_data,
  aes(
    Day,
    polyp_number,
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
    name = "Number of Polyps",
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
    y = "Number of Polyps",
    colour = "Replicate"
  ) +
  
  theme_bw()

  


p <- ggplot(
  plot_data,
  aes(
    x = Day,
    y = polyp_number,
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
    data = mean_data,
    aes(
      x = Day,
      ymin = mean_polyps - se,
      ymax = mean_polyps + se,
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
    data = mean_data,
    aes(
      x = Day,
      y = mean_polyps,
      group = 1
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.5
  ) +
  
  geom_point(
    data = mean_data,
    aes(
      x = Day,
      y = mean_polyps
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
    y = "Number of Polyps"
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
    "output/polyp_results_by_treatment.png",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)


# Medusa ------------------------------------------------------------------



# combine into a single var for plotting

plot_data_m <- medusa_data %>%
  left_join(
    treatment_data,
    by = c(
      "salinity_level" = "salinity_level",
      "Day" = "Day"
    )
  )

# plot the medusae diameter for the individual replicates 

replicates_m <- plot_data_m %>%
  distinct(salinity_level, replicate)

for(i in seq_len(nrow(replicates_m))) {
  
  treatment_i <- replicates_m$salinity_level[i]
  replicate_i <- replicates_m$replicate[i]
  
  dat <- plot_data_m %>%
    filter(
      salinity_level == treatment_i,
      replicate == replicate_i
    )
  
  
  scale_factor <- max(dat$bell_diam, na.rm = TRUE) / 40
  
  print(
    dat %>%
      select(Day, salinity_level, replicate, bell_diam)
  )
  
  p <- ggplot(dat, aes(x = Day)) +
    
    geom_line(
      aes(y = bell_diam),
      colour = "#0072B2",
      linewidth = 1
    ) +
    
    geom_point(
      aes(y = bell_diam),
      colour = "#0072B2",
      size = 2
    ) +
    
    geom_step(
      aes(y = Salinity * scale_factor),
      colour = "#D55E00",
      linewidth = 1
    ) +
    
    scale_y_continuous(
      name = "Bell diameter (mm)",
      
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
      "bell_diam_rep",
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


for(treatment_i in unique(plot_data_m$salinity_level)) {
  
  dat <- plot_data_m %>%
    filter(salinity_level == treatment_i) %>%
    arrange(replicate, Day)
  
  scale_factor <- max(dat$bell_diam, na.rm = TRUE) / 40
  
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
        y = bell_diam,
        colour = factor(replicate),
        group = replicate
      ),
      linewidth = 1
    ) +
    
    geom_point(
      aes(
        y = bell_diam,
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
      name = "Bell diameter (mm)",
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
      "_medusae_bell_diam.png"
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

scale_factor <- max(plot_data_m$bell_diam, na.rm = TRUE) / 40

# mean bell diam values per treatment
mean_data_m <- plot_data_m %>%
  group_by(salinity_level, Day) %>%
  summarise(
    mean_diam = mean(bell_diam, na.rm = TRUE),
    se = sd(bell_diam, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot(
  plot_data_m,
  aes(
    Day,
    bell_diam,
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
    name = "Bell diameter (mm)",
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
    y = "Bell diameter (mm)",
    colour = "Replicate"
  ) +
  
  theme_bw()




p <- ggplot(
  plot_data_m,
  aes(
    x = Day,
    y = bell_diam,
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
    data = mean_data_m,
    aes(
      x = Day,
      ymin = mean_diam - se,
      ymax = mean_diam + se,
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
    data = mean_data_m,
    aes(
      x = Day,
      y = mean_diam,
      group = 1
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.5
  ) +
  
  geom_point(
    data = mean_data_m,
    aes(
      x = Day,
      y = mean_diam
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
    y = "Bell diameter (mm)"
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
    "output/med_bell_diam_results_by_treatment.png",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)


# plot the medusae YII for the individual replicates 

replicates_m <- plot_data_m %>%
  distinct(salinity_level, replicate)

for(i in seq_len(nrow(replicates_m))) {
  
  treatment_i <- replicates_m$salinity_level[i]
  replicate_i <- replicates_m$replicate[i]
  
  dat <- plot_data_m %>%
    filter(
      salinity_level == treatment_i,
      replicate == replicate_i
    )
  
  
  scale_factor <- max(dat$YII, na.rm = TRUE) / 40
  
  print(
    dat %>%
      select(Day, salinity_level, replicate, YII)
  )
  
  p <- ggplot(dat, aes(x = Day)) +
    
    geom_line(
      aes(y = YII),
      colour = "#0072B2",
      linewidth = 1
    ) +
    
    geom_point(
      aes(y = YII),
      colour = "#0072B2",
      size = 2
    ) +
    
    geom_step(
      aes(y = Salinity * scale_factor),
      colour = "#D55E00",
      linewidth = 1
    ) +
    
    scale_y_continuous(
      name = "YII",
      
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
      "_YII_rep",
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


for(treatment_i in unique(plot_data_m$salinity_level)) {
  
  dat <- plot_data_m %>%
    filter(salinity_level == treatment_i) %>%
    arrange(replicate, Day)
  
  scale_factor <- max(dat$YII, na.rm = TRUE) / 40
  
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
        y = YII,
        colour = factor(replicate),
        group = replicate
      ),
      linewidth = 1
    ) +
    
    geom_point(
      aes(
        y = YII,
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
      name = "YII",
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
      "_medusae_YII.png"
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

scale_factor <- max(plot_data_m$YII, na.rm = TRUE) / 40

# mean bell diam values per treatment
mean_data_m <- plot_data_m %>%
  group_by(salinity_level, Day) %>%
  summarise(
    mean_diam = mean(bell_diam, na.rm = TRUE),
    se_bd = sd(bell_diam, na.rm = TRUE) / sqrt(n()),
    mean_YII = mean(YII, na.rm = TRUE),
    se_YII = sd(YII, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot(
  plot_data_m,
  aes(
    Day,
    YII,
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
    name = "YII",
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
    y = "YII",
    colour = "Replicate"
  ) +
  
  theme_bw()




p <- ggplot(
  plot_data_m,
  aes(
    x = Day,
    y = YII,
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
    data = mean_data_m,
    aes(
      x = Day,
      ymin = mean_YII - se_YII,
      ymax = mean_YII + se_YII,
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
    data = mean_data_m,
    aes(
      x = Day,
      y = mean_YII,
      group = 1
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.5
  ) +
  
  geom_point(
    data = mean_data_m,
    aes(
      x = Day,
      y = mean_YII
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
    y = "Effective Quantum Yield (YII)"
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
    "output/med_YII_results_by_treatment.png",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)




