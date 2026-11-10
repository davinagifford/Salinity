# data_wrangle.R
#
# A file to manipulate data files and create key variables. 
#
# Created: 13 October 2026
# Created by: Davina Gifford
# Updated: 10 November 2026
# Updated by: Davina Gifford

# load libraries
library(tidyr)
library(tidyverse)
library(lubridate)

# import data
raw_data <- read.csv("data/Salinity_experiment_data.csv")
treatments <- read.csv("data/treatment_data.csv")

# rename columns
raw_data <- raw_data %>%
  rename(life_stage = Life.stage)

# Extract polyp data
polyps <- raw_data %>%
  filter(life_stage == "Polyp")

# Extract medusa data
medusae <- raw_data %>%
  filter(life_stage == "Medusa")

# tidy polyp data
polyps <- polyps %>%
  select(where(~ !all(is.na(.) | trimws(as.character(.)) == ""))) %>%
  rename(salinity_level = Salinity.level., polyp_number = Number.of.polyps.present, buds = How.many.polyps.with.buds., dead = Number.of.dead.polyps, strobilation = Stobilation.) 
  
head(polyps)


# tidy medusa data
medusae <- medusae %>%
  select(where(~ !all(is.na(.) | trimws(as.character(.)) == ""))) %>%
  rename(salinity_level = Salinity.level.2, bell_diam = Bell.diameter..mm., YII = YII.value, replicate = Replicate2 ) 

head(medusae)

medusae <- medusae %>%
  mutate(
    Date = as.Date(mdy_hms(Start.time)),
    Day = as.numeric(Date - min(Date, na.rm = TRUE)) + 1
  )

medusae %>%
  arrange(salinity_level, replicate, Day) %>%
  select(
    Day, Date, salinity_level, replicate,
    Alive., bell_diam, YII, Any.notes
  ) %>%
  head(30)

# Because medusae where measured twice, once for YII and once for daim, need to combine the data into a single row for each replicate

medusae_combined <- medusae %>%
  group_by(Date, Day, salinity_level, replicate) %>%
  summarise(
    bell_diam = first(na.omit(bell_diam)),
    YII = first(na.omit(YII)),
    Alive = first(na.omit(Alive.)),
    Changes_in_colour = first(na.omit(Changes.in.colour.)),
    Colour_description = first(na.omit(If.changes.in.colour..describe.)),
    Any_notes = paste(unique(na.omit(Any.notes)), collapse = "; "),
    .groups = "drop"
  )

