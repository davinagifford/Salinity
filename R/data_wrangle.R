# data_wrangle.R
#
# A file to manipulate data files and create key variables. 
#
# Created: 13 October 2026
# Created by: Davina Gifford
# Updated: 13 November 2026
# Updated by: Davina Gifford

# load libraries
library(tidyr)
library(tidyverse)
library(lubridate)

# import data
raw_data <- read.csv("data/Salinity_experiment_data.csv")
treatments <- read.csv("data/treatment_data.csv")

# update treatment data 

treatments <- treatments %>% 
  mutate(Treatment = case_when(Treatment == "10" ~ "10 ppt",
                               Treatment == "20" ~ "20 ppt",
                               Treatment == "15" ~ "15 ppt",
                               Treatment == "5" ~ "5 ppt",
                               Treatment == "Control" ~ "Control",
                               Treatment == "Procedural control" ~ "Procedural control",
                               ))

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
  rename(salinity_level = Salinity.level., replicate = Replicate, polyp_number = Number.of.polyps.present, buds = How.many.polyps.with.buds., dead = Number.of.dead.polyps, strobilation = Stobilation.) 

polyps <- polyps %>% 
  mutate(
    Date = as.Date(mdy_hms(Start.time)),
    Day = as.numeric(Date - min(Date, na.rm = TRUE)) + 1,
    salinity_level = case_when(salinity_level == "Black" ~ "Procedural control",
                               salinity_level == "Blue" ~ "Control",
                               salinity_level == "White" ~ "5 ppt",
                               salinity_level == "Green" ~ "10 ppt",
                               salinity_level == "Red" ~ "15 ppt",
                               salinity_level == "Yellow" ~ "20 ppt",
                               TRUE ~ "10 ppt"), # a row was being changed to NA instead of green so this is to fix that.
    # fix a duplicate 
    replicate = if_else(
        Day == 4 &
        salinity_level == "Control" &
        replicate == "5" &
        polyp_number == 4,
      "3",
      replicate
    )
  ) 

# insert row of missing data
missing_dat <- polyps[18,]

missing_dat$replicate <- "1"
missing_dat$polyp_number <- 4

polyps <- rbind(polyps, missing_dat)

# update values in columns with NA to 0

polyps <- polyps %>% 
  replace_na(list(dead = 0)) %>% 
  replace_na(list(buds = 0))


polyps_ordered <- polyps %>% 
  filter(replicate != "6 (spare)") %>%
  mutate(replicate = as.numeric(replicate)) %>%
  mutate(phase = case_when(Day <= 2 ~ "Acclimation",
                           Day <= 6 ~ "Exposure",
                           Day <= 10 ~ "Recovery",
                           Day < 16 ~ "Post-Recovery")) %>% 
  arrange(Date, Day, salinity_level, replicate) %>%
  select(
    Date, Day, phase, salinity_level, replicate,
    polyp_number, buds, dead, strobilation, Any.notes.
  ) 
  
head(polyps_ordered)
summary(polyps_ordered)

# tidy medusa data
medusae <- medusae %>%
  select(where(~ !all(is.na(.) | trimws(as.character(.)) == ""))) %>%
  rename(salinity_level = Salinity.level.2, bell_diam = Bell.diameter..mm., YII = YII.value, replicate = Replicate2 ) 

head(medusae)

medusae <- medusae %>%
  mutate(
    Date = as.Date(mdy_hms(Start.time)),
    Day = as.numeric(Date - min(Date, na.rm = TRUE)) + 1,
    salinity_level = case_when(salinity_level == "Black" ~ "Procedural control",
                               salinity_level == "Blue" ~ "Control",
                               salinity_level == "White" ~ "5 ppt",
                               salinity_level == "Green" ~ "10 ppt",
                               salinity_level == "Red" ~ "15 ppt",
                               salinity_level == "Yellow" ~ "20 ppt"),
    
  ) 

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

# save tidied data for ease of use
saveRDS(medusae_combined, file = "var/medusae_tidy.rds")
saveRDS(polyps_ordered, file = "var/polyp_tidy.rds")
saveRDS(treatments, file = "var/treatment_tidy.rds")

