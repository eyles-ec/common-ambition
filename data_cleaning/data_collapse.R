library(dplyr)
library(tidyverse)
library(data.table)
library(rlang)
library(lubridate)
library(collapse)

#collapse function

collapse_weekly <- function(df,
                            time_col = "time",
                            year_col = "year",
                            period_col = "period",
                            group_vars = c("ethn_simple", "location"),
                            sum_vars = c("hiv_test", "declined_hiv_test", "current_prep", "declined_prep", "stopped_prep", "extant_hiv", "new_hiv", "new_late_hiv", "sti_test", "pep"),
                            count_vars = c("weekly_episode_count", "sti_test_count")) {
  df %>%
    group_by(.data[[year_col]], 
             .data[[time_col]], 
             .data[[period_col]], 
             !!!syms(group_vars)) %>%
    summarise(
      across(all_of(sum_vars), ~ sum(.x == TRUE, na.rm = TRUE)),
      across(all_of(count_vars), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}

#rate calculation function

add_rates <- function(df,
                      count_vars = c("hiv_test", "declined_hiv_test", "current_prep", "declined_prep", "stopped_prep", "extant_hiv", "new_hiv", "new_late_hiv", "sti_test", "weekly_episode_count", "sti_test_count"),
                      population_col = "population",
                      scale = 1000) {
  df %>%
    mutate(across(all_of(count_vars),
                  ~ (.x / .data[[population_col]]) * scale,
                  .names = "{.col}_rate"))
}

#set working directory

setwd(YOURWD)

#load analysis dataset

cab <- read.csv("./subdirectory/Processed/combined_episodes.csv")

#ensure type matching for date column
cab$EventDate<- as.Date(cab$EventDate)

#collapse to weekly count data

cab_weekly <- collapse_weekly(cab)

#generate population offsets 
#first build a reference table for the population of each
pop_table <- tibble(
  location = c("Bristol", "Croydon"),
  total_pop = c(479000, 390800)
)

#then ACHC/non ACHC proportions for each 
group_props <- tibble(
  ethn_simple = c("ACHC", "non ACHC", "ACHC", "non ACHC"),
  location = c("Bristol", "Bristol", "Croydon", "Croydon"),
  proportion = c(0.08, 0.92, 0.264, 0.736)
)

#calculate population of each group
group_pop <- group_props %>%
  left_join(pop_table, by = "location") %>%
  mutate(population = proportion * total_pop,
         offset = log(population))

#housekeeping
rm(pop_table, group_props)

#join in population and offset to weekly table

cab_weekly <- cab_weekly %>%
  left_join(group_pop, by = c("ethn_simple", "location"))

#calculate rates per 1000 (can specify scale = to change rate). this is by the group population

cab_weekly <- add_rates(cab_weekly)

#generate grouping variables for modelling 

cab_weekly <- cab_weekly %>%
  mutate(
    group_bristol = case_when(
      location == "Bristol" & ethn_simple == "ACHC" ~ "Bristol ACHC",
      location == "Bristol" & ethn_simple == "non ACHC" ~ "Bristol non ACHC",
      location == "Croydon" & ethn_simple == "ACHC" ~ "Croydon ACHC",
      location == "Croydon" & ethn_simple == "non ACHC" ~ "Croydon non ACHC",
      TRUE ~ NA_character_
    ),
    group_bristol = factor(group_bristol, levels = c(
      "Bristol ACHC", "Bristol non ACHC", "Croydon ACHC", "Croydon non ACHC"
    ))
  )

#export weekly data 
write.csv(cab_weekly, "./subdirectory/Analysis/weekly_combined.csv")
