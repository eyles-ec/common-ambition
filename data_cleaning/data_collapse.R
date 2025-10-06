library(dplyr)
library(tidyverse)
library(data.table)
library(rlang)
library(lubridate)
library(collapse)

#episode count function, including mitigating for repeated episode numbers on different dates

episode_count <- function(df, patient = "PatientIdentifier", episode = "EpisodeNumber", date = "EventDate") {
  # Create a summary table of unique episode-date combinations per person, and count them
  episode_summary <- df %>%
    distinct(.data[[patient]], .data[[episode]], .data[[date]]) %>%
    group_by(.data[[patient]]) %>%
    summarise(episode_count = n(), .groups = "drop")
  
  # Join the count back to original dataframe with a left join
  df <- df %>%
    left_join(episode_summary, by = patient)
  
  return(df)
}

#collapse function

collapse_weekly <- function(df,
                            week_col = "week",
                            year_col = "year",
                            group_vars = c("ethn_simple", "location"),
                            sum_vars = c("hiv_test", "declined_hiv_test", "current_prep", "declined_prep", "stopped_prep")) {
  df %>%
    group_by(.data[[year_col]], .data[[week_col]], !!!syms(group_vars)) %>%
    summarise(across(all_of(sum_vars), ~ sum(.x == TRUE, na.rm = TRUE)), .groups = "drop")
}

#rate calculation function

add_rates <- function(df,
                      count_vars = c("hiv_test", "declined_hiv_test", "current_prep", "declined_prep", "stopped_prep"),
                      population_col = "population",
                      scale = 1000) {
  df %>%
    mutate(across(all_of(count_vars),
                  ~ (.x / .data[[population_col]]) * scale,
                  .names = "{.col}_rate"))
}


#set working directory

setwd("YOUR WD")

#load analysis dataset

cab <- read.csv("./subdirectory/Processed/combined_episodes.csv")

#ensure type matching for date column
cab$EventDate<- as.Date(cab$EventDate)

#generate episode count

cab <- episode_count(cab)

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

#calculate rates per 1000 (can specify scale = to change rate)

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
