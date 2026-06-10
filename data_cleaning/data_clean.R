library(dplyr)
library(tidyverse)
library(data.table)
library(rlang)
library(lubridate)

#function for creating codelist flags easily. you can name the new column anything

flag_cl <- function(df, code_list, new_col_name = "has_any_code") {
  new_col <- sym(new_col_name) #language wrangling
  df %>%
    mutate(
      !!new_col := rowSums(across(starts_with("DiagCode"), ~ . %in% code_list)) > 0 #checks across all rows of DiagCode etc to search for the codelist.
    )
}

#function for checking and consolidating duplicates as new croydon extract has overlapping dates with cab1
dedup <- function(df, key_vars){

    dup_temp<- df %>% 
                  group_by(across(all_of(key_vars))) %>%
                  summarise(n = n(), .groups = "drop")
          
    #sum up all the T in this
    n_dup <- sum(dup_temp$n > 1)
          #write a message to the console
          message("There are ", n_dup, " duplicate episodes")
          #return original data with . 
    
    df<- df %>%
      #group by variables, keep only the first of any duplicates
      group_by(across(all_of(key_vars))) %>%
      slice(1) %>%
      ungroup() #return data in episodic format
    
  return(df)
  
}

#function for cleaning the data and creating filters
data_processing <- function(df, hiv_codes, sti_codes) {
  #add a flag for whether an episode includes hiv diagnosis, hiv testing and/or prep or sti testing
  df <- df %>%
    mutate(
      is_hiv = if_any(starts_with("DiagCode"), ~ . %in% hiv_codes),
      is_sti = if_any(starts_with("DiagCode"), ~ . %in% sti_codes),
      is_sti_hiv = is_sti | is_hiv
    ) 
  
  #Simplified ethnicity
  #This includes two and one digit ethnicity codes 
  df <- df %>%
    mutate(
      ethn_simple = case_when(
        EthnicGroupCode %in% c("M", "D", "N", "P") ~ "ACHC",
        EthnicGroupCode %in% c("MO", "NO", "PD", "DO", "PO", "EO") ~ "ACHC",
        TRUE ~ "non ACHC"
      )
    )
  
  # HIV test flags
  df <- flag_cl(df, c("P1B", "P1C"), "declined_hiv_test")
  df <- flag_cl(df, c("P1A", "T4", "T7", "T-HIV", "S2", "RHTR", "HA12", "RHTN"), "hiv_test")
  
  # HIV diagnosis flags
  df <- flag_cl(df, c("H", "H1X", "H1AX", "H1BX"), "extant_hiv")
  df <- flag_cl(df, c("H1", "H1A", "H1B"), "new_hiv")
  df <- flag_cl(df, "H1B", "new_late_hiv")
  
  # PrEP flags
  df <- flag_cl(df, c("O41", "O42", "O43", "O51", "O52", "O53", "O54", "PREP", "PREP3", "PREP4", "PREP1", "PREP2"), "current_prep")
  df <- flag_cl(df, c("O44", "PREP0"), "declined_prep")
  df <- flag_cl(df, "O45", "stopped_prep")
  df <- flag_cl(df, c("HA13", "PEPS"), "pep")
  
  #STI flag with HIV, plus count
  df <- flag_cl(df, c("T1", "T12", "T2", "T3", "T5", "T6", "TT", "P1A", "T4", "T7", "T-HIV", "S2", "RHTR", "HA12", "RHTN"), "sti_test_hiv")
  df <- df %>%
    mutate(
      sti_test_count_hiv = rowSums(across(starts_with("DiagCode"), ~ . %in% c("T1", "T12", "T2", "T3", "T5", "T6", "TT", "P1A", "T4", "T7", "T-HIV", "S2", "RHTR", "HA12", "RHTN")))
    )
  
  # STI flag no HIV, plus count
  df <- flag_cl(df, c("T1", "T12", "T2", "T3", "T5", "T6", "TT"), "sti_test_no_hiv")
  df <- df %>%
    mutate(
      sti_test_count_no_hiv = rowSums(across(starts_with("DiagCode"), ~ . %in% c("T1", "T12", "T2", "T3", "T5", "T6", "TT")))
    )
  
  # Date breakdown
  df$EventDate <- as.Date(df$EventDate, format = "%d/%m/%Y")
  df <- df %>%
    mutate(
      #break down components
      year = year(EventDate),
      month = month(EventDate),
      week = isoweek(EventDate),
      quarter = quarter(EventDate),
      #create pre/post intervention timeline
      time = as.integer(difftime(EventDate, as.Date("2022-04-25"), units = "weeks")),
      period = if_else(time > 0, 1, 0)
    )
  
  return(df)
}

#episode count function, including mitigating for repeated episode numbers on different dates

episode_count <- function(df, 
                          patient = "PatientIdentifier", 
                          episode = "EpisodeNumber", 
                          date = "EventDate", 
                          time_col = "time") {
  # Count total episodes per patient
  episode_summary <- df %>%
    distinct(.data[[patient]], .data[[episode]], .data[[date]]) %>%
    group_by(.data[[patient]]) %>%
    summarise(episode_count = n(), .groups = "drop")
  
  # Count weekly episodes per patient
  weekly_summary <- df %>%
    distinct(.data[[patient]], .data[[episode]], .data[[date]], .data[[time_col]]) %>%
    group_by(.data[[patient]], .data[[time_col]]) %>%
    summarise(weekly_episode_count = n(), .groups = "drop")
  
  # Join both summaries back to original dataframe
  df <- df %>%
    left_join(episode_summary, by = patient) %>%
    left_join(weekly_summary, by = c(patient, time_col))
  
  return(df)
}

#function for combining dataframes. can be more than one
#change master_index = 1 to change which is the 'master' df  in the list
#this includes type coercion to the master df, in case of type mismatch
append_dfs <- function(..., master_index = 1) {
  dfs <- list(...)
  master <- dfs[[master_index]]
  master_cols <- names(master)
  master_types <- sapply(master, class)
  
  aligned_dfs <- lapply(seq_along(dfs), function(i) {
    df <- dfs[[i]]
    
    # Add missing columns
    missing_cols <- setdiff(master_cols, names(df))
    if (length(missing_cols) > 0) {
      df[missing_cols] <- NA
      message("Added missing columns to df", i, ": ", paste(missing_cols, collapse = ", "))
    }
    
    # Reorder columns to match master
    df <- df[, master_cols, drop = FALSE]
    
    # Coerce types to match master
    for (col in master_cols) {
      if (!inherits(df[[col]], master_types[[col]])) {
        suppressWarnings({
          df[[col]] <- tryCatch({
            as(df[[col]], master_types[[col]])
          }, error = function(e) {
            warning("Could not coerce column '", col, "' in df", i, " to type ", master_types[[col]])
            df[[col]]
          })
        })
      }
    }
    
    return(df)
  })
  
  bind_rows(aligned_dfs)
}


#collect paths - save your own paths to a file called paths.R that is ignored by git (.gitignore) 

source("../paths.R")

#set working directory 

setwd(wd)

#load datasets

unity_cab1 <- read.csv("./RawData/unity.csv")
croydon_cab1 <- read.csv("./RawData/croydon.csv")
ck_cab1 <- read.csv("./RawData/CK_csv.csv")

unity_cab2 <- read.csv("./RawData_CAB2/unity.csv")
croydon_cab2_p1 <- read.csv("./RawData_CAB2/croydon_new_extract/Croydon_oct23_dec23.csv")
croydon_cab2_p2 <- read.csv("./RawData_CAB2/croydon_new_extract/Croydon_jan24_june24.csv")
croydon_cab2_p3 <- read.csv("./RawData_CAB2/croydon_new_extract/Croydon_july24_dec24.csv")
croydon_cab2 <- append_dfs(croydon_cab2_p1, croydon_cab2_p2, croydon_cab2_p3)

#date admin for corydon
croydon_cab1 <- croydon_cab1 %>%
  mutate(EventDate = as.Date(EventDate, format = "%d/%m/%Y"))

croydon_cab1<- croydon_cab1 %>%  filter(EventDate <= as.Date("2023-10-27")) #ensure there's no duplicates between cab1 and cab2 extracts 

croydon_cab2 <- croydon_cab2 %>%
  mutate(EventDate = as.Date(EventDate, format = "%d/%m/%Y"))


#append datasets using append_data function
unity <- append_dfs(unity_cab1, ck_cab1, unity_cab2)

#date admin
unity <- unity %>% mutate(EventDate = as.Date(EventDate, format = "%d/%m/%Y"))

#croydon 2 onto croydon 1

croydon <- append_dfs(croydon_cab1, croydon_cab2)

croydon <- croydon %>% mutate(EventDate = as.Date(EventDate, format = "%d/%m/%Y"))

#set key vars for duplicate checking
key_vars <- c("PatientIdentifier", "EventDate")

#check/remove any duplicates 

unity <- dedup(unity, key_vars)
croydon <- dedup(croydon, key_vars) 

#load shappt reference codes for HIV and PrEP

hiv_filter <- read.csv("./Analysis/codelists/hiv_shappt.csv")
sti_filter <- read.csv("./Analysis/codelists/sti_shappt.csv")

#use data_processing function to process croydon and unity datasets
unity_hiv <- data_processing(unity, hiv_filter$code, sti_filter$code)
croydon_hiv <- data_processing(croydon, hiv_filter$code, sti_filter$code)

#add episode count
unity_hiv <- episode_count(unity_hiv)
croydon_hiv <- episode_count(croydon_hiv)

#add croydon or bristol flag
unity_hiv$location <- "Bristol"
croydon_hiv$location <- "Croydon"

#sort by date
unity_hiv <- unity_hiv %>% arrange(EventDate)
croydon_hiv <- croydon2_hiv %>% arrange(EventDate)

#save intermediate step as CSV
write.csv(unity_hiv, "./Analysis/Processed/unity_episodes.csv")
write.csv(croydon_hiv, "./Analysis/Processed/croydon_episodes_NEWEST.csv")

#append datasets together to create one analytic dataset

combined_data <- append_dfs(unity_hiv, croydon_hiv)

combined_data <- combined_data %>% arrange(EventDate)

#add in groupings
combined_data <- combined_data %>%
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

#save combined dataset
write.csv(combined_data, "./Analysis/Processed/combined_episodes_corrected.csv")

