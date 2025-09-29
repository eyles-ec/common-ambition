library(dplyr)
library(tidyverse)
library(data.table)

#function for creating codelist flags easily. you can name the new column anything

flag_cl <- function(df, code_list, new_col_name = "has_any_code") {
  new_col <- sym(new_col_name) #language wrangling
  df %>%
    mutate(
      !!new_col := rowSums(across(starts_with("DiagCode"), ~ . %in% code_list)) > 0 #checks across all rows of DiagCode etc to search for the codelist.
    )
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

#set working directory

setwd("Own path")
      
#load datasets

unity_cab1 <- read.csv("./subdirectory/unity1.csv")
unity_cab2 <- read.csv("./subdirectory/unity2.csv")
ck_cab1 <- read.csv("./subdirectory/ck.csv")
croydon_cab1 <- read.csv("./subdirectory/croydon1.csv")
croydon_cab2 <- read.csv("./subdirectory/croydon2.csv")

#append datasets using append_data function

unity <- append_dfs(unity_cab1, ck_cab1, unity_cab2)

#croydon 2 onto croydon 1

croydon <- append_dfs(croydon_cab1, croydon_cab2)

#load shappt reference codes for HIV and PrEP

shappt_filter <- read.csv("./subdirectory/codelist/hiv_shappt.csv")

#pull out diagnostic codes to filter data
diag_codes<-shappt_filter$code

#filter only diagnostic columns with HIV-related codes
#there are ten separate diagnostic columns 
unity_cab1_hiv <- unity_cab1 %>%
  filter(if_any(starts_with("DiagCode"), ~ . %in% diag_codes))

#Generate simplified ethnicity column
unity_cab1_hiv <- unity_cab1_hiv %>%
  mutate(
    ethn_simple = case_when(
      EthnicGroupDescription %in% c("Caribbean",
                                    "White and Black Caribbean",
                                    "African",
                                    "Any other Black background",
                                    "White and Black African") ~ EthnicGroupDescription,
      TRUE ~ "non ACHC"
    )
  )

#generate HIV test variables

#variable if test declined or inappropriate (Shappt P1B, P1C)
cl_hiv_dec = c("P1B", "P1C")

unity_cab1_hiv <- flag_cl(unity_cab1_hiv, cl_hiv_dec, "declined_hiv_test")

#variable if test conducted (Shappt P1A, T4, T7, T-HIV)
cl_hiv_test = c("P1A", "T4", "T7", "T-HIV")

unity_cab1_hiv <- flag_cl(unity_cab1_hiv, cl_hiv_test, "hiv_test")

#generate HIV diagnosis variables

#extant HIV (Shappt H, H1X H1AX, H1BX (last three are previoulsy diagnosis elsewhere diagnostic codes))
cl_extant = c("H", "H1X", "H1AX", "H1BX")

unity_cab1_hiv <- flag_cl(unity_cab1_hiv, cl_extant, "extant_hiv")

#New diagnosis (Shappt H1, H1A, H1B)
cl_new = c("H1", "H1A", "H1B")

unity_cab1_hiv <- flag_cl(unity_cab1_hiv, cl_new, "new_hiv")

#New diagnosis late (Shappt H1B)

unity_cab1_hiv <- flag_cl(unity_cab1_hiv, "H1B", "new_late_hiv")

#generate PrEP variables

#Starting or continuing PrEP even through other source, including prescriptions (Shappt O41, O42, O43, O51, O52, O53)
cl_curr_prep <- c("O41", "O42", "O43", "O51", "O52", "O53")

unity_cab1_hiv <- flag_cl(unity_cab1_hiv, cl_curr_prep, "current_prep")

#Declined PrEP (Shappt O44)

unity_cab1_hiv <- flag_cl(unity_cab1_hiv, "O44", "declined_prep")

#Stopped PrEP (Shappt O45)

unity_cab1_hiv <- flag_cl(unity_cab1_hiv, "O45", "stopped_prep")

#generate year, month, and (ISO) week variables
unity_cab1_hiv$EventDate <- as.Date(unity_cab1_hiv$EventDate, format = "%d/%m/%Y")

unity_cab1_hiv <- unity_cab1_hiv %>%
  mutate(
    year = year(EventDate),
    month = month(EventDate),
    week = isoweek(EventDate)  # ISO week number (1–53)
  )

  