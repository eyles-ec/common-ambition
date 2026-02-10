#This is to check for matching descriptions and their corresponding codes as the data are messy
#and codes don't correspond to the national SHAPPT codelists

library(dplyr)
library(tidyr)
library(stringr)

check_matches <- function(df, include_keywords, exclude_keywords = NULL, flag_name = "keyword_flag") {
  include_pattern <- paste(include_keywords, collapse = "|")
  
  # Only build exclusion pattern if exclusion keywords are provided
  exclude_pattern <- if (!is.null(exclude_keywords)) paste(exclude_keywords, collapse = "|") else NULL
  
  df <- df %>%
    mutate(
      !!flag_name := apply(
        select(., starts_with("DiagDescription")),
        1,
        function(row) {
          has_include <- any(stringr::str_detect(row, regex(include_pattern, ignore_case = TRUE)))
          has_exclude <- if (!is.null(exclude_pattern)) {
            any(stringr::str_detect(row, regex(exclude_pattern, ignore_case = TRUE)))
          } else {
            FALSE
          }
          has_include & !has_exclude
        }
      )
    )
  
  #pull out rows that match
  matched_rows <- df %>% filter(!!sym(flag_name) == TRUE)
  
  #count rows that match
  match_count <- nrow(matched_rows)
  
  #pivot longer to match codes and descriptions
  long_df <- matched_rows %>%
    select(starts_with("DiagCode"), starts_with("DiagDescription")) %>%
    pivot_longer(
      cols = everything(),
      names_to = c(".value", "index"),
      names_pattern = "(DiagCode|DiagDescription)(\\d+)"
    ) %>%
    filter(stringr::str_detect(DiagDescription, regex(include_pattern, ignore_case = TRUE)))
  
  #return a list with matched dataframe and count of matches
  return(list(matches = long_df, count = match_count))
}


#collect paths - save your own paths to a file called paths.R that is ignored by git (.gitignore) 

source("../paths.R")

#set working directory 

setwd(wd)

cab <- read.csv("./Processed/combined_episodes.csv")

# list your included keywords, can be partial matches
# exclusion criteria are also allowed
# below example is looking for hiv but excluding tests and screening

include <- c("hiv")
exclude <- c("test", "scr")

#call the function - exclusion criteria are not mandatory

matches <- check_matches(cab, include, exclude)


#print out matches in the data to check n

matches$matches %>%
  distinct(DiagCode, DiagDescription) %>%
  print(n = Inf)
