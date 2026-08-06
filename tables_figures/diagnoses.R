library(dplyr)
library(tidyr)

#function to pull test episodes and diagnosis episodes from episodic data
#use test and diagnostic codes to check subsequent diagnosis after test
#as this is episodic, if someone has more than one test/diagnosis pair they will be counted however
#many times they have it
test_to_diagnosis <- function(cab,
                             test_codes,
                             diagnosis_codes,
                             output_file = NULL) {
  
  diag_cols <- paste0("DiagCode", 1:10)
  
  #pull out matching test episodes
  test_eps <- cab %>%
    filter(
      if_any(all_of(diag_cols), ~ .x %in% test_codes)
    )
  
  #pull out matching diagnosis episodes, but keep just minimum needed for counts
  diag_eps <- cab %>%
    filter(
      if_any(all_of(diag_cols), ~ .x %in% diagnosis_codes)
    ) %>%
    select(
      PatientIdentifier,
      EventDate
    ) %>%
    rename(diag_date = EventDate)
  
  #Summarise by group_bristol (aka ACHC non ACHC and location)
  results <- test_eps %>%
    select(
      PatientIdentifier,
      EventDate,
      group_bristol,
      time
    ) %>%
    rename(test_date = EventDate) %>%
    mutate(test_id = row_number()) %>%
    left_join(
      diag_eps,
      by = "PatientIdentifier",
      relationship = "many-to-many"
    ) %>%
    group_by(test_id) %>%
    summarise(
      group_bristol = first(group_bristol),
      period = ifelse(first(time) <= 0, "Pre-CAB", "CAB"),
      diagnosed_later = any(diag_date >= first(test_date), na.rm = TRUE), #if diagnosis happens after a test
      .groups = "drop"
    ) %>%
    group_by(group_bristol, period) %>%
    summarise(
      n_tests = n(),
      n_diag = sum(diagnosed_later),
      pct_diag = round(100 * mean(diagnosed_later), 1),
      .groups = "drop"
    )
  
  if (!is.null(output_file)) {
    write.csv(results, output_file, row.names = FALSE)
  }
  
  return(results)
}


#collect paths - save your own paths to a file called paths.R that is ignored by git (.gitignore) 

source("../paths.R")

#set working directory 

setwd(wd)

#load clean dataset at episode level
cab<- read.csv("./Analysis/Processed/combined_episodes_corrected.csv")

#remove croydon non ACHC
cab <- cab %>% 
  filter(group_bristol != "Croydon non ACHC")

#load codelists
hiv_filter <- read.csv("./Analysis/codelists/hiv_shappt.csv")
sti_filter <- read.csv("./Analysis/codelists/sti_shappt.csv")
diag_filter <- read.csv("./Analysis/codelists/diagnosis_shappt.csv")

#pull out hiv diagnosis
hiv_diagnosis<-c("H", "H1X", "H1AX", "H1BX",
  "H1", "H1A", "H1B"
)

#remove diagnosis codes from hiv filter to get just tests 
hiv_filter <- hiv_filter %>%
  filter(!code %in% c(
    "H", "H1X", "H1AX", "H1BX",
    "H1", "H1A", "H1B"
  ))

#create just STI diagnosis
sti_diagnosis_tmp <- diag_filter %>%
  filter(!code %in% c(
    "H", "H1X", "H1AX", "H1BX",
    "H1", "H1A", "H1B"
  ))

#diagnosis to vector
sti_diagnosis <- sti_diagnosis_tmp$code

#remove tmp file 
rm(sti_diagnosis_tmp)

#create test codelist
test_codes <- c(hiv_filter$code, sti_filter$code)

#hiv test codelist
hiv_test <- hiv_filter$code

#sti test codelist
sti_test <- sti_filter$code

results_overall <- test_to_diagnosis(cab = cab,
                                    test_codes = test_codes,
                                    diagnosis_codes = diag_filter$code,
                                    output_file = "./Analysis/all_tests_to_diagnosis.csv"
                                    )

results_sti <- test_to_diagnosis(cab = cab,
                                     test_codes = sti_test,
                                     diagnosis_codes = sti_diagnosis$code,
                                     output_file = "./Analysis/sti_tests_to_diagnosis.csv"
                                )

results_hiv <- test_to_diagnosis(cab = cab,
                                 test_codes = hiv_test,
                                 diagnosis_codes = hiv_diagnosis,
                                 output_file = "./Analysis/hiv_tests_to_diagnosis.csv"
                                )



