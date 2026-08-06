

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

#remove diagnosis codes from hiv filter to get just tests 
hiv_filter <- hiv_filter %>%
  filter(!code %in% c(
    "H", "H1X", "H1AX", "H1BX",
    "H1", "H1A", "H1B"
  ))

#create test codelist
test_codes <- c(hiv_filter, sti_filter)

#get names for diagnosis/coded columns 
diag_cols <- paste0("DiagCode", 1:10)

#check for test codes, pull out the episodes 
test_eps <- cab %>%
filter(
  if_any(all_of(diag_cols), ~ .x %in% test_codes)
)

#check for diagnostic codes, pull out the episodes 
diag_eps <- cab %>%
  filter(
    if_any(all_of(diag_cols), ~ .x %in% diag_filter)
  )