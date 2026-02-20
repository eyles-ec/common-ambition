library(dplyr)
library(tidyr)
library(gt) 

#helper function to redact counts of (by default, less than 6, but modifiable with threshold argument)
#call it in generate_table1 function after categorical summary
redact_small_counts <- function(df, threshold = 6) {
  df %>%
    mutate(across(
      where(is.character),
      ~ ifelse(
        !is.na(value) &
          # redact only 1–5, not 0
          grepl(paste0("^([1-", threshold - 1, "])\\s*\\("), .),
        paste0("<", threshold),
        .
      )
    ))
}

#generate patient level data from episodic for input into table 1 function
make_patient_data <- function(data) {
  data %>%
    group_by(PatientIdentifier) %>%
    summarise(
      #demographics
      Sex = first(Sex),
      AgeAtAttendance = mean(AgeAtAttendance, na.rm = TRUE),
      
      #create counts from logical variables (e.g. how many HIV tests/person)
      across(
        .cols = c(hiv_test, declined_hiv_test,
                  sti_test_no_hiv, sti_test_hiv),
        .fns = ~ sum(. == TRUE, na.rm = TRUE),
        .names = "{.col}"
      ),
      
      #keep HIV status as categorical (rare outcome)
      extant_hiv     = any(extant_hiv == TRUE),
      new_hiv        = any(new_hiv == TRUE),
      new_late_hiv   = any(new_late_hiv == TRUE),
      
      #STI test counts to summaries
      sti_test_count_no_hiv = sum(sti_test_count_no_hiv, na.rm = TRUE),
      sti_test_count_hiv    = sum(sti_test_count_hiv, na.rm = TRUE),
      
      #PrEP treated the same as HIV
      current_prep  = any(current_prep == TRUE),
      declined_prep = any(declined_prep == TRUE),
      stopped_prep  = any(stopped_prep == TRUE),
      
      #summary of service use, e.g. episodes per patient
      episode_count = n(),
      
      #assign to the right group
      group_bristol = first(group_bristol),
      
      .groups = "drop" #ungroup output
    )
}

#function to generate table 1 (Descriptives) for any input data. Requires vector of which 
#variables are categorical, continuous, labels for the variables for the table, a 'row grouping,'
#e.g. ordering by topic, and a title (defaults to table 1)
generate_table1 <- function(data, group_var, categorical_vars, continuous_vars,
                            variable_labels, row_groups, title = "Table 1: Descriptive Summary by Group") {
  
  #categorical summary 
  cat_summary <- data %>%
    mutate(across(all_of(categorical_vars), as.factor)) %>% #factorise in case they're not
    dplyr::select(all_of(c(group_var, categorical_vars))) %>% #subset only cat variables
    pivot_longer( #long table is easier to count
      cols = -all_of(group_var),
      names_to = "variable",
      values_to = "value"
    ) %>%
    group_by(variable, value, .data[[group_var]]) %>%  
    summarise(n = n(), .groups = "drop") %>% #produce count
    group_by(variable, .data[[group_var]]) %>%
    mutate(percent = round(100 * n / sum(n), 1)) %>% #produce %
    ungroup() %>%
    mutate(label = paste0(n, " (", percent, "%)")) %>% #create a label
    dplyr::select(any_of(c("variable", "value", group_var, "label"))) %>% #select only what we need to keep
    pivot_wider(names_from = all_of(group_var), values_from = label) #make it wide again
  
  #redact any count <6 (use threshold argument if you want a different one)
  cat_summary <- redact_small_counts(cat_summary)
  
  #continuous summary
  cont_summary <- data %>%
    dplyr::select(all_of(c(group_var, continuous_vars))) %>% #keep only desired variables
    pivot_longer( #long format means it's easier to summarise
      cols = -all_of(group_var),
      names_to = "variable",
      values_to = "value"
    ) %>%
    group_by(variable, .data[[group_var]]) %>%
    summarise(
      mean = round(mean(value, na.rm = TRUE), 1), #get means and sds
      sd   = round(sd(value, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    mutate(label = paste0(mean, " (", sd, ")")) %>% #create labels
    dplyr::select(any_of(c("variable", group_var, "label"))) %>% #keep only what we need for table
    pivot_wider(names_from = all_of(group_var), values_from = label) #make wide again
  
  #combine the two summaries together
  table1 <- bind_rows(cat_summary, cont_summary)
  
  #add in the variable labels to the table
  table1$variable <- variable_labels[table1$variable]
  
  #reorder the variables by your row_groups argument
  #can keep variables of the same topic together
  ordered_vars <- unlist(row_groups)
  ordered_labels <- variable_labels[ordered_vars]
  
  table1 <- table1 %>%
    mutate(variable = factor(variable, levels = ordered_labels)) %>%
    arrange(variable)
  
  #build a gt table, which is nicer to look at in R, exportable to LaTeX
  gt_tbl <- table1 %>%
    gt(rowname_col = "variable") %>%
    tab_header(title = title)
  
  #add the row groups to the gt table
  for (grp in names(row_groups)) {
    gt_tbl <- gt_tbl %>%
      tab_row_group(
        label = grp,
        rows = variable %in% variable_labels[row_groups[[grp]]]
      )
  }
  
  #returns the data exportable to csv, and the gt table (exportable to LaTeX etc)
  return(list(data = table1, gt = gt_tbl))
}

#pull wd from paths.R (put in .gitignore)
source("../paths.R")

#set working directory 

setwd(wd)

#load analysis dataset

cab <- read.csv("./Processed/combined_episodes.csv")

#create combined location x ethnicity 

cab <- cab %>%
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

#remove the non ACHC from Croydon as this group is not used
cab_subset <- cab[cab$group_bristol != "Croydon non ACHC",]

#define variables for table
group_var <- "group_bristol"

#episodic categorical variables
categorical_vars <- c(
  "Sex", "hiv_test", "declined_hiv_test", "extant_hiv",
  "new_hiv", "new_late_hiv", "current_prep", "declined_prep",
  "stopped_prep", "sti_test_no_hiv", "sti_test_hiv"
)

#episodic continuous variables
continuous_vars <- c("AgeAtAttendance", "episode_count", "sti_test_count_no_hiv", "sti_test_count_hiv")

#variable labels for table (nicer than variable names)
variable_labels <- c(
  #demographics
  Sex = "Sex",
  AgeAtAttendance = "Age (years)",
  
  #HIV and HIV testing
  hiv_test = "HIV test",
  declined_hiv_test = "Declined HIV test",
  extant_hiv = "Known HIV",
  new_hiv = "New HIV diagnosis",
  new_late_hiv = "Late HIV diagnosis",
  
  #PrEP
  current_prep = "Current PrEP use",
  declined_prep = "Declined PrEP",
  stopped_prep = "Stopped PrEP",
  
  #STI testing (with and without HIV)
  sti_test_no_hiv = "STI test (without HIV)",
  sti_test_hiv = "STI test (with HIV)",
  sti_test_count_no_hiv = "STI test count (without HIV)",
  sti_test_count_hiv = "STI test count (with HIV)",
  
  #service use
  episode_count = "Episode count",
  
  #grouping variable
  group_bristol = "Group"
)

#ordering list, so the table will be in this order
row_groups <- list(
  "Demographics" = c("Sex", "AgeAtAttendance"),
  "Testing" = c("hiv_test", "declined_hiv_test", "extant_hiv", "new_hiv", "new_late_hiv", "sti_test_no_hiv", "sti_test_hiv", "sti_test_count_no_hiv", "sti_test_count_hiv"),
  "PrEP" = c("current_prep", "declined_prep", "stopped_prep"),
  "Service use" = c("episode_count")
)

#patient level categorical variables 
categorical_vars_patient <- c(
  "Sex",
  "extant_hiv", "new_hiv", "new_late_hiv",
  "current_prep", "declined_prep", "stopped_prep"
)

#patient level continuous variables 
continuous_vars_patient <- c(
  "AgeAtAttendance",
  "episode_count",
  "hiv_test", "declined_hiv_test",
  "sti_test_count_no_hiv", "sti_test_count_hiv"
)

#subset pre and post periods
cab_pre  <- cab_subset %>% filter(time < 0)
cab_post <- cab_subset %>% filter(time >= 0)

#generate patient data for overall and pre/post period
patient_all  <- make_patient_data(cab_subset)
patient_pre  <- make_patient_data(cab_pre)
patient_post <- make_patient_data(cab_post)

#generate tables for episodic data
table1_all  <- generate_table1(cab_subset,  group_var, categorical_vars, continuous_vars, variable_labels, row_groups, title = "Table 1: Overall Episodic Summary by Group")
table1_pre  <- generate_table1(cab_pre,  group_var, categorical_vars, continuous_vars, variable_labels, row_groups, title = "Table 1: Pre-CAB Episodic Summary by Group")
table1_post <- generate_table1(cab_post, group_var, categorical_vars, continuous_vars, variable_labels, row_groups, title = "Table 1: Post-CAB Episodic Summary by Group")

#generate tables for patient data
table1_patient_all  <- generate_table1(patient_all,  group_var, categorical_vars_patient, continuous_vars_patient, variable_labels, row_groups, title = "Table 1: Overall Patient Summary by Group")
table1_patient_pre  <- generate_table1(patient_pre,  group_var, categorical_vars_patient, continuous_vars_patient, variable_labels, row_groups, title = "Table 1: Pre-CAB Patient Summary by Group")
table1_patient_post <- generate_table1(patient_post, group_var, categorical_vars_patient, continuous_vars_patient, variable_labels, row_groups, title = "Table 1: Post-CAB Patient Summary by Group")

#create a directory for table 1 outputs
dir.create("./Table1", showWarnings = FALSE)

#write episodic tables to csv
write.csv(table1_all$data, "./Table1/table1_overall_episodic.csv", row.names = FALSE)
write.csv(table1_pre$data, "./Table1/table1_pre_intervention_episodic.csv", row.names = FALSE)
write.csv(table1_post$data, "./Table1/table1_post_intervention_episodic.csv", row.names = FALSE)

#write patient tables to csv
write.csv(table1_patient_all$data, "./Table1/table1_overall_patient.csv", row.names = FALSE)
write.csv(table1_patient_pre$data, "./Table1/table1_pre_intervention_patient.csv",row.names = FALSE)
write.csv(table1_patient_post$data, "./Table1/table1_post_intervention_patient.csv", row.names = FALSE)
