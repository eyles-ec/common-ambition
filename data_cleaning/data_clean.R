library(dplyr)
library(tidyverse)
library(data.table)

#set working directory

setwd("Own path")
      
#load dataset

unity_cab1 <- read.csv("./subdirectory/unity.csv")

#load shappt reference and hiv codes

shappt <- read.csv("./subdirectory/codelist/shappt_codelist.csv")
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

#save dataset 

write.csv(unity_hiv.csv)


  