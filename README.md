# Common Ambition Bristol
This repository contains code relevant to the quantitative evlaluation of Common Ambition Bristol (CAB) by NIHR ARC West. 

A project page with a plain English summary can be found [here](https://arc-w.nihr.ac.uk/research/projects/common-ambition-bristol-addressing-hiv-stigma-and-testing-in-partnership-with-african-caribbean-communities/)

## CAB Project Aims

The project aims to reverse HIV health inequalities experienced by people of African and Caribbean heritage living in Bristol and the surrounding area by:

- Increasing HIV knowledge
- Reducing HIV stigma
- Increasing uptake of wider sexual health services​
- Increasing uptake of HIV testing and PrEP (a pill that can stop you getting HIV)

## Specific quantitative objectives

## Quantitative Analysis

### Statistical Analysis

[@eyles-ec](https://github.com/eyles-ec) will lead the statistical analysis, which aims to study the impact of the Common Ambition Bristol (CAB) on HIV testing and use of sexual health services in the African Caribbean Heritage Community (ACHC) in Bristol. This will be done using anonymised patient level Electronic Health Records from Bristol and Croydon sexual health services. 

The following outcomes are to be analysed as weekly counts, using Poisson or negative binomial regression:

- HIV tests
- HIV diagnosis
- PrEP prescriptions
- Episodes of care to assess general use of sexual health services
- STI tests for any infection to assess STI testing incidence 

Bristol ACHC are compared to Bristol non-ACHC and also, separately, Croydon ACHC. Population offsets are calculated in order to make it possible to compare rates between the main and control series. 

### Health Economics Analysis

## Guide to the Repo

Data cleaning steps can be found in [/data_cleaning](https://github.com/eyles-ec/common-ambition/tree/main/data_cleaning). 

- [data_clean.R](https://github.com/eyles-ec/common-ambition/blob/main/data_cleaning/data_clean.R) shows the process from the raw data extract CSVs from Unity in Bristol and Croydon services. It includes functions for generating flag variables from codelists, appending multiple dataframes (including type coercion), and processing the data with new variables. 
- [data_collapse.R](https://github.com/eyles-ec/common-ambition/blob/main/data_cleaning/data_collapse.R) shows the collapse to weekly time series data of the cleaned data from the above. It also includes calculating population offset, and calculating rates for each outcome. It adds a grouping variable for analysis (Bristol ACHC, Bristol non ACHC, Croydon ACHC, Croydon non ACHC {not used, but preserved in the data}).

Data analysis steps can be found in[/data analysis](https://github.com/eyles-ec/common-ambition/tree/main/data_analysis).

- [loess_graphs.R](https://github.com/eyles-ec/common-ambition/blob/main/data_analysis/loess_graphs.R) shows the process for plotting the weekly data by the grouping varaibles and a chosen outcome variable, including a save function with overwrite protection. It also includes a LOESS line of trend. 
- [models.R](https://github.com/eyles-ec/common-ambition/blob/main/data_analysis/models.R) includes a function to fit CITS models, with several user-defined parameters. It also includes a function to generate counterfactuals based on a CITS model generated there, and a summary table which is exportable, including a counterfactual-observed comparison, with confidence intervals estimated with the delta method. It also has a function to make a classic CITS plot comparing two groups, including a line based on the generated counterfactual for Bristol ACHC. This is saveable, with overwrite protection. The functions are run in a for loop so outcomes can be added or removed easily from the piepline

Table creation steps can be found in [/tables](https://github.com/eyles-ec/common-ambition/tree/main/tables).

- [table_1.R](https://github.com/eyles-ec/common-ambition/blob/main/tables/table_1.R) allows for the creation of descriptive summary tables on both continuous and categorical variables, for both episodic level and patient level data, to be exported to CSV. There is also a redaction helper function to replace counts of <6 with '<6.' It also creates a gt table, which can be exported to LaTeX if desired. 
