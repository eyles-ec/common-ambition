library(tidyverse)

#helper function for summary tables
#set file names and then comparisons via the file names into a tibble
#includes the same for censored analysis
parse_meta <- function(file) {
  name <- basename(file)
  
  tibble(
    file = file,
    type = case_when(
      str_detect(name, "^CScits|^CSfisher") ~ "52_week", #detect if it's censored
      TRUE ~ "full"
    ),
    analysis = case_when(
      str_detect(name, "model") ~ "model",
      str_detect(name, "summary_table") ~ "summary",
      str_detect(name, "fisher") ~ "fisher"
    ),
    comparison = str_extract(name, "bristol_within|bristol_vs_croydon"),
    outcome = str_extract(name, "(hiv_test|new_hiv|sti_test_count_hiv|sti_test_count_no_hiv|weekly_episode_count|current_prep)")
  )
}

#return a formatted string for the estimate and confidence interval
fmt_ci <- function(est, low, high, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"), est, low, high)
}

#function for pulling out a model coefficient from the 'term' column, returning the first match 
get_row <- function(df, pattern) {
  df %>% filter(str_detect(term, pattern)) %>% slice(1)
}

#set results directory and list files that are available there
results_dir <- "./Analysis/results_corrected"
files <- list.files(results_dir, pattern = "\\.csv$", full.names = TRUE)

#run the parseing function over all files in the list using map_dfr from purrr
#repeats for each file then binds the results together

all_files <- map_dfr(files, parse_meta)

#pull out fisher results, including rounded p
fisher_tbl <- all_files %>%
  filter(analysis == "fisher") %>%
  mutate(data = map(file, read_csv)) %>%
  mutate(
    p_value = map_dbl(data, ~ .x %>%
                        filter(group == "p_value") %>%
                        select(where(is.numeric)) %>%
                        pull(1))
  ) %>%
  mutate(
    p_value_fmt = sprintf("%.3f", p_value)
  ) %>%
  select(type, comparison, outcome, p_value, p_value_fmt)

#pull out summary results for reporting from counterfactual/model summary table
summary_tbl <- all_files %>%
  filter(analysis == "summary") %>%
  mutate(data = map(file, read_csv)) %>%
  unnest(data) %>%
  mutate(
    difference_ci = fmt_ci(total_difference, #cf difference across the whole post period
                           total_difference_lower,
                           total_difference_upper),
    
    weekly_diff_ci = fmt_ci(average_weekly_difference, #avg weekly diff calculated elsewhere
                            average_weekly_difference_lower,
                            average_weekly_difference_upper)
  ) %>%
  transmute( #keep only what's listed in the table
    type,
    comparison,
    outcome,
    
    #total difference across all periods
    total_difference,
    total_difference_lower,
    total_difference_upper,
    difference_ci,
    
    #weekly difference across all period
    weekly_diff = average_weekly_difference,
    weekly_diff_lower = average_weekly_difference_lower,
    weekly_diff_upper = average_weekly_difference_upper,
    weekly_diff_ci
  )

#clean table of model coefficients

coef_tbl <- all_files %>%
  filter(analysis == "model") %>%
  mutate(data = map(file, read_csv)) %>%
  mutate( #pull out relevant coefficients for table
    baseline = map(data, ~ get_row(.x, "^group_bristol")),
    baseline_trend = map(data, ~ get_row(.x, "^time:group_bristol")),
    step = map(data, ~ get_row(.x, "group_bristol.*:period$")),
    post_trend = map(data, ~ get_row(.x, "time:group_bristol.*:period"))
    
  ) %>%
  mutate(
    #extract log-coefficients
    baseline_est = map_dbl(baseline, ~ .x$estimate),
    baseline_low = map_dbl(baseline, ~ .x$estimate - 1.96 * .x$std.error),
    baseline_high = map_dbl(baseline, ~ .x$estimate + 1.96 * .x$std.error),
    baseline_p = map_dbl(baseline, ~.x$p.value),
    
    trend_est = map_dbl(baseline_trend, ~ .x$estimate),
    trend_low = map_dbl(baseline_trend, ~ .x$estimate - 1.96 * .x$std.error),
    trend_high = map_dbl(baseline_trend, ~ .x$estimate + 1.96 * .x$std.error),
    trend_p = map_dbl(baseline_trend, ~.x$p.value),
    
    step_est = map_dbl(step, ~ .x$estimate),
    step_low = map_dbl(step, ~ .x$estimate - 1.96 * .x$std.error),
    step_high = map_dbl(step, ~ .x$estimate + 1.96 * .x$std.error),
    step_p = map_dbl(step, ~.x$p.value),
    
    post_est = map_dbl(post_trend, ~ .x$estimate),
    post_low = map_dbl(post_trend, ~ .x$estimate - 1.96 * .x$std.error),
    post_high = map_dbl(post_trend, ~ .x$estimate + 1.96 * .x$std.error),
    post_p = map_dbl(post_trend, ~.x$p.value)
  ) %>%
  mutate(
    #exponentiate into rate ratios 
    baseline_rr = exp(baseline_est),
    baseline_rr_low = exp(baseline_low),
    baseline_rr_high = exp(baseline_high),
    
    trend_rr = exp(trend_est),
    trend_rr_low = exp(trend_low),
    trend_rr_high = exp(trend_high),
    
    step_rr = exp(step_est),
    step_rr_low = exp(step_low),
    step_rr_high = exp(step_high),
    
    post_rr = exp(post_est),
    post_rr_low = exp(post_low),
    post_rr_high = exp(post_high)
  ) %>%
  mutate(
    baseline_rr_ci = fmt_ci(baseline_rr, baseline_rr_low, baseline_rr_high),
    trend_rr_ci = fmt_ci(trend_rr, trend_rr_low, trend_rr_high),
    step_rr_ci = fmt_ci(step_rr, step_rr_low, step_rr_high),
    post_rr_ci = fmt_ci(post_rr, post_rr_low, post_rr_high),
    baseline_p_fmt = sprintf("%.3f", baseline_p),
    trend_p_fmt    = sprintf("%.3f", trend_p),
    step_p_fmt     = sprintf("%.3f", step_p),
    post_p_fmt     = sprintf("%.3f", post_p)
    )%>%
  select( #keep only relevant columns
    type, comparison, outcome,
    
    baseline_rr, baseline_rr_ci, baseline_p, baseline_p_fmt,
    trend_rr, trend_rr_ci, trend_p, trend_p_fmt,
    step_rr, step_rr_ci, step_p, step_p_fmt,
    post_rr, post_rr_ci, post_p, post_p_fmt
  )


#write all summary tables to CSV within results directory
write.csv(fisher_tbl,
          file.path(results_dir, "table_summaries/table_fisher_clean.csv"),
          row.names = FALSE)

write.csv(summary_tbl,
          file.path(results_dir, "table_summaries/table_cf_summary_clean.csv"),
          row.names = FALSE)

write.csv(coef_tbl,
          file.path(results_dir, "table_summaries/table_coefficients_clean.csv"),
          row.names = FALSE)

