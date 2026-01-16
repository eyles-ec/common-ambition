library(stats)
library(tidyverse)
library(MASS)
library(broom)

#function to fit CITS models, with a 6 week 'burn in period' that can be changed as we have data back to -59.5
#the 0 in the time variable in this case is the CAB intervention date in April 2022, represented in time_var
#time_cutoff thus is giving us the burn in period
#outcome_var can be changed to which outcome you like, default is hiv test count
#group_var is the grouping of location and ACHC
#exposure_var is the population of the specific group selected, calculated in the collapse code
#groups_to_include are which comparison you're making with group_var, default is Bristol ACHC v Bristol non ACHC
#reference_group allows you to set the ref group for the model pf above
#robust_if_overdispersed allows the function to try negative binomial by default (but there's also an option for quasipoisson)
#saves the model output, dispersion, family used, and predicted values for use in the counterfactual in a list
#also prints the final selected model if display = TRUE (by default)

#due to the additional data, i have added a time_since_intervention term to the models
#this can tell us effectively the weekly % increase relative to the comparison group
#make sure to exponentiate the coefficient to interpret it 

#also restricts data to the end of Bristol data in week 140 from intervention


fit_cits_model <- function(df, outcome_var = "hiv_test", 
                           group_var = "group_bristol", 
                           time_var = "time", 
                           period_var = "period", 
                           exposure_var = "population", 
                           groups_to_include = c("Bristol ACHC", "Bristol non ACHC"),
                           reference_group = "Bristol non ACHC",
                           time_cutoff = -53,
                           robust_if_overdispersed = TRUE,
                           fallback_family = "neg_binomial",
                           display = TRUE) {
  
  if (fallback_family == "neg_binomial" && !requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for negative binomial modeling. Please install it.")
  }
  
  # Exclude all data from and after the week of 31 Dec 2024 (time index 140)
  #as no Bristol data past that point 
  df <- df %>%
    filter(.data[[time_var]] < 140)
  
  # Filter and factorise group variable
  df_filtered <- df %>%
    filter(.data[[time_var]] > time_cutoff,
           .data[[group_var]] %in% groups_to_include) %>%
    mutate(
      !!group_var := factor(.data[[group_var]], levels = groups_to_include),
      !!group_var := relevel(.data[[group_var]], ref = reference_group),
      time_since_intervention = ifelse(.data[[period_var]] == 1, .data[[time_var]], 0)
    )
  
  # Build formula with time_since_intervention
  formula <- as.formula(paste0(outcome_var, " ~ ",
                               time_var, "*", group_var, "*", period_var, " + ",
                               "time_since_intervention *", group_var, " + ",
                               "offset(log(", exposure_var, "))"))
  
  # Fit initial Poisson model
  poisson_model <- glm(formula,
                       family = poisson(link = "log"),
                       data = df_filtered)
  
  # Check for overdispersion
  dispersion <- sum(residuals(poisson_model, type = "pearson")^2) / poisson_model$df.residual
  
  # Switch to robust family if needed
  if (dispersion > 1.5 && robust_if_overdispersed) {
    warning(paste("Overdispersion detected (dispersion =", round(dispersion, 2), 
                  "). Switching to", fallback_family, "model."))
    
    model <- switch(fallback_family,
                    quasipoisson = glm(formula,
                                       family = quasipoisson(link = "log"),
                                       data = df_filtered),
                    neg_binomial = MASS::glm.nb(formula,
                                                data = df_filtered),
                    stop("Unsupported fallback_family. Choose 'quasipoisson' or 'neg_binomial'."))
  } else {
    model <- poisson_model
  }
  
  # Add predicted values and confidence intervals
  pred <- predict(model, type = "link", se.fit = TRUE)
  df_filtered <- df_filtered %>%
    mutate(
      yhat = exp(pred$fit),
      yhat_lower = exp(pred$fit - 1.96 * pred$se.fit),
      yhat_upper = exp(pred$fit + 1.96 * pred$se.fit)
    )
  
  # Print model summary
  if (display) {
    cat("\n--- Final Model Summary (", family(model)$family, ") ---\n", sep = "")
    print(summary(model))
  }
  
  return(list(
    data = df_filtered,
    model = model,
    dispersion = dispersion,
    family_used = family(model)$family
  ))
}



#function to generate the counterfactual (if there was no CAB)
#needs the dataframe and model from CITS function above
#we generate this, then append it on the list output from the CITS function for the rest of the pipeline

generate_counterfactual <- function(df, model, 
                                    group_name = "Bristol ACHC", 
                                    outcome_var = "hiv_test", 
                                    time_var = "time",
                                    group_var = "group_bristol",
                                    period_var = "period",
                                    offset_var = "population") {
  
  coefs <- coef(model)
  
  # Safely extract coefficients with fallback to 0 if missing
  get_coef <- function(name) ifelse(name %in% names(coefs), coefs[[name]], 0)
  
  # Extract relevant coefficients for counterfactual
  intercept      <- get_coef("(Intercept)")
  group_effect   <- get_coef(paste0(group_var, group_name))
  time_effect    <- get_coef(time_var)
  time_group     <- get_coef(paste0(time_var, ":", group_var, group_name))
  
  # Compute slopes
  slope_reference_group <- time_effect
  slope_target_group <- time_effect + time_group
  slope_difference <- slope_target_group - slope_reference_group
  
  # Apply counterfactual only to post-intervention period for target group
  df <- df %>%
    mutate(cf = if_else(.data[[group_var]] == group_name & .data[[period_var]] == 1,
                        exp(intercept + group_effect +
                              .data[[time_var]] * slope_reference_group) * .data[[offset_var]],
                        NA_real_),
           cf_diff = if_else(.data[[group_var]] == group_name & .data[[period_var]] == 1,
                             .data[[outcome_var]] - cf,
                             NA_real_))
  
  # Create summary table for outputs
  summary_df <- df %>%
    filter(.data[[group_var]] == group_name, .data[[period_var]] == 1) %>%
    summarise(
      total_observed = sum(.data[[outcome_var]], na.rm = TRUE),
      total_counterfactual = sum(cf, na.rm = TRUE),
      total_difference = sum(cf_diff, na.rm = TRUE),
      total_predicted = sum(yhat, na.rm = TRUE),
      total_pred_diff = total_predicted - total_counterfactual,
      pred_perc_increase = 100 * total_pred_diff / total_counterfactual,
      weeks_post = n(),
      annualised_difference = total_difference * 52 / weeks_post,
      average_weekly_difference = total_difference / weeks_post,
      percent_increase = 100 * total_difference / total_counterfactual,
      slope_target_group = slope_target_group,
      slope_reference_group = slope_reference_group,
      slope_difference = slope_difference
    )
  
  return(list(
    data = df,
    summary_table = summary_df
  ))
}

#create a standardised CITS plot, including the counterfactual line for Bristol ACHC
#can save the plot if desired with own filename, overwrite protection. you also can scale it
#height and width are in inches by default 

plot_cits <- function(cits_model,
                      outcome_var = "hiv_test",        # NEW: outcome column
                      outcome_label = "HIV tests",     # NEW: label for plot
                      group_var = "group_bristol",
                      groups_to_include = c("Bristol ACHC", "Bristol non ACHC"),
                      save_plot = FALSE,
                      filename = "./subdirectory/plots/cits_plot.png",
                      overwrite = FALSE,
                      width = 10,
                      height = 6,
                      dpi = 300) {
  
  df_plot <- cits_model$data
  
  # Dynamically extract the group variable
  df_plot <- df_plot %>%
    mutate(group = .data[[group_var]]) %>%
    filter(group %in% groups_to_include)
  
  # Counterfactual line only for Bristol ACHC
  cf_data <- df_plot %>%
    filter(group == "Bristol ACHC", !is.na(cf))
  
  # Build plot
  p <- ggplot(df_plot, aes(x = time, y = .data[[outcome_var]], color = group)) +
    geom_point(alpha = 0.5) +
    geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper, fill = group),
                alpha = 0.2, color = NA) +
    geom_line(aes(y = yhat), size = 1.2) +
    geom_line(data = cf_data,
              aes(y = cf), linetype = "dotted", color = "grey40", size = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = paste("Observed vs Predicted", outcome_label, "Over Time"),
      subtitle = "Dashed line = CAB interventions (Week of April 25, 2022)\nDotted line = counterfactual prediction for Bristol ACHC",
      x = "Weeks Since Intervention",
      y = outcome_label,
      color = "Group",
      fill = "Group"
    ) +
    facet_wrap(~ group, scales = "free_y") +
    theme_minimal()
  
  # Save plot if requested
  if (save_plot) {
    if (file.exists(filename) && !overwrite) {
      warning(paste("File", filename, "already exists. Set overwrite = TRUE to replace it."))
    } else {
      ggsave(filename, plot = p, width = width, height = height, dpi = dpi)
    }
  }
  
  return(p)
}


#creates a plot comparing the slopes of the two groups, up to x_end of 140 weeks (end of Bristol data)
#this is customisable to however many weeks you desire, include calculation of slope difference
#uses hardcoded names of the coefficients as they don't change between models, but check this each time


plot_slope_comparison <- function(cits_model,
                                  outcome_label = "HIV tests",     # NEW: single argument used everywhere
                                  group_labels = c("Bristol ACHC", "Bristol non ACHC"),
                                  x_end = 140,
                                  save_plot = FALSE,
                                  filename = "./subdirectory/plots/slope_comparison.png",
                                  overwrite = FALSE,
                                  width = 8,
                                  height = 6,
                                  dpi = 300) {
  # Extract coefficients (hardcoded as dynamic way did not work)
  coefs <- coef(cits_model$model)
  slope_control <- coefs["time_since_intervention"]
  slope_achc <- slope_control + coefs["group_bristolBristol ACHC:time_since_intervention"]
  
  # Calculate endpoints and difference for the slopes
  y_achc <- slope_achc * x_end
  y_control <- slope_control * x_end
  slope_gap <- y_achc - y_control
  y_midpoint <- (y_achc + y_control) / 2
  
  # Create slope data with dynamic labels
  slope_df <- tibble(
    group = group_labels,
    x = 0,
    xend = x_end,
    y = 0,
    yend = c(y_achc, y_control)
  )
  
  # Build plot using gg plot, with information describing how to interpret it.
  p <- ggplot() +
    geom_segment(data = slope_df, aes(x = x, xend = xend, y = y, yend = yend, color = group), size = 1.5) +
    annotate("segment", x = x_end, xend = x_end, y = y_control, yend = y_achc,
             linetype = "dotted", color = "black", size = 1) +
    annotate("text", x = x_end - 0.3, y = y_midpoint,
             label = paste0("Difference: ", round(slope_gap, 2)),
             hjust = 1, size = 4, color = "black") +
    labs(
      title = "Post-Intervention Slope Comparison",
      subtitle = paste0(
        "Each line shows the estimated increase in ", outcome_label,
        " over time since the intervention.\n",
        "Y-axis reflects total additional ", outcome_label,
        ", not the weekly rate."
      ),
      x = "Weeks Since Intervention",
      y = paste("Estimated Additional", outcome_label),
      color = "Group"
    ) +
    theme_minimal()
  
  # Save plot if requested
  if (save_plot) {
    if (file.exists(filename) && !overwrite) {
      warning(paste("File", filename, "already exists. Set overwrite = TRUE to replace it."))
    } else {
      ggsave(filename, plot = p, width = width, height = height, dpi = dpi)
    }
  }
  
  return(p)
}

#set working directory
setwd(YOURWD)

#load analysis dataset
cab<- read.csv("./Analysis/weekly_combined.csv")

#define outcomes

outcomes <- list(
  hiv_test = list(
    outcome_var = "hiv_test",
    label = "HIV tests"
  ),
  new_hiv = list(
    outcome_var = "new_hiv",
    label = "New HIV diagnoses"
  ),
  sti_test_count = list(
    outcome_var = "sti_test_count",
    label = "STI tests"
  ),
  weekly_episode_count = list(
    outcome_var = "weekly_episode_count",
    label = "Weekly episodes"
  )
)

#comparisons 
comparisons <- list(
  bristol_within = list(
    groups = c("Bristol ACHC", "Bristol non ACHC"),
    ref = "Bristol non ACHC"
  ),
  bristol_vs_croydon = list(
    groups = c("Bristol ACHC", "Croydon ACHC"),
    ref = "Croydon ACHC"
  )
)

# output folders creation
dir.create("./Analysis/results", recursive = TRUE, showWarnings = FALSE)
dir.create("./Analysis/plots",   recursive = TRUE, showWarnings = FALSE)

# TRUE to save plots as PNGs through your plotting functions (false if unneeded)
save_plots <- TRUE

# Store outputs here
results <- list()

#run it in a loop for easy changes to outcomes etc
for (outcome_name in names(outcomes)) {
  specification <- outcomes[[outcome_name]]
  
  for (comparison_name in names(comparisons)) {
    comparison <- comparisons[[comparison_name]]
    
    message("Processing: ", outcome_name, " — ", comparison_name)
    
    # Fit the CITs model
    cits <- fit_cits_model(
      df = cab,
      outcome_var = specification$outcome_var,
      groups_to_include = comparison$groups,
      reference_group = comparison$ref
    )
    
    # Generate counterfactuals and attach back to object
    cf <- generate_counterfactual(df = cits$data, model = cits$model)
    cits$data <- cf$data
    cits$summary_table <- cf$summary_table
    
    # Store in nested list: results[[outcome]][[comparison]]
    if (!outcome_name %in% names(results)) results[[outcome_name]] <- list()
    results[[outcome_name]][[comparison_name]] <- cits
    
    #create and save plots (CITS and slope comparison)
    #plots are also printed 
    p1 <- plot_cits(
      cits,
      outcome_var = specification$outcome_var,
      outcome_label = specification$label,
      group_var = "group_bristol",
      groups_to_include = comparison$groups,
      save_plot = save_plots,
      filename  = paste0("./Analysis/plots/newcits_",
                         comparison_name, "_", outcome_name, ".png")
    )
    print(p1)
    
    p2 <- plot_slope_comparison(
      cits,
      outcome_label = specification$label,
      group_labels = comparison$groups,
      save_plot = save_plots,
      filename  = paste0("./Analysis/plots/newslope_",
                         comparison_name, "_", outcome_name, ".png")
    )
    print(p2)
    
    # tidy and save model specs 
    tidy_df <- broom::tidy(cits$model)
    write.csv(
      tidy_df,
      paste0("./Analysis/results/newcits_", comparison_name, "_",
             outcome_name, "_model.csv"),
      row.names = FALSE
    )
    
    write.csv(
      cits$summary_table,
      paste0("./Analysis/results/newcits_", comparison_name, "_",
             outcome_name, "_summary.csv"),
      row.names = FALSE
    )
  }
}


