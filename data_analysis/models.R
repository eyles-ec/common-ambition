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
                           origin_date = as.Date("2022-04-25"),
                           robust_if_overdispersed = TRUE,
                           fallback_family = "neg_binomial",
                           display = TRUE) {
  
  # Exclude data beyond December 31 2024 (Croydon tranche has further data)
  df <- df %>%
    dplyr::filter(.data[[time_var]] < 140)
  
  # Filter and factorise group variable; ensure period is 0/1 integer
  df_filtered <- df %>%
    dplyr::filter(.data[[time_var]] > time_cutoff,
                  .data[[group_var]] %in% groups_to_include) %>%
    dplyr::mutate(
      !!group_var := factor(.data[[group_var]], levels = groups_to_include),
      !!group_var := stats::relevel(.data[[group_var]], ref = reference_group),
      # Ensure period is numeric 0/1, coerce it
      !!period_var := as.integer(.data[[period_var]]),
      # Reconstruct a dates for each weekly time point to do month dummies
      # Later to do - fix cleaning part of pipeline to not require reconstruction
      week_date = origin_date + (.data[[time_var]] * 7),
      month_f = factor(format(week_date, "%m"),
                       levels = sprintf("%02d", 1:12),
                       labels = month.abb)
    )
  
  # CITS mean model + month fixed effects
  formula <- stats::as.formula(paste0(
    outcome_var, " ~ ",
    time_var, "*", group_var, "*", period_var, " + ",
    "month_f + ",
    "offset(log(", exposure_var, "))"
  ))
  
  #fit initial Poisson model
  poisson_model <- stats::glm(formula,
                              family = stats::poisson(link = "log"),
                              data = df_filtered)
  
  #Pearson check for overdispersion 
  dispersion <- sum(stats::residuals(poisson_model, type = "pearson")^2) / poisson_model$df.residual
  
  #Switch to NB (or quasi) if needed
  if (dispersion > 1.5 && robust_if_overdispersed) {
    warning(sprintf("Overdispersion detected (dispersion = %.2f). Switching to %s model.",
                    dispersion, fallback_family))
    
    model <- switch(fallback_family,
                    quasipoisson = stats::glm(formula,
                                              family = stats::quasipoisson(link = "log"),
                                              data = df_filtered),
                    neg_binomial = MASS::glm.nb(formula, data = df_filtered),
                    stop("Unsupported fallback_family. 'quasipoisson' or 'neg_binomial' only in function call."))
  } else {
    model <- poisson_model
  }
  
  #Add fitted values + pointwise CI for mean (link-scale SE)
  pred <- stats::predict(model, type = "link", se.fit = TRUE)
  df_filtered <- df_filtered %>%
    dplyr::mutate(
      yhat = exp(pred$fit),
      yhat_lower = exp(pred$fit - 1.96 * pred$se.fit),
      yhat_upper = exp(pred$fit + 1.96 * pred$se.fit)
    )
  
  if (display) {
    cat("\n--- Final Model Summary (", stats::family(model)$family, ") ---\n", sep = "")
    print(summary(model))
  }
  
  list(
    data = df_filtered,
    model = model,
    dispersion = dispersion,
    family_used = stats::family(model)$family
  )
}


#function to generate the counterfactual (if there was no CAB)
#needs the dataframe and model from CITS function above
#we generate this, then append it on the list output from the CITS function for the rest of the pipeline
#also possible to look at the counterfactual summary with the summary option on (by default, or off (summary =FALSE))
#because of autocorrelation we use Newey-West calculation for SE/CIs, nw_lag sets the lag for the calcualtion (in weeks)

generate_counterfactual <- function(df, model,
                                    group_name  = "Bristol ACHC",
                                    outcome_var = "hiv_test",
                                    time_var    = "time",
                                    group_var   = "group_bristol",
                                    period_var  = "period",
                                    level       = 0.95,
                                    nw_lag    = 4,
                                    prewhite  = TRUE,
                                    adjust    = TRUE,
                                    summary = TRUE) {
  
  #check arguments and columns since we are doing multiple outcomes/comparisons
  needed_cols <- c(outcome_var, time_var, group_var, period_var)
  missing_cols <- setdiff(needed_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("generate_counterfactual(): df is missing column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  
  #calculate constant for CIs 
  alpha <- 1 - level
  z <- stats::qnorm(1 - alpha/2)
  
  #evaluate offset() terms on newdata (works for glm / glm.nb)
  compute_offset <- function(model, newdata) {
    mf  <- stats::model.frame(stats::terms(model), data = newdata, na.action = stats::na.pass)
    off <- stats::model.offset(mf)
    if (is.null(off)) rep(0, nrow(newdata)) else as.numeric(off)
  }
  
  #Add row id to merge in counterfactuals later to the correct weeks
  df2 <- df |>
    dplyr::mutate(.row_id = dplyr::row_number())
  
  #Post-period rows for Bristol ACHC
  post <- df2 |>
    dplyr::filter(.data[[group_var]] == group_name,
                  .data[[period_var]] == 1)
  
  if (nrow(post) == 0) {
    out0 <- list(
      data = df2 |> dplyr::select(-.row_id),
      summary_table = dplyr::tibble(
        weeks_post = 0L,
        vcov_used = NA_character_
      ),
      cf_data = post
    )
    if (glimpse_summary) dplyr::glimpse(out0$summary_table)
    return(out0)
  }
  
  #Offsets (same exposure, same rows)
  offset_vec <- compute_offset(model, post)
  
  #Model matrices (aligned with coefficients), in order to calculate counterfactuals
  tt <- stats::delete.response(stats::terms(model))
  X  <- stats::model.matrix(tt, post)
  
  #Generate counterfactual design matrix, e.g. 
  #What would Bristol ACHC look like post-CAB if they followed the same post change
  #as the other group + common effects, without the Bristol ACHC specific trend
  # Keep post-period rows, but remove the treated group's post differential:
  # 1) Bristol ACHC group level change in post: group:period
  # 2) Bristol ACHC group slope change in post: time:group:period
  X_cf <- X
  
  group_token <- paste0(group_var, group_name)
  col_level_diff <- paste0(group_token, ":", period_var)
  col_slope_diff <- paste0(time_var, ":", group_token, ":", period_var)
  
  cols_to_zero <- c(col_level_diff, col_slope_diff)
  
  #check for missing columns in case of error with names
  missing_cols2 <- setdiff(cols_to_zero, colnames(X_cf))
  if (length(missing_cols2) > 0) {
    stop(
      "Controlled counterfactual could not find these columns in model matrix:\n  ",
      paste(missing_cols2, collapse = "\n  "),
      "\n\nAvailable columns include:\n  ",
      paste(colnames(X_cf), collapse = ", ")
    )
  }
  
  X_cf[, cols_to_zero] <- 0
  
  #generate variance-covariance matrix using Newey West technique,
  #which is more robust for autocorrelation present in data. (default for this is 4 weeks)
  V_full <- sandwich::NeweyWest(model, lag = nw_lag,
                                prewhite = prewhite,
                                adjust   = adjust)
  vcov_used <- paste0("NW_lag_", nw_lag)
  
  beta_full <- stats::coef(model)
  
  #ensure all the columns match between the cf and the vcov matrix
  common <- intersect(colnames(X), names(beta_full))
  common <- intersect(common, rownames(V_full))
  
  X    <- X[, common, drop = FALSE]
  X_cf <- X_cf[, common, drop = FALSE]
  beta <- beta_full[common]
  V2   <- V_full[common, common, drop = FALSE]
  
  #including the offset, calculate the mean and counterfacual mean
  #exponentiated due to model output
  eta_hat_post <- as.numeric(X %*% beta) + offset_vec
  yhat_post    <- exp(eta_hat_post)
  
  eta_cf_post <- as.numeric(X_cf %*% beta) + offset_vec
  cf_post     <- exp(eta_cf_post)
  
  y_obs      <- post[[outcome_var]]
  weeks_post <- nrow(post)
  
  #calculate totals for summary table 
  total_observed       <- sum(y_obs, na.rm = TRUE)
  total_counterfactual <- sum(cf_post, na.rm = TRUE)
  total_difference     <- sum(y_obs - cf_post, na.rm = TRUE)
  
  total_predicted  <- sum(yhat_post, na.rm = TRUE)
  total_pred_diff  <- total_predicted - total_counterfactual
  
  percent_increase   <- if (total_counterfactual > 0) 100 * total_difference / total_counterfactual else NA_real_
  pred_perc_increase <- if (total_counterfactual > 0) 100 * total_pred_diff / total_counterfactual else NA_real_
  
  annualised_difference     <- if (weeks_post > 0) total_difference * 52 / weeks_post else NA_real_
  average_weekly_difference <- if (weeks_post > 0) total_difference / weeks_post else NA_real_
  
  #use delta method to calculate slopes with matrix crossproducts
  g_cf   <- as.numeric(crossprod(X_cf, cf_post))    # d sum(cf) / d beta
  g_pred <- as.numeric(crossprod(X,    yhat_post))  # d sum(pred) / d beta
  
  se_total_cf   <- sqrt(drop(t(g_cf)   %*% V2 %*% g_cf))
  se_total_pred <- sqrt(drop(t(g_pred) %*% V2 %*% g_pred))
  
  #observed total is non-random, so the uncertainty comes from modelled counterfactual
  se_total_diff <- se_total_cf
  
  g_pred_diff <- g_pred - g_cf
  se_total_pred_diff <- sqrt(drop(t(g_pred_diff) %*% V2 %*% g_pred_diff))
  
  #the percent increase depends on the counterfactual 
  if (is.finite(total_counterfactual) && total_counterfactual > 0) {
    g_percent_inc <- 100 * (-total_observed / (total_counterfactual^2)) * g_cf
    se_percent_inc <- sqrt(drop(t(g_percent_inc) %*% V2 %*% g_percent_inc))
    
    g_pred_percent_inc <- 100 * ((1/total_counterfactual) * g_pred -
                                   (total_predicted/(total_counterfactual^2)) * g_cf)
    se_pred_percent_inc <- sqrt(drop(t(g_pred_percent_inc) %*% V2 %*% g_pred_percent_inc))
  } else {
    se_percent_inc <- NA_real_
    se_pred_percent_inc <- NA_real_
  }
  
  se_annualised_diff <- if (weeks_post > 0) se_total_diff * 52 / weeks_post else NA_real_
  se_avg_weekly_diff <- if (weeks_post > 0) se_total_diff / weeks_post else NA_real_
  
  #mini helper function to calculate cis for everything
  ci <- function(est, se, lower_bound = -Inf) {
    if (!is.finite(est) || !is.finite(se)) return(c(NA_real_, NA_real_))
    c(max(lower_bound, est - z * se), est + z * se)
  }
  
  #generate cis fof cf, differences, etc. 
  ci_total_cf        <- ci(total_counterfactual, se_total_cf, lower_bound = 0)
  ci_total_diff      <- ci(total_difference, se_total_diff)
  ci_total_pred      <- ci(total_predicted, se_total_pred, lower_bound = 0)
  ci_total_pred_diff <- ci(total_pred_diff, se_total_pred_diff)
  
  ci_percent_inc  <- ci(percent_increase, se_percent_inc)
  ci_pred_percent <- ci(pred_perc_increase, se_pred_percent_inc)
  ci_annualised   <- ci(annualised_difference, se_annualised_diff)
  ci_avg_weekly   <- ci(average_weekly_difference, se_avg_weekly_diff)
  
  #join per-week counterfactual back into the data using the row ids we made earlier
  cf_tbl <- dplyr::tibble(.row_id = post$.row_id, cf = cf_post)
  
  #add cf back into the original data, makes it easier to plot 
  out <- df2 |>
    dplyr::left_join(cf_tbl, by = ".row_id") |>
    dplyr::mutate(
      cf = dplyr::if_else(.data[[group_var]] == group_name & .data[[period_var]] == 1, cf, NA_real_),
      cf_diff = dplyr::if_else(.data[[group_var]] == group_name & .data[[period_var]] == 1,
                               .data[[outcome_var]] - cf, NA_real_)
    ) |>
    dplyr::select(-.row_id)
  
  #also return a tidy summary table for reporting 
  summary_df <- dplyr::tibble(
    total_observed = total_observed,
    total_counterfactual = total_counterfactual,
    total_counterfactual_lower = ci_total_cf[1],
    total_counterfactual_upper = ci_total_cf[2],
    total_difference = total_difference,
    total_difference_lower = ci_total_diff[1],
    total_difference_upper = ci_total_diff[2],
    total_predicted = total_predicted,
    total_predicted_lower = ci_total_pred[1],
    total_predicted_upper = ci_total_pred[2],
    total_pred_diff = total_pred_diff,
    total_pred_diff_lower = ci_total_pred_diff[1],
    total_pred_diff_upper = ci_total_pred_diff[2],
    percent_increase = percent_increase,
    percent_increase_lower = ci_percent_inc[1],
    percent_increase_upper = ci_percent_inc[2],
    pred_perc_increase = pred_perc_increase,
    pred_perc_increase_lower = ci_pred_percent[1],
    pred_perc_increase_upper = ci_pred_percent[2],
    weeks_post = weeks_post,
    annualised_difference = annualised_difference,
    annualised_difference_lower = ci_annualised[1],
    annualised_difference_upper = ci_annualised[2],
    average_weekly_difference = average_weekly_difference,
    average_weekly_difference_lower = ci_avg_weekly[1],
    average_weekly_difference_upper = ci_avg_weekly[2],
    vcov_used = vcov_used
  )
  
  #stick it all together in a list
  cf_data <- out |>
    dplyr::filter(.data[[group_var]] == group_name,
                  .data[[period_var]] == 1) |>
    dplyr::select(dplyr::all_of(c(time_var, group_var, period_var)), cf)
  
  result <- list(data = out, summary_table = summary_df, cf_data = cf_data)
  
  #print an easy summary of the cf table if you want 
  if (summary) dplyr::glimpse(result$summary_table)
  
  return(result)
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
    geom_point(alpha = 0.5)+
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
  ),
  current_prep = list(
    outcome_var = "current_prep",
    label = "current PrEP prescriptions"
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

# # TRUE to save plots as PNGs through your plotting functions (false if unneeded)
# save_plots <- TRUE

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
    # write a message with info for debugging 
    
    cf <- generate_counterfactual(
      df = cits$data,
      model = cits$model,
      group_name  = "Bristol ACHC",
      outcome_var = specification$outcome_var,
      time_var    = "time",
      group_var   = "group_bristol",
      period_var  = "period",
      summary = TRUE
    )
    
    cits_out <- list(
      model = cits$model,
      family_used = cits$family_used,
      dispersion = cits$dispersion,
      
      data = cf$data,                        # contains yhat + cf etc.
      summary_table = cf$summary_table,      # keep numeric internally (recommended)

      meta = list(
        outcome_name = outcome_name,
        outcome_var = specification$outcome_var,
        outcome_label = specification$label,
        comparison_name = comparison_name,
        groups = comparison$groups,
        reference_group = comparison$ref
      )
    )
    
    if (!outcome_name %in% names(results)) results[[outcome_name]] <- list()
    results[[outcome_name]][[comparison_name]] <- cits_out
    #create and save plots (CITS and slope comparison)
   
     #plots are also printed 
    p1 <- plot_cits(
      list(data = cits_out$data),
      outcome_var = specification$outcome_var,
      outcome_label = specification$label,
      group_var = "group_bristol",
      groups_to_include = comparison$groups,
      save_plot = FALSE,
      filename  = paste0("./Analysis/plots/newcits_",
                         comparison_name, "_", outcome_name, ".png")
    )
    print(p1)
  
    #tidy and save model specs and summary
    tidy_df <- broom::tidy(cits_out$model)
    write.csv(
      tidy_df,
      paste0("./Analysis/results/newcits_", comparison_name, "_",
             outcome_name, "_model.csv"),
      row.names = FALSE
    )

    write.csv(
      cits_out$summary_table,
      paste0("./Analysis/results/newcits_", comparison_name, "_",
             outcome_name, "_summary_table.csv"),
      row.names = FALSE
    )
  }
}


