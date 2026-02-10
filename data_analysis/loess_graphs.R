library(dplyr)
library(ggplot2)
library(rlang)
library(lubridate)

#plot function for all on one plot. Generates a plot either with everything by default on one plot
#setting facet to = TRUE will give you the plots for each group on separate facets 
#setting y_label to a string will label the y axis
#setting plot_title to a string will give you a plot title
#save_plot will save it with a default filename if TRUE
#filename can be set in the options
#overwrite option lets you overwrite files, otherwise you get a warning
#width and height are by default measured in inches
#dpi is dots per square inch 

plot_outcomes <- function(df, outcome_var, 
                          group_var = "group_bristol", 
                          time_var = "time", 
                          facet = FALSE, 
                          y_label = NULL, 
                          plot_title = NULL, 
                          save_plot = FALSE, 
                          filename = "./plots/plot_output.png",
                          overwrite = FALSE,
                          width = 8, height = 5, dpi = 300) {
  
  # Use the outcome variable name as the default if y_label isn't provided
  if (is.null(y_label)) {
    y_label <- outcome_var
  }
  
  # Build base plot
  p <- ggplot(df, aes_string(x = time_var, y = outcome_var, color = group_var)) +
    geom_line() +
    geom_smooth(method = "loess", se = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = paste("Trend of", outcome_var, if (facet) "by Group (Faceted)" else "by Group"),
      x = "Time",
      y = y_label,
      color = "Group"
    ) +
    theme_minimal()
  
  # Add title if one is specified
  if (!is.null(plot_title)) {
    p <- p + ggtitle(plot_title)
  }
  
  # Facet the plot if requested 
  if (facet) {
    p <- p + facet_wrap(as.formula(paste("~", group_var)))
  }
  
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

#collect paths - save your own paths to a file called paths.R that is ignored by git (.gitignore) 

source("../paths.R")

#set working directory 

setwd(wd)

#load analysis dataset
cab<- read.csv("./Analysis/weekly_combined.csv")


#call plot function
plot_outcomes(cab[cab$location == "Bristol",], "new_hiv_rate", y_label  = "Weekly episode rate per 1000", plot_title = "Weekly episode rate by group")
