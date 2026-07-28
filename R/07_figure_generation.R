###############################################################################
# 07_figure_generation.R
# Figure generation (Figures 1-8, Supplementary Figures)
# Author: Wen Zhou
# Date: 2026-07-24
# R version: 4.2.2
###############################################################################

library(ggplot2)
library(forestplot)
library(patchwork)

source("R/functions.R")

# Load data
forest_data <- readRDS("output/forest_data.rds")

#------------------------------------------------------------------------------
# 1. CREATE FOREST PLOTS (Figures 2-8)
#------------------------------------------------------------------------------

# Pre-defined color palette
colors_fresh <- c("#4A90D9", "#2ECC71", "#5BC0DE", "#F5B7B1", "#48C9B0", 
                  "#F8C471", "#85C1E9", "#82E0AA", "#A3E4D7", "#D7BDE2", 
                  "#F9E79F", "#AED6F1", "#A9DFBF", "#FADBD8")

create_forest_plot <- function(data, outcome_name) {
  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }
  
  data <- data[order(data$Beta), ]
  data$Subgroup <- factor(data$Subgroup, levels = data$Subgroup)
  
  p <- ggplot(data, aes(x = Beta, y = Subgroup)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), 
                   height = 0.2, color = "#4A90D9", linewidth = 0.8) +
    geom_point(aes(size = N), color = "#2C5F8A", shape = 16) +
    scale_size_continuous(range = c(2, 5), guide = guide_legend(title = "N")) +
    labs(x = "Coefficient (95% CI)", y = NULL, 
         title = outcome_name) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  return(p)
}

outcomes <- c("ESS", "REM", "GDS", "STAI_total", "STAI_state", 
              "STAI_trait", "SCOPA_AUT")

figures <- list()

for (outcome in outcomes) {
  data <- forest_data[[outcome]]
  if (!is.null(data) && nrow(data) > 0) {
    p <- create_forest_plot(data, outcome)
    figures[[outcome]] <- p
    # Save individual figure
    # ggsave(paste0("output/Figure_", match(outcome, outcomes) + 1, "_", outcome, ".pdf"), 
    #        p, width = 8, height = 6)
  }
}

cat("Forest plots generated (Figures 2-8).\n")

#------------------------------------------------------------------------------
# 2. CREATE TRAJECTORY PLOTS (Supplementary Figures)
#------------------------------------------------------------------------------

trajectory_models <- readRDS("output/trajectory_models.rds")
longitudinal_data <- readRDS("data/longitudinal_cleaned.rds")

create_trajectory_plot <- function(model, data, outcome_name) {
  if (is.null(model)) return(NULL)
  
  # Extract predictions
  pred <- tryCatch({
    predict_lcmm(model, newdata = data, var.time = "visit_time")
  }, error = function(e) {
    warning("Prediction failed for ", outcome_name, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(pred)) {
    # Fallback: simple plot using observed data
    p <- ggplot(data, aes(x = visit_time, y = .data[[outcome_name]])) +
      geom_smooth(aes(group = PATNO), se = FALSE, alpha = 0.1, linewidth = 0.3) +
      labs(x = "Follow-up (years)", y = outcome_name,
           title = paste(outcome_name, "Trajectories")) +
      theme_minimal()
    return(p)
  }
  
  pred_df <- as.data.frame(pred)
  p <- ggplot(pred_df, aes(x = time, y = pred, color = factor(class))) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(aes(ymin = pred - 1.96 * se, ymax = pred + 1.96 * se, 
                    fill = factor(class)), alpha = 0.2, color = NA) +
    labs(x = "Follow-up (years)", y = outcome_name,
         title = paste(outcome_name, "Trajectory Classes"),
         color = "Class", fill = "Class") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

for (outcome in names(trajectory_models)) {
  p <- create_trajectory_plot(trajectory_models[[outcome]], 
                              longitudinal_data, outcome)
  if (!is.null(p)) {
    # ggsave(paste0("output/Supplementary_", outcome, "_trajectory.pdf"), 
    #        p, width = 8, height = 6)
  }
}

cat("Trajectory plots generated (Supplementary Figures).\n")

#------------------------------------------------------------------------------
# 3. FLOWCHART (Figure 1) - Text representation
#------------------------------------------------------------------------------

cat("\n=== Figure 1: Flowchart (Text Representation) ===\n")
cat("
+--------------------------------------+
| PPMI database initial screening       |
| (N = X,XXX)                          |
+--------------------------------------+
                |
                v
+--------------------------------------+
| Excluded (n = 226):                  |
| - Missing motor subtype (n = 91)     |
| - No follow-up data (n = 135)        |
+--------------------------------------+
                |
                v
+--------------------------------------+
| Final analytic sample                 |
| (N = 1,384)                          |
+--------------------------------------+
                |
    +-----------+-----------+
    |                       |
    v                       v
+------------------+  +------------------+
| Cross-sectional  |  | Longitudinal     |
| analysis         |  | analysis         |
| (N = 1,384)      |  | (N = 1,249)      |
+------------------+  +------------------+
")

# Save session info
writeLines(capture.output(sessionInfo()), "output/session_info_07.txt")

cat("\nAll figure generation complete.\n")
