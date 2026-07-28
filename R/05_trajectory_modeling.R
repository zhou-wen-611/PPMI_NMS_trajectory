###############################################################################
# 05_trajectory_modeling.R
# Latent class trajectory modeling (Supplementary Figures)
# Author: Wen Zhou
# Date: 2026-07-24
# R version: 4.2.2
###############################################################################

library(lcmm)
library(dplyr)

source("R/functions.R")

# Load longitudinal data
longitudinal_data <- readRDS("data/longitudinal_cleaned.rds")

#------------------------------------------------------------------------------
# 1. DEFINE OUTCOMES FOR TRAJECTORY ANALYSIS
#------------------------------------------------------------------------------

trajectory_outcomes <- c("ESS", "REM", "GDS", "STAI_total", "STAI_state", 
                         "STAI_trait", "SCOPA_AUT", "MoCA")

#------------------------------------------------------------------------------
# 2. FIT TRAJECTORY MODELS
#------------------------------------------------------------------------------

fit_trajectory_model <- function(data, outcome, max_classes = 5) {
  # Remove missing values for the outcome
  data <- data[!is.na(data[[outcome]]), ]
  
  if (length(unique(data$PATNO)) < 50) {
    warning("Not enough participants for trajectory modeling")
    return(NULL)
  }
  
  # Fit models with different number of classes (2 to max_classes)
  models <- list()
  for (ng in 2:max_classes) {
    tryCatch({
      model <- hlme(
        fixed = as.formula(paste(outcome, "~ visit_time + I(visit_time^2)")),
        random = ~ visit_time,
        mixture = ~ visit_time + I(visit_time^2),
        ng = ng,
        data = data,
        subject = "PATNO",
        maxiter = 100
      )
      models[[as.character(ng)]] <- model
    }, error = function(e) {
      warning("Failed to fit ", ng, "-class model for ", outcome, ": ", e$message)
    })
  }
  
  if (length(models) == 0) {
    return(NULL)
  }
  
  # Select best model based on BIC and other criteria
  bics <- sapply(models, function(m) {
    if (!is.null(m)) summary(m)$BIC else NA
  })
  
  # Choose model with lowest BIC, or if BIC is similar, prefer parsimonious
  selected <- models[[which.min(bics)]]
  
  # Calculate APP and entropy for the selected model
  if (!is.null(selected)) {
    post_prob <- selected$pprob[, -1]  # Remove first column (subject ID)
    class_assign <- apply(post_prob, 1, which.max)
    app_values <- calculate_app(post_prob, class_assign)
    entropy <- calculate_entropy(post_prob)
    
    cat(outcome, ":", "Selected", selected$ng, "classes,",
        "Entropy =", round(entropy, 3), "\n")
    cat("  APP:", paste(round(app_values, 3), collapse = ", "), "\n")
  }
  
  return(selected)
}

# Fit models for all outcomes
trajectory_models <- list()

for (outcome in trajectory_outcomes) {
  cat("\nFitting trajectory models for", outcome, "...\n")
  trajectory_models[[outcome]] <- fit_trajectory_model(
    longitudinal_data, outcome, max_classes = 5
  )
}

#------------------------------------------------------------------------------
# 3. SAVE RESULTS
#------------------------------------------------------------------------------

# Save trajectory models
saveRDS(trajectory_models, "output/trajectory_models.rds")

# Extract and save class assignments
class_assignments <- list()

for (outcome in names(trajectory_models)) {
  model <- trajectory_models[[outcome]]
  if (!is.null(model)) {
    post_prob <- model$pprob
    class_assign <- apply(post_prob[, -1], 1, which.max)
    class_assignments[[outcome]] <- data.frame(
      PATNO = post_prob[, 1],
      Class = class_assign
    )
  }
}

saveRDS(class_assignments, "output/class_assignments.rds")

cat("\nTrajectory modeling complete.\n")

# Save session info
writeLines(capture.output(sessionInfo()), "output/session_info_05.txt")
