###############################################################################
# 03_subgroup_analysis.R
# Subgroup analyses (Figures 2-8)
# Author: Wen Zhou
# Date: 2026-07-24
# R version: 4.2.2
###############################################################################

library(dplyr)
library(ggplot2)
library(tidyr)

source("R/functions.R")

# Load cleaned data
baseline_data <- readRDS("data/baseline_cleaned.rds")

# Define covariates (same as in cross-sectional analysis)
covariates <- c("age", "sex", "age_at_onset", "race", "BMI", 
                "education", "SBR", "HY_stage", "UPDRS3")

#------------------------------------------------------------------------------
# 1. DEFINE SUBGROUPS
#------------------------------------------------------------------------------

subgroups <- list(
  age = list(var = "age", cuts = c(65), labels = c("<65", "≥65")),
  sex = list(var = "sex", cuts = NULL, labels = c("Female", "Male")),
  BMI = list(var = "BMI", cuts = c(25), labels = c("<25", "≥25")),
  UPDRS3 = list(var = "UPDRS3", cuts = c(33), labels = c("<33", "≥33")),
  genetic = list(var = "genetic_status", cuts = NULL, labels = c("Sporadic", "Genetic_carrier")),
  age_at_onset = list(var = "age_at_onset", cuts = c(60), labels = c("<60", "≥60"))
)

outcomes <- c("ESS", "REM", "GDS", "STAI_total", "STAI_state", 
              "STAI_trait", "SCOPA_AUT")

#------------------------------------------------------------------------------
# 2. RUN SUBGROUP ANALYSIS
#------------------------------------------------------------------------------

run_subgroup_analysis <- function(data, outcome, subgroup_def, covariates) {
  results <- data.frame()
  
  if (is.null(subgroup_def$cuts)) {
    # Categorical variable
    levels <- unique(data[[subgroup_def$var]])
    for (lv in levels) {
      subset_data <- data[data[[subgroup_def$var]] == lv, ]
      if (nrow(subset_data) > 10) {
        model <- lm(as.formula(paste(outcome, "~ motor_group +", 
                                     paste(covariates, collapse = " + "))), 
                    data = subset_data)
        coef <- extract_nonTD_coef(model)
        results <- rbind(results, data.frame(
          Subgroup = lv,
          N = nrow(subset_data),
          Beta = coef["beta"],
          CI_lower = coef["ci_lower"],
          CI_upper = coef["ci_upper"],
          P = coef["p"]
        ))
      }
    }
  } else {
    # Continuous variable with cut points
    for (i in seq_along(subgroup_def$cuts)) {
      cut_val <- subgroup_def$cuts[i]
      label <- subgroup_def$labels[i]
      
      if (i == 1) {
        subset_data <- data[data[[subgroup_def$var]] < cut_val, ]
      } else {
        subset_data <- data[data[[subgroup_def$var]] >= subgroup_def$cuts[i-1] & 
                              data[[subgroup_def$var]] < cut_val, ]
      }
      # Handle last interval
      if (i == length(subgroup_def$cuts)) {
        subset_data <- data[data[[subgroup_def$var]] >= cut_val, ]
      }
      
      if (nrow(subset_data) > 10) {
        model <- lm(as.formula(paste(outcome, "~ motor_group +", 
                                     paste(covariates, collapse = " + "))), 
                    data = subset_data)
        coef <- extract_nonTD_coef(model)
        results <- rbind(results, data.frame(
          Subgroup = label,
          N = nrow(subset_data),
          Beta = coef["beta"],
          CI_lower = coef["ci_lower"],
          CI_upper = coef["ci_upper"],
          P = coef["p"]
        ))
      }
    }
  }
  
  return(results)
}

# Run for all outcomes and subgroups
subgroup_results <- list()

for (outcome in outcomes) {
  cat("Analyzing", outcome, "subgroups...\n")
  outcome_results <- list()
  for (sg_name in names(subgroups)) {
    outcome_results[[sg_name]] <- run_subgroup_analysis(
      baseline_data, outcome, subgroups[[sg_name]], covariates
    )
    if (nrow(outcome_results[[sg_name]]) > 0) {
      outcome_results[[sg_name]]$Subgroup_var <- sg_name
    }
  }
  subgroup_results[[outcome]] <- outcome_results
}

#------------------------------------------------------------------------------
# 3. SAVE RESULTS FOR FIGURES
#------------------------------------------------------------------------------

forest_data <- list()
for (outcome in outcomes) {
  combined <- data.frame()
  for (sg_name in names(subgroups)) {
    if (!is.null(subgroup_results[[outcome]][[sg_name]]) && 
        nrow(subgroup_results[[outcome]][[sg_name]]) > 0) {
      combined <- rbind(combined, subgroup_results[[outcome]][[sg_name]])
    }
  }
  forest_data[[outcome]] <- combined
}

saveRDS(forest_data, "output/forest_data.rds")

cat("\nSubgroup analysis complete. Results saved for Figures 2-8.\n")

# Save session info
writeLines(capture.output(sessionInfo()), "output/session_info_03.txt")
