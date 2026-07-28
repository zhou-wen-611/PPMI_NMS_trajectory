###############################################################################
# 06_trajectory_logistic_regression.R
# Trajectory-motor subtype association (Tables 4-5)
# Author: Wen Zhou
# Date: 2026-07-24
# R version: 4.2.2
###############################################################################

library(dplyr)

source("R/functions.R")

# Load data
baseline_data <- readRDS("data/baseline_cleaned.rds")
class_assignments <- readRDS("output/class_assignments.rds")

#------------------------------------------------------------------------------
# 1. DEFINE UNFAVORABLE TRAJECTORIES (for binary outcomes)
#------------------------------------------------------------------------------

unfavorable_class <- list(
  ESS = 2,           # Rapid increase
  REM = 2,           # Rapid increase
  GDS = 3,           # High-level progressive increase
  STAI_total = 2,    # High-level progressive increase
  STAI_state = 2,    # High-level progressive increase
  STAI_trait = 2,    # High-level progressive increase
  SCOPA_AUT = 3,     # High-level rapid increase
  MoCA = 3           # Low-level progressive decline
)

#------------------------------------------------------------------------------
# 2. LOGISTIC REGRESSION FOR BINARY TRAJECTORIES (Table 4)
#------------------------------------------------------------------------------

run_binary_trajectory_logistic <- function(class_data, outcome, unfavorable_class, baseline) {
  merged <- merge(class_data, baseline, by = "PATNO")
  merged$unfavorable <- ifelse(merged$Class == unfavorable_class, 1, 0)
  merged$motor_group <- factor(merged$motor_group, levels = c("TD", "non-TD"))
  
  model <- glm(
    unfavorable ~ motor_group + age + sex + age_at_onset + race + BMI + 
      education + SBR + HY_stage + UPDRS3,
    data = merged,
    family = binomial()
  )
  
  summary_model <- summary(model)
  coef_idx <- grep("motor_groupnon-TD", rownames(summary_model$coefficients))
  
  if (length(coef_idx) > 0) {
    OR <- exp(summary_model$coefficients[coef_idx, "Estimate"])
    ci_lower <- exp(summary_model$coefficients[coef_idx, "Estimate"] - 
                      1.96 * summary_model$coefficients[coef_idx, "Std. Error"])
    ci_upper <- exp(summary_model$coefficients[coef_idx, "Estimate"] + 
                      1.96 * summary_model$coefficients[coef_idx, "Std. Error"])
    p_val <- summary_model$coefficients[coef_idx, "Pr(>|z|)"]
    
    return(data.frame(
      Outcome = outcome,
      N = nrow(merged),
      N_unfavorable = sum(merged$unfavorable),
      OR = OR,
      CI_lower = ci_lower,
      CI_upper = ci_upper,
      P = p_val
    ))
  } else {
    return(NULL)
  }
}

binary_logistic_results <- list()

for (outcome in names(class_assignments)) {
  if (outcome %in% names(unfavorable_class)) {
    cat("Running binary logistic regression for", outcome, "...\n")
    binary_logistic_results[[outcome]] <- run_binary_trajectory_logistic(
      class_assignments[[outcome]], 
      outcome, 
      unfavorable_class[[outcome]], 
      baseline_data
    )
  }
}

binary_results_table <- do.call(rbind, binary_logistic_results)

cat("\n=== Binary Trajectory Logistic Regression Results (Table 4) ===\n")
print(binary_results_table)

# write_csv(binary_results_table, "output/table4_binary_trajectory.csv")

#------------------------------------------------------------------------------
# 3. LOGISTIC REGRESSION FOR MULTI-CLASS TRAJECTORIES (Table 5)
#------------------------------------------------------------------------------

run_multi_class_trajectory_logistic <- function(class_data, outcome, baseline) {
  merged <- merge(class_data, baseline, by = "PATNO")
  merged$motor_group <- factor(merged$motor_group, levels = c("TD", "non-TD"))
  
  # Use multinomial logistic regression (nnet package required)
  # For simplicity, we use separate binary logistic regressions
  # comparing each higher class to the lowest (reference) class
  
  classes <- sort(unique(merged$Class))
  results <- data.frame()
  
  for (i in 2:length(classes)) {
    ref_class <- classes[1]
    target_class <- classes[i]
    
    # Subset to only these two classes
    subset_data <- merged[merged$Class %in% c(ref_class, target_class), ]
    subset_data$binary_class <- ifelse(subset_data$Class == target_class, 1, 0)
    
    model <- glm(
      binary_class ~ motor_group + age + sex + age_at_onset + race + BMI + 
        education + SBR + HY_stage + UPDRS3,
      data = subset_data,
      family = binomial()
    )
    
    summary_model <- summary(model)
    coef_idx <- grep("motor_groupnon-TD", rownames(summary_model$coefficients))
    
    if (length(coef_idx) > 0) {
      OR <- exp(summary_model$coefficients[coef_idx, "Estimate"])
      ci_lower <- exp(summary_model$coefficients[coef_idx, "Estimate"] - 
                        1.96 * summary_model$coefficients[coef_idx, "Std. Error"])
      ci_upper <- exp(summary_model$coefficients[coef_idx, "Estimate"] + 
                        1.96 * summary_model$coefficients[coef_idx, "Std. Error"])
      p_val <- summary_model$coefficients[coef_idx, "Pr(>|z|)"]
      
      results <- rbind(results, data.frame(
        Outcome = outcome,
        Comparison = paste0("Class ", target_class, " vs Class ", ref_class),
        N = nrow(subset_data),
        OR = OR,
        CI_lower = ci_lower,
        CI_upper = ci_upper,
        P = p_val
      ))
    }
  }
  
  return(results)
}

multi_class_results <- list()

for (outcome in names(class_assignments)) {
  if (outcome %in% c("GDS", "SCOPA_AUT", "MoCA")) {
    cat("Running multi-class logistic regression for", outcome, "...\n")
    multi_class_results[[outcome]] <- run_multi_class_trajectory_logistic(
      class_assignments[[outcome]], outcome, baseline_data
    )
  }
}

multi_class_results_table <- do.call(rbind, multi_class_results)

cat("\n=== Multi-Class Trajectory Logistic Regression Results (Table 5) ===\n")
print(multi_class_results_table)

# write_csv(multi_class_results_table, "output/table5_multi_class_trajectory.csv")

# Save session info
writeLines(capture.output(sessionInfo()), "output/session_info_06.txt")
