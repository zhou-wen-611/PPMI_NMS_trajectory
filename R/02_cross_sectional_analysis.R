###############################################################################
# 02_cross_sectional_analysis.R
# Cross-sectional analysis: Motor subtype vs. non-motor symptoms (Table 2)
# Author: Wen Zhou
# Date: 2026-07-24
# R version: 4.2.2
###############################################################################

library(dplyr)
library(mice)
library(lmtest)
library(sandwich)

source("R/functions.R")

# Load cleaned data
baseline_data <- readRDS("data/baseline_cleaned.rds")

#------------------------------------------------------------------------------
# 1. DEFINE VARIABLES
#------------------------------------------------------------------------------

outcomes <- c("ESS", "REM", "GDS", "STAI_total", "STAI_state", 
              "STAI_trait", "SCOPA_AUT", "MoCA")

# Covariates as specified in Methods
covariates <- c("age", "sex", "age_at_onset", "race", "BMI", 
                "education", "SBR", "HY_stage", "UPDRS3")

#------------------------------------------------------------------------------
# 2. MULTIPLE IMPUTATION
#------------------------------------------------------------------------------

# Variables to include in imputation model
vars_for_imputation <- c("PATNO", "motor_group", outcomes, covariates)
imputation_data <- baseline_data %>% select(all_of(vars_for_imputation))

# Perform multiple imputation (5 repetitions as specified in Methods)
set.seed(20260724)
imp <- mice(imputation_data, m = 5, method = "pmm", maxit = 20, printFlag = FALSE)

# Complete imputed datasets
imp_data_list <- complete(imp, action = "all")

#------------------------------------------------------------------------------
# 3. LINEAR REGRESSION (WITH POOLED RESULTS)
#------------------------------------------------------------------------------

run_linear_regression_pooled <- function(imp_list, outcome, covariates) {
  results <- list()
  
  for (i in 1:length(imp_list)) {
    data <- imp_list[[i]]
    # Ensure motor_group is factor with TD as reference
    data$motor_group <- factor(data$motor_group, levels = c("TD", "non-TD"))
    
    formula <- as.formula(paste(outcome, "~", paste(c("motor_group", covariates), 
                                                     collapse = " + ")))
    model <- lm(formula, data = data)
    results[[i]] <- model
  }
  
  # Pool results using Rubin's rules
  pooled <- pool(results)
  summary_pooled <- summary(pooled)
  
  # Extract non-TD coefficient
  td_idx <- grep("motor_groupnon-TD", rownames(summary_pooled))
  
  if (length(td_idx) > 0) {
    beta <- summary_pooled[td_idx, "estimate"]
    se <- summary_pooled[td_idx, "std.error"]
    p_val <- summary_pooled[td_idx, "p.value"]
    ci_lower <- beta - 1.96 * se
    ci_upper <- beta + 1.96 * se
    return(c(beta = beta, ci_lower = ci_lower, ci_upper = ci_upper, p = p_val))
  } else {
    return(c(beta = NA, ci_lower = NA, ci_upper = NA, p = NA))
  }
}

#------------------------------------------------------------------------------
# 4. RUN ANALYSIS FOR ALL OUTCOMES
#------------------------------------------------------------------------------

results_table <- data.frame(
  Outcome = outcomes,
  Beta = NA, CI_lower = NA, CI_upper = NA, P = NA
)

for (i in seq_along(outcomes)) {
  cat("Analyzing", outcomes[i], "...\n")
  res <- run_linear_regression_pooled(imp_data_list, outcomes[i], covariates)
  results_table[i, c("Beta", "CI_lower", "CI_upper", "P")] <- res
}

cat("\n=== Cross-sectional Analysis Results (Table 2) ===\n")
print(results_table)

# Save results
# write_csv(results_table, "output/table2_cross_sectional.csv")

# Save session info
writeLines(capture.output(sessionInfo()), "output/session_info_02.txt")
