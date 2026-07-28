###############################################################################
# 04_sensitivity_analysis.R
# Sensitivity analyses: Complete-case analysis
# Author: Wen Zhou
# Date: 2026-07-24
# R version: 4.2.2
###############################################################################

library(dplyr)

source("R/functions.R")

# Load cleaned data
baseline_data <- readRDS("data/baseline_cleaned.rds")

outcomes <- c("ESS", "REM", "GDS", "STAI_total", "STAI_state", 
              "STAI_trait", "SCOPA_AUT", "MoCA")
covariates <- c("age", "sex", "age_at_onset", "race", "BMI", 
                "education", "SBR", "HY_stage", "UPDRS3")

#------------------------------------------------------------------------------
# 1. COMPLETE CASE ANALYSIS
#------------------------------------------------------------------------------

# Exclude observations with missing data (complete-case analysis)
complete_data <- baseline_data[complete.cases(baseline_data[, c("motor_group", covariates)]), ]
cat("Complete cases:", nrow(complete_data), "of", nrow(baseline_data), "\n")

sensitivity_results <- data.frame(
  Outcome = outcomes,
  Complete_Beta = NA,
  Complete_CI_lower = NA,
  Complete_CI_upper = NA,
  Complete_P = NA
)

for (i in seq_along(outcomes)) {
  # Complete case model
  model <- lm(as.formula(paste(outcomes[i], "~ motor_group +", 
                               paste(covariates, collapse = " + "))), 
              data = complete_data)
  coef <- extract_nonTD_coef(model)
  sensitivity_results[i, c("Complete_Beta", "Complete_CI_lower", 
                           "Complete_CI_upper", "Complete_P")] <- coef
}

cat("\n=== Sensitivity Analysis Results ===\n")
print(sensitivity_results)

# Save results
# write_csv(sensitivity_results, "output/sensitivity_results.csv")

# Save session info
writeLines(capture.output(sessionInfo()), "output/session_info_04.txt")
