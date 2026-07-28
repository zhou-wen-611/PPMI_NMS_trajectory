###############################################################################
# 01_data_preparation.R
# Data cleaning and preparation for PPMI analysis
# Author: Wen Zhou
# Date: 2026-07-24
# R version: 4.2.2
###############################################################################

library(dplyr)
library(tidyr)
library(readr)

source("R/functions.R")

set.seed(20260724)

#------------------------------------------------------------------------------
# 1. SIMULATE DATA (FOR TESTING)
#------------------------------------------------------------------------------
# Note: In actual analysis, replace with PPMI data downloaded from:
# https://www.ppmi-info.org/access-data-specimens/download-data

cat("Generating simulated PPMI data for testing...\n")
sim_data <- simulate_ppmi_data(n = 1384)

baseline_data <- sim_data$baseline
longitudinal_data <- sim_data$longitudinal

cat("Baseline data:", nrow(baseline_data), "participants\n")
cat("Longitudinal data:", nrow(longitudinal_data), "observations\n")

#------------------------------------------------------------------------------
# 2. LOAD ACTUAL PPMI DATA (UNCOMMENT WHEN DATA ARE AVAILABLE)
#------------------------------------------------------------------------------
# baseline_data <- read_csv("data/ppmi_baseline.csv")
# longitudinal_data <- read_csv("data/ppmi_longitudinal.csv")

#------------------------------------------------------------------------------
# 3. DATA CLEANING
#------------------------------------------------------------------------------

# Check for missing data
missing_baseline <- colSums(is.na(baseline_data))
cat("\nMissing data in baseline:\n")
print(missing_baseline[missing_baseline > 0])

# Define motor subtype groups (TD vs non-TD)
baseline_data <- baseline_data %>%
  mutate(
    motor_group = case_when(
      motor_subtype == "TD" ~ "TD",
      motor_subtype %in% c("PIGD", "indeterminate") ~ "non-TD",
      TRUE ~ NA_character_
    )
  )

# Remove participants with missing motor subtype
baseline_data <- baseline_data %>%
  filter(!is.na(motor_group))

cat("\nMotor subtype distribution:\n")
print(table(baseline_data$motor_group))

#------------------------------------------------------------------------------
# 4. CREATE LONGITUDINAL DATASET FOR TRAJECTORY ANALYSIS
#------------------------------------------------------------------------------

# Identify participants with follow-up data
participants_with_fu <- unique(longitudinal_data$PATNO)
baseline_with_fu <- baseline_data %>%
  filter(PATNO %in% participants_with_fu)

cat("\nParticipants with follow-up:", nrow(baseline_with_fu), "\n")
cat("Participants excluded (no follow-up):", 
    nrow(baseline_data) - nrow(baseline_with_fu), "\n")

# For trajectory modeling, require at least baseline + 1 follow-up
# (already satisfied by having any longitudinal data)

#------------------------------------------------------------------------------
# 5. EXPORT CLEANED DATA
#------------------------------------------------------------------------------

# Uncomment to export
# write_csv(baseline_data, "data/baseline_cleaned.csv")
# write_csv(baseline_with_fu, "data/baseline_longitudinal.csv")
# write_csv(longitudinal_data, "data/longitudinal_cleaned.csv")

cat("\nData preparation complete.\n")
cat("Total baseline sample:", nrow(baseline_data), "\n")
cat("Longitudinal sample:", nrow(baseline_with_fu), "\n")
cat("Total observations:", nrow(longitudinal_data), "\n")

# Save session info
writeLines(capture.output(sessionInfo()), "output/session_info_01.txt")
