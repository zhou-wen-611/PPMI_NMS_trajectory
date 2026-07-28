###############################################################################
# R/functions.R
# Custom helper functions for the PPMI analysis
# Author: Wen Zhou
# Date: 2026-07-24
# R version: 4.2.2
###############################################################################

#' Simulate PPMI-like data for testing
#'
#' Creates synthetic data similar to PPMI structure for testing code without
#' accessing actual PPMI data. All values are randomly generated and do not
#' represent real patient data.
#'
#' @param n Number of participants to simulate
#' @return List containing baseline and longitudinal data frames
simulate_ppmi_data <- function(n = 1384) {
  set.seed(20260724)
  
  # Baseline data
  data_baseline <- data.frame(
    PATNO = paste0("P", sprintf("%04d", 1:n)),
    age = round(rnorm(n, 65, 10), 1),
    sex = sample(c("Male", "Female"), n, replace = TRUE, prob = c(0.6, 0.4)),
    BMI = round(rnorm(n, 26, 4), 1),
    education = round(rnorm(n, 15, 3), 1),
    age_at_onset = round(rnorm(n, 58, 12), 1),
    race = sample(c("White", "Black", "Asian", "Other"), n, replace = TRUE, 
                  prob = c(0.85, 0.05, 0.05, 0.05)),
    SBR = round(runif(n, 0.5, 3), 2),
    HY_stage = sample(1:3, n, replace = TRUE, prob = c(0.4, 0.4, 0.2)),
    UPDRS3 = round(rnorm(n, 28, 12), 1),
    motor_subtype = sample(c("TD", "non-TD"), n, replace = TRUE, 
                           prob = c(0.67, 0.33)),
    genetic_status = sample(c("Sporadic", "Genetic_carrier"), n, 
                            replace = TRUE, prob = c(0.8, 0.2))
  )
  
  # Baseline non-motor symptoms
  data_baseline$ESS <- round(rnorm(n, 8, 4), 1)
  data_baseline$REM <- round(rnorm(n, 5, 3), 1)
  data_baseline$GDS <- round(rnorm(n, 4, 3), 1)
  data_baseline$STAI_total <- round(rnorm(n, 40, 10), 1)
  data_baseline$STAI_state <- round(rnorm(n, 35, 8), 1)
  data_baseline$STAI_trait <- round(rnorm(n, 38, 9), 1)
  data_baseline$SCOPA_AUT <- round(rnorm(n, 12, 6), 1)
  data_baseline$MoCA <- round(rnorm(n, 26, 4), 1)
  
  # Longitudinal data
  n_visits <- sample(3:10, n, replace = TRUE, 
                     prob = c(0.1, 0.15, 0.2, 0.2, 0.15, 0.1, 0.05, 0.05))
  
  data_long <- list()
  for (i in 1:n) {
    visits <- n_visits[i]
    times <- sort(runif(visits, 0.5, 14))
    for (j in 1:visits) {
      data_long[[length(data_long) + 1]] <- data.frame(
        PATNO = data_baseline$PATNO[i],
        visit_time = times[j],
        ESS = round(data_baseline$ESS[i] + rnorm(1, 0.2 * times[j], 1), 1),
        REM = round(data_baseline$REM[i] + rnorm(1, 0.1 * times[j], 0.8), 1),
        GDS = round(data_baseline$GDS[i] + rnorm(1, 0.15 * times[j], 0.8), 1),
        STAI_total = round(data_baseline$STAI_total[i] + rnorm(1, 0.3 * times[j], 2), 1),
        STAI_state = round(data_baseline$STAI_state[i] + rnorm(1, 0.2 * times[j], 1.5), 1),
        STAI_trait = round(data_baseline$STAI_trait[i] + rnorm(1, 0.25 * times[j], 1.5), 1),
        SCOPA_AUT = round(data_baseline$SCOPA_AUT[i] + rnorm(1, 0.2 * times[j], 1), 1),
        MoCA = round(data_baseline$MoCA[i] + rnorm(1, -0.1 * times[j], 0.8), 1)
      )
    }
  }
  data_long <- do.call(rbind, data_long)
  
  return(list(baseline = data_baseline, longitudinal = data_long))
}

#' Extract non-TD coefficient from linear regression
#'
#' @param model A linear regression model object (lm)
#' @return Vector with beta, CI lower, CI upper, and p-value
extract_nonTD_coef <- function(model) {
  coef_summary <- summary(model)$coefficients
  td_idx <- grep("motor_subtypenon-TD", rownames(coef_summary))
  
  if (length(td_idx) > 0) {
    beta <- coef_summary[td_idx, "Estimate"]
    se <- coef_summary[td_idx, "Std. Error"]
    p_val <- coef_summary[td_idx, "Pr(>|t|)"]
    ci_lower <- beta - 1.96 * se
    ci_upper <- beta + 1.96 * se
    return(c(beta = beta, ci_lower = ci_lower, ci_upper = ci_upper, p = p_val))
  } else {
    return(c(beta = NA, ci_lower = NA, ci_upper = NA, p = NA))
  }
}

#' Calculate entropy for latent class trajectory models
#'
#' @param post_prob Matrix of posterior probabilities (rows = subjects, columns = classes)
#' @return Entropy value
calculate_entropy <- function(post_prob) {
  n <- nrow(post_prob)
  entropy <- 0
  for (i in 1:n) {
    pi <- post_prob[i, ]
    entropy <- entropy + sum(-pi * log(pi))
  }
  return(1 - entropy / (n * log(ncol(post_prob))))
}

#' Calculate average posterior probability for each class
#'
#' @param post_prob Matrix of posterior probabilities
#' @param class_assign Vector of class assignments
#' @return Vector of APP values per class
calculate_app <- function(post_prob, class_assign) {
  n_classes <- ncol(post_prob)
  app <- numeric(n_classes)
  for (k in 1:n_classes) {
    idx <- which(class_assign == k)
    if (length(idx) > 0) {
      app[k] <- mean(post_prob[idx, k])
    } else {
      app[k] <- NA
    }
  }
  return(app)
}
