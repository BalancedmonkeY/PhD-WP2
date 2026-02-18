# Function for describing whether the user-specified threshold has been met, and how #
Threshold_Description <- function(meta = NULL, est = NULL, ci_lb = NULL, ci_ub = NULL, pvalue = NULL, outcome = "RR",
                                  sig_level, zero = NA, lower = NA, upper = NA, est_pos = NA, est_neg = NA, new_or_og, SSnew) {
  

  # Where focus is on statistical significance (where alpha = 0.05)
  if (sig_level == 0.05 & !is.na(zero) & is.na(lower) & is.na(upper) & is.na(est_pos) & is.na(est_neg)) {
    meta_pvalue <- ifelse(is.null(pvalue), meta$pval, pvalue)
    if (meta_pvalue < 0.05) {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " give a statistically significant (alpha = 0.05) pooled estimate")
    } else {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " not give a statistically significant (alpha = 0.05) pooled estimate")
    }
    # Where focus is on statistical significance (where alpha != 0.05)
  } else if (sig_level != 0.05 & !is.na(zero) & is.na(lower) & is.na(upper) & is.na(est_pos) & is.na(est_neg)) {
    meta_lower <- ifelse(is.null(ci_lb), meta$ci.lb, ci_lb)
    meta_upper <- ifelse(is.null(ci_ub), meta$ci.ub, ci_ub)
    if (meta_upper < 0 | meta_lower > 0) {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " give a statistically significant (alpha = ", sig_level, ") pooled estimate")
    } else {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " not give a statistically significant (alpha = ", sig_level, ") pooled estimate")
    }
    # Where focus is on a specific upper or lower confidence band
  } else if (is.na(zero) & !is.na(lower) & is.na(upper) & is.na(est_pos) & is.na(est_neg)) {
    meta_lower <- ifelse(is.null(ci_lb), meta$ci.lb, ci_lb)
    if (meta_lower > lower) {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " give a ", round(100*(1 - sig_level),1), "% CI that is higher than ", ifelse(outcome %in% c('OR', 'RR'), exp(lower), lower))
    } else {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " not give a ", round(100*(1 - sig_level),1), "% CI that is higher than ", ifelse(outcome %in% c('OR', 'RR'), exp(lower), lower))
    }
  } else if (is.na(zero) & is.na(lower) & !is.na(upper) & is.na(est_pos) & is.na(est_neg)) {
    meta_upper <- ifelse(is.null(ci_ub), meta$ci.ub, ci_ub)
    if (meta_upper < upper) {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " give a ", round(100*(1 - sig_level),1), "% CI that is lower than ", ifelse(outcome %in% c('OR', 'RR'), exp(upper), upper))
    } else {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " not give a ", round(100*(1 - sig_level),1), "% CI that is lower than ", ifelse(outcome %in% c('OR', 'RR'), exp(upper), upper))
    }
    # Where focus is on clinical significance
  } else if (is.na(zero) & is.na(lower) & is.na(upper) & (!is.na(est_pos) | !is.na(est_neg))) {
    meta_est <- as.numeric(ifelse(is.null(est), meta$beta, est))
    if (!is.na(est_pos) & meta_est > est_pos) {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " give a pooled estimate that is higher than ", ifelse(outcome %in% c('OR', 'RR'), exp(est_pos), est_pos))
    } else if (!is.na(est_neg) & meta_est < est_neg) {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " give a pooled estimate that is lower than ", ifelse(outcome %in% c('OR', 'RR'), exp(est_neg), est_neg))
    } else if (!is.na(est_pos) & is.na(est_neg) & meta_est <= est_pos) {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " not give a pooled estimate that is higher than ", ifelse(outcome %in% c('OR', 'RR'), exp(est_pos), est_pos))
    } else if (!is.na(est_neg) & is.na(est_pos) & meta_est >= est_neg) {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " not give a pooled estimate that is lower than ", ifelse(outcome %in% c('OR', 'RR'), exp(est_neg), est_neg))
    } else if (!is.na(est_neg) & !is.na(est_pos) & meta_est <= est_pos & meta_est >= est_neg) {
      threshold_result <- paste0(ifelse(new_or_og == "new", paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will"), "Current studies do"), " not give a pooled estimate that is lower than ", ifelse(outcome %in% c('OR', 'RR'), exp(est_neg), est_neg), " or higher than ", ifelse(outcome %in% c('OR', 'RR'), exp(est_pos), est_pos))
    }
  } else {
    print ("(Only) one of zero, lower, upper, or est_pos/est_neg need to be given a value")
  }
}