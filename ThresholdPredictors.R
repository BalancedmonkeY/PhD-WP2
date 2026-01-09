#-------------------------------------------------------#
# Tools for estimating thresholds used by other reviews #
#-------------------------------------------------------#

#--------------#
# Risk of Bias #
#--------------#

#' @param prev_RoB_avg Weighted average risk of bias score from the previous version
#' @param prev_levels The number of levels the evidence was downgraded due to risk of bias in the previous version
#' @return list containing: suggested_threshold_1 - The suggested threshold for downgrading one level
#'                          suggested_threshold_2 - The suggested threshold for downgrading two levels

RoB_threshold_finder <- function(
    prev_RoB_avg,
    prev_levels
) {
  
  suggested_threshold_1 = 1.5
  suggested_threshold_2 = 2.5

  # increase thresholds if not previously downgraded, but previous score was 1.5 or higher
  if (prev_levels == 0 & prev_RoB_avg >= 1.5) {
    suggested_threshold_1 <- ceiling((prev_RoB_avg+0.01)*20)/20 # puts threshold at next highest value, on 0.05 scale
    if (prev_RoB_avg >= 2.5) {
      suggested_threshold_2 <- ceiling((prev_RoB_avg+((3-prev_RoB_avg)/2))*20)/20  # moves threshold 2 to in-between
    }
  }
  
  # decrease one_level_threshold if previously downgraded 1 level, but previous score was lower than 1.5
  else if (prev_levels == 1 & prev_RoB_avg < 1.5) {
    suggested_threshold_1 <- floor(prev_RoB_avg*20)/20
  }
  
  # increase two_level_threshold if previously downgraded 1 level, but previous score was 2.5 or higher
  else if (prev_levels == 1 & prev_RoB_avg >= 2.5) {
    suggested_threshold_2 <- ceiling((prev_RoB_avg+0.01)*20)/20 # puts threshold at next highest value, on 0.05 scale
  }
  
  # decrease threshold if previously downgraded 2 levels, but previous score was less than 2.5
  else if (prev_levels == 2 & prev_RoB_avg < 2.5) {
    suggested_threshold_2 <- floor(prev_RoB_avg*20)/20
    if (prev_RoB_avg < 1.5) {
      suggested_threshold_1 <- floor((prev_RoB_avg-((prev_RoB_avg-1)/2))*20)/20 # moves threshold 1 to in-between
    }
  }
  
  return(list(suggested_threshold_1 = suggested_threshold_1,
              suggested_threshold_2 = suggested_threshold_2))
  
}

#------------------#
# Publication bias #
#------------------#

#' @param prev_stat Funnel plot statistics from previous update
#' @param prev_industry Presence of industry influence from previous update
#' @param prev_search Presence of non-comprehensive search from previous update
#' @param prev_levels The number of levels the evidence was downgraded due to publication bias in the previous version
#' @return suggested_threshold The suggested threshold for downgrading one level

Pubbias_threshold_finder <- function(
    prev_stat,
    prev_industry,
    prev_search,
    prev_levels
) {
  
  suggested_threshold = NULL
  
  # decrease threshold if not previously downgraded, but p-value was less than 0.1
  if (prev_levels == 0 & prev_stat < 0.1 & prev_industry == FALSE & prev_search == FALSE) {
    suggested_threshold <- floor((prev_stat)*20)/20
  } 
  
  # increase threshold if previously downgraded, but p-value was 0.1 or higher
  else if (prev_levels == 1 & prev_stat < threshold & prev_industry == FALSE & prev_search == FALSE) {
    suggested_threshold <- ceiling((prev_stat+0.01)*20)/20
  }
  
  if (!is.null(suggested_threshold)) {
    return(suggested_threshold = suggested_threshold)
  } else {
    message("No new threshold is needed")
  }
  
}

#---------------#
# Inconsistency #
#---------------#

#' @param prev_levels The number of levels the evidence was downgraded due to inconsistency in the previous version
#' @param prev_Jaccard Jaccard index for previous update
#' @param prev_est_var Variation in estimates score for previous update
#' @return list containing: suggested_variation = The suggested threshold for variation in estimates
#'                          suggested_Jaccard = The suggested threshold for the Jaccard index

Inconsistency_threshold_finder <- function(
    prev_Jaccard,
    prev_est_var,
    prev_levels
) {
  
  suggested_variation = NULL
  suggested_Jaccard = NULL

  # increase threshold if previously downgraded, but est_var was 0.8 or higher
  if (prev_levels == 1) {
    if (prev_est_var >= 0.8) {
      suggested_variation <- ceiling((prev_est_var+0.01)*20)/20
      if (suggested_variation >= 1) { # may happen if the previous est_var equaled 1 exactly
        suggested_variation = 1
      }
      if (prev_est_var >= 0.95) { # i.e., threshold will become 1 which will never lead to a downgrading
        
        # increase threshold if previously downgraded, but variation threshold has become 1, and Jaccard was 0.4 or higher
        if (prev_Jaccard >= 0.4) {
          suggested_Jaccard <- ceiling((prev_Jaccard+0.01)*20)/20
        }
      }
      
    # increase threshold if previously downgraded, est_var is less than 0.8, but Jaccard is 0.4 or higher  
    } else {
      if (prev_Jaccard >= 0.4) {
        suggested_Jaccard <- ceiling((prev_Jaccard+0.01)*20)/20
      }
    }
    
  # decrease threshold if didn't previously downgrade, but est_var was less than 0.8  
  } else if (prev_levels == 0) {
    if (prev_est_var < 0.8) {
      suggested_variation <- floor((prev_est_var)*20)/20
      if (prev_est_var < 0.55) { # i.e., threshold will become 0.5 which is is the lowest possible and indicative of checking the Jaccard
        
        # decrease threshold if didn't previously downgrade, but variation threshold have become 0.5, and Jaccard was lower than 0.4
        if (prev_Jaccard < 0.4) {
          suggested_Jaccard <- floor((prev_Jaccard)*20)/20
        }
      } 
    }
  }
  
  if (!is.null(suggested_variation) | !is.null(suggested_Jaccard)) {
    return(list(
      suggested_variation = suggested_variation,
      suggested_Jaccard = suggested_Jaccard
    ))
  } else {
    message("No new thresholds are needed")
  }
  
}

#-------------#
# Imprecision #
#-------------#

#' @param outcome Outcome measure (RR or MD)
#' @param prev_events Total number of events across the previous version
#' @param prev_CI_ub Upper bound of 95% CI of previous pooled effect (log scale if ratio)
#' @param prev_CI_lb Lower bound of 95% CI of previous pooled effect (log scale if ratio)
#' @param prev_levels The number of levels the evidence was downgraded due to imprecision in the previous version
#' @return list containing: suggested_event_threshold_1 = The suggested threshold for number of events (1 level)
#'                          suggested_event_threshold_2 = The suggested threshold for number of events (2 levels)
#'                          suggested_CI_threshold_pos = The suggested threshold for forming meaningful +ve effect region to assess CI width
#'                          suggested_CI_threshold_neg = The suggested threshold for forming meaningful -ve effect region to assess CI width

Imprecision_threshold_finder <- function(
    outcome,
    prev_events,
    prev_CI_lb,
    prev_CI_ub,
    prev_levels
) {
  
  suggested_event_threshold_1 = NULL
  suggested_event_threshold_2 = NULL
  suggested_CI_threshold_pos = NULL
  suggested_CI_threshold_neg = NULL

  # Obtain downgrading levels due to each measure using default values
  # Total number of events
  if (outcome == "RR") {
    if (prev_events < 100) {
      prev_event_downgrade <- 2
    } else if (prev_events < 300) {
      prev_event_downgrade <- 1
    } else {
      prev_event_downgrade <- 0
    }
  } else {
    prev_event_downgrade <- "NULL"
  }
  # CI width
  prev_interpretation_df <- CI_interpretation(
    CI_lb = prev_CI_lb,
    CI_ub = prev_CI_ub,
    outcome = outcome
  )
  prev_interpretation_downgrade <- sum(prev_interpretation_df) - 1
  
  # Using imprecision_autoadjust_matrix, look-up action needed
  imprecision_autoadjust_matrix <- read.csv("imprecision_autoadjust_matrix.csv", header = TRUE)
  idx <- which(imprecision_autoadjust_matrix$Level == prev_levels & 
                 imprecision_autoadjust_matrix$Events == prev_event_downgrade &
                 imprecision_autoadjust_matrix$CI == prev_interpretation_downgrade)
  
  # Apply corresponding action (if needed)
  if (length(idx) > 0) {
    # Event threshold 1
    if (imprecision_autoadjust_matrix$Event.threshold.1[idx] == "Decrease") {
      event_threshold_1 <- prev_events
      event_threshold_2 <- round(prev_events/2, 1)
    } 
    # Event threshold 2
    if (imprecision_autoadjust_matrix$Event.threshold.2[idx] == "Increase") {
      event_threshold_2 <- prev_events + 1
    } else if (imprecision_autoadjust_matrix$Event.threshold.2[idx] == "Decrease") {
      event_threshold_2 <- prev_events
    }
    # CI threshold
    if (imprecision_autoadjust_matrix$CI.threshold[idx] == "Increase" & prev_levels == 0) {
      threshold <- ceiling(max(-prev_CI_lb, prev_CI_ub)*20)/20
    } else if (imprecision_autoadjust_matrix$CI.threshold[idx] == "Increase" & prev_levels == 1) {
      threshold <- ceiling(min(-prev_CI_lb, prev_CI_ub)*20)/20
    } else if (imprecision_autoadjust_matrix$CI.threshold[idx] == "Decrease" & prev_levels == 1) {
      threshold <- floor((max(-prev_CI_lb, prev_CI_ub)-0.01)*20)/20
    } else if (imprecision_autoadjust_matrix$CI.threshold[idx] == "Decrease" & prev_levels == 2) {
      threshold <- floor((min(-prev_CI_lb, prev_CI_ub)-0.01)*20)/20
    }
  }
  
  if (!is.null(suggested_event_threshold_1) | !is.null(suggested_event_threshold_2) | !is.null(suggested_CI_threshold_pos) | !is.null(suggested_CI_threshold_neg)) {
    return(list(suggested_event_threshold_1 = event_threshold_1,
                suggested_event_threshold_2 = event_threshold_2,
                suggested_CI_threshold_pos = threshold,
                suggested_CI_threshold_neg = -threshold))
  } else {
    message("No new thresholds are needed")
  }
  
  
  
}

  

