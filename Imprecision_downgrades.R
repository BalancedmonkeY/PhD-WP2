#---------------------------------------------------------#
# Tools for estimating any downgrading due to imprecision #
#---------------------------------------------------------#

library(dplyr)

#' @param data Dataset containing information
#' @param event_cols Names of columns in data that refer to the number of events
#' @return total_events = total number of events across the review

sum_events <- function(
    data, 
    event_cols
) {
  
  # Create numerical table of all event columns
  event_data <- sapply(data %>% select(all_of(event_cols)), as.numeric )
  
  # Sum up columns and rows
  total_events = sum(colSums(event_data))
  
  return(total_events)
  
}


#' @param CI_ub Upper bound of 95% CI of pooled effect (on relative risk scale)
#' @param CI_lb Lower bound of 95% CI of pooled effect (on relative risk scale)
#' @param threshold Threshold to represent meaningful effect, in terms of a percentage change in relative risk
#' @return interpretation_df -> A dataframe of the interpretation, made up of three columns with logical indicators
CI_interpretation_bin <- function(
    CI_ub,
    CI_lb,
    threshold = 0.25
) {
  
  # Initialise dataframe
  interpretation_df <- data.frame(neg.effect = 0, no.effect = 0, pos.effect = 0)
  # Add 1 to columns where appropriate based on CI
  if (CI_lb < 1-threshold) {interpretation_df$neg.effect <- 1} 
  if (CI_ub > 1+threshold) {interpretation_df$pos.effect <- 1} 
  if ((CI_ub >= 1-threshold) & (CI_lb <= 1+threshold)) {interpretation_df$no.effect <- 1}
  
  return(interpretation_df)
}


#' @param CI_ub Upper bound of 95% CI of pooled effect (on mean difference scale)
#' @param CI_lb Lower bound of 95% CI of pooled effect (on mean difference scale)
#' @param threshold Threshold to represent meaningful effect, in terms of a mean difference (+/-)
#' @return interpretation_df -> A dataframe of the interpretation, made up of three columns with logical indicators
CI_interpretation_cont <- function(
    CI_ub,
    CI_lb,
    threshold = 0.5
) {
  
  # Initialise dataframe
  interpretation_df <- data.frame(neg.effect = 0, no.effect = 0, pos.effect = 0)
  # Add 1 to columns where appropriate based on CI
  if (CI_lb < 0-threshold) {interpretation_df$neg.effect <- 1} 
  if (CI_ub > 0+threshold) {interpretation_df$pos.effect <- 1} 
  if ((CI_ub >= 0-threshold) & (CI_lb <= 0+threshold)) {interpretation_df$no.effect <- 1}
  
  return(interpretation_df)
}


#' @param total_events Total number of events across the review
#' @param event_threshold_2 Threshold for number of events that will lead to downgrading 2 levels
#' @param event_threshold_1 Threshold for number of events that will lead to downgrading 1 level
#' @param CI_ub Upper bound of 95% CI of pooled effect
#' @param CI_lb Lower bound of 95% CI of pooled effect
#' @param CI_threshold Threshold to represent meaningful effect
#' @param auto_adjust Logical parameter to indicate whether to automatically adjust the thresholds based on previous versions of the review
#' @param prev_events Total number of events across the previous version
#' @param prev_CI_ub Upper bound of 95% CI of previous pooled effect
#' @param prev_CI_lb Lower bound of 95% CI of previous pooled effect
#' @param prev_interpretation_df Interpretation dataframe from previous version
#' @param prev_levels The number of levels the evidence was downgraded due to imprecision in the previous version
#' @return list containing: levels = Number of levels the evidence is likely to be downgraded due to imprecision (0, 1, or 2)
#'                          event_threshold_2 = threshold value used for number of events that will lead to downgrading 2 levels (may have changed from autoadjust)
#'                          event_threshold_1 = threshold value used for number of events that will lead to downgrading 1 level (may have changed from autoadjust)
#'                          CI_threshold =  threshold value used represent meaningful effect (may have changed from autoadjust)

imprecision_downgrades <- function(
    total_events = NULL,
    event_threshold_2 = 100,
    event_threshold_1 = 300,
    CI_lb,
    CI_ub,
    CI_threshold = NULL,
    auto_adjust = FALSE,
    prev_events = NULL,
    prev_CI_lb,
    prev_CI_ub,
    prev_interpretation_df,
    prev_levels
) {
  
  # Obtain respective threshold value for the interpretation data frame
  if (!is.null(CI_threshold)) { # if user has specified
    threshold <- CI_threshold
  } else if (!is.null(total_events) & !is.na(total_events)) { # if dichotomous
    threshold <- 0.25
  } else { # if continuous
    threshold <- 0.5
  }
  
  # Adjust thresholds if needed, based on previous version
  if (auto_adjust) {
    # Calculate levels from each measure for the previous version
    if (!is.null(prev_events) & !is.na(prev_events)) {
      if (prev_events < event_threshold_2) {
        prev_event_downgrade <- 2
      } else if (prev_events < event_threshold_1) {
        prev_event_downgrade <- 1
      } else {
        prev_event_downgrade <- 0
      }
      prev_event_downgrade <- "NULL"
    }
    if (sum(prev_interpretation_df) == 3) {
      prev_interpretation_downgrade <- 2
    } else if (sum(prev_interpretation_df) == 2) {
      prev_interpretation_downgrade <- 1
    } else {
      prev_interpretation_downgrade <- 0
    }
    # Using imprecision_autoadjust_matrix, look-up action needed
    imprecision_autoadjust_matrix <- read.csv("imprecision_autoadjust_matrix.csv", header = TRUE)
    idx <- which(imprecision_autoadjust_matrix$Level == prev_levels & 
                  imprecision_autoadjust_matrix$Events == prev_event_downgrade &
                  imprecision_autoadjust_matrix$CI == prev_interpretation_downgrade)
    # Apply corresponding action
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
    if (!is.null(prev_events) & !is.na(prev_events)) {
      null <- 1
    } else {
      null <- 0
    }
    if (imprecision_autoadjust_matrix$CI.threshold[idx] == "Increase" & prev_levels == 0) {
      threshold <- ceiling(max(null - prev_CI_lb, prev_CI_ub - null)*20)/20
    } else if (imprecision_autoadjust_matrix$CI.threshold[idx] == "Increase" & prev_levels == 1) {
      threshold <- ceiling(min(null - prev_CI_lb, prev_CI_ub - null)*20)/20
    } else if (imprecision_autoadjust_matrix$CI.threshold[idx] == "Decrease" & prev_levels == 1) {
      threshold <- floor((max(null - prev_CI_lb, prev_CI_ub - null)-0.01)*20)/20
    } else if (imprecision_autoadjust_matrix$CI.threshold[idx] == "Decrease" & prev_levels == 2) {
      threshold <- floor((min(null - prev_CI_lb, prev_CI_ub - null)-0.01)*20)/20
    }
    
  }
  
  # Obtain interpretation data frame
  if (!is.null(total_events) & !is.na(total_events)) {
    interpretation_df <- CI_interpretation_bin(
      CI_lb = CI_lb,
      CI_ub = CI_ub,
      threshold = threshold
    )
  } else {
    interpretation_df <- CI_interpretation_cont(
      CI_lb = CI_lb,
      CI_ub = CI_ub,
      threshold = threshold
    )
  }
  
  # If dichotomous assess downgrading levels due to small number of events
  
  if (!is.null(total_events) & !is.na(total_events)) {
    
    if (total_events < event_threshold_2) {
      event_downgrade <- 2
    } else if (total_events < event_threshold_1) {
      event_downgrade <- 1
    } else {
      event_downgrade <- 0
    }
  }
  
  # Assess downgrading levels due to interpretation coverage
  if (sum(interpretation_df) == 3) {
    interpretation_downgrade <- 2
  } else if (sum(interpretation_df) == 2) {
    interpretation_downgrade <- 1
  } else {
    interpretation_downgrade <- 0
  }
  
  # Take maximum number of levels downgraded of the two approaches (where present)
  if (!is.null(total_events) & !is.na(total_events)) {
    levels = max(event_downgrade, interpretation_downgrade)
  } else {
    levels = interpretation_downgrade
  }
  
  return(list(levels = levels,
              event_threshold_2 = event_threshold_2,
              event_threshold_1 = event_threshold_1,
              CI_threshold = CI_threshold))
}
