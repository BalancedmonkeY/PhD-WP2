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


#' @param outcome Outcome measure (RR or MD)
#' @param total_events Total number of events across the review
#' @param event_threshold_2 Threshold for number of events that will lead to downgrading 2 levels
#' @param event_threshold_1 Threshold for number of events that will lead to downgrading 1 level
#' @param CI_ub Upper bound of 95% CI of pooled effect
#' @param CI_lb Lower bound of 95% CI of pooled effect
#' @param CI_threshold Threshold to represent meaningful effect
#' @return levels = Number of levels the evidence is likely to be downgraded due to imprecision (0, 1, or 2)

imprecision_downgrades <- function(
    outcome,
    total_events = NULL,
    event_threshold_2 = 100,
    event_threshold_1 = 300,
    CI_lb,
    CI_ub,
    CI_threshold = NULL
) {
  
  # Obtain respective threshold value for the interpretation data frame
  if (!is.null(CI_threshold)) { # if user has specified
    threshold <- CI_threshold
  } else if (outcome == "RR") { # if dichotomous
    threshold <- 0.25
  } else { # if continuous
    threshold <- 0.5
  }
  
  # Obtain interpretation data frame
  if (outcome == "RR") {
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
  
  if (outcome == "RR") {
    
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
  if (outcome == "RR") {
    levels = max(event_downgrade, interpretation_downgrade)
  } else {
    levels = interpretation_downgrade
  }
  
  return(levels = levels)
}
