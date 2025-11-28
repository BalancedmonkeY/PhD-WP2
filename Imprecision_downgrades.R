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


#' @param CI_ub Upper bound of 95% CI of pooled effect (log scale if ratio)
#' @param CI_lb Lower bound of 95% CI of pooled effect (log scale if ratio)
#' @param threshold_pos Threshold to represent meaningful positive effect (log scale if ratio)
#' @param threshold_neg Threshold to represent meaningful negative effect (log scale if ratio)
#' @param outcome Type of outcome measure (RR or MD)
#' @return interpretation_df -> A dataframe of the interpretation, made up of three columns with logical indicators
CI_interpretation <- function(
    CI_ub,
    CI_lb,
    threshold_pos = NULL,
    threshold_neg = NULL,
    outcome
) {
  
  # Set defaults
  if (is.null(threshold_pos)) {
    if (outcome == "RR") {
      threshold_pos <- log(1.25)
    } else {
      threshold_pos <- 0.5
    }
  } 
  if (is.null(threshold_neg)) {
    if (outcome == "RR") {
      threshold_neg <- log(0.75)
    } else {
      threshold_neg <- -0.5
    }
  }

  
  # Initialise dataframe
  interpretation_df <- data.frame(neg.effect = 0, no.effect = 0, pos.effect = 0)
  
  # Add 1 to columns where appropriate based on CI
  if (CI_lb < threshold_neg) {interpretation_df$neg.effect <- 1} 
  if (CI_ub > threshold_pos) {interpretation_df$pos.effect <- 1} 
  if ((CI_ub >= threshold_neg) & (CI_lb <= threshold_pos)) {interpretation_df$no.effect <- 1}
  
  return(interpretation_df)
}


#' @param outcome Outcome measure (RR or MD)
#' @param total_events Total number of events across the review
#' @param event_threshold_2 Threshold for number of events that will lead to downgrading 2 levels
#' @param event_threshold_1 Threshold for number of events that will lead to downgrading 1 level
#' @param CI_ub Upper bound of 95% CI of pooled effect (log scale if ratio)
#' @param CI_lb Lower bound of 95% CI of pooled effect (log scale if ratio)
#' @param CI_threshold_pos Threshold to represent meaningful +ve effect
#' @param CI_threshold_neg Threshold to represent meaningful -ve effect
#' @return levels = Number of levels the evidence is likely to be downgraded due to imprecision (0, 1, or 2)

imprecision_downgrades <- function(
    outcome,
    total_events = NULL,
    event_threshold_2 = 100,
    event_threshold_1 = 300,
    CI_lb,
    CI_ub,
    CI_threshold_pos = NULL,
    CI_threshold_neg = NULL
) {
  
  # Obtain interpretation data frame
  interpretation_df <- CI_interpretation(
    CI_lb = CI_lb,
    CI_ub = CI_ub,
    threshold_pos = CI_threshold_pos,
    threshold_neg = CI_threshold_neg,
    outcome = outcome
  )
  
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
  interpretation_downgrade <- sum(interpretation_df) - 1
  
  # Take maximum number of levels downgraded of the two approaches (where present)
  if (outcome == "RR") {
    levels = max(event_downgrade, interpretation_downgrade)
  } else {
    levels = interpretation_downgrade
  }
  
  return(levels = levels)
}
