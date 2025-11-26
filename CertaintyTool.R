#------------------------------------------------------------------------------------------------#
# Tool to ascertain if new studies for a meta-analysis will change the certainty of the evidence #
#------------------------------------------------------------------------------------------------#
#-------------------------------#
# Clareece Nevill November 2025 #
#-------------------------------#

library(dplyr)
library(metafor)
source("RoB_downgrades.R")
source("Inconsistency_downgrades.R")
source("Imprecision_downgrades.R")
source("Publication_bias_downgrades.R")

#' @param data Data frame where each row contains a study, and the columns include: (i) RoB ratings, (ii) estimates & CI bounds & vi (in log form if ratio), (iii) Weight, (iv) Number of events (if dichotomous)
#' @param overall_rob Overall rating of risk of bias
#' @param random_selection Random sequence generation (selection bias) rating (RoB 1 only)
#' @param allocation_selection Allocation concealment (selection bias) rating (RoB 1 only)
#' @param performance Performance bias rating (RoB 1 only)
#' @param detection Detection bias rating (RoB 1 only)
#' @param attrition Attrition bias rating (RoB 1 only)
#' @param reporting Reporting bias rating (RoB 1 only)
#' @param other Other bias rating (RoB 1 only)
#' @param weights Meta-analysis weighting for each study
#' @param event_cols Column names for events (if dichotomous)
#' @param CI_lb_col Column name that refer to the CI lower bound (in log-form for ratios)
#' @param CI_ub_col Column name that refers to the CI upper bound (in log forms for ratios)
#' @param estimates Column name that refers to the point estimates of all studies (in log forms for ratios)
#' @param variances Column name that refers to the variances of each study
#' @param rob_tool What type of risk of bias tool was used (1 or 2)
#' @param outcome Outcome measure
#' @param null_effect Value for which to assess the point estimates against
#' @param discrepancies Logical argument for whether there exist discrepancies between published and unpublished data
#' @param industry Logical argument for whether there exists industry influence
#' @param search Logical argument for whether there are concerns regarding the search integrity
#' @param RoB_1_threshold Threshold for rating down by 1 from RoB
#' @param RoB_2_threshold Threshold for rating down by 2 from RoB
#' @param events_1_threshold Threshold for rating down by 1 from total number of events
#' @param events_2_threshold Threshold for rating down by 2 from total number of events
#' @param CI_threshold Threshold for meaningful effect to then downgrade for imprecision
#' @param Jaccard_threshold Threshold for downgrading by 1 from Jaccard index
#' @param variation_threshold Threshold for downgrading by 1 from variation in estimates
#' @param Eggers_threshold Threshold for downgrading by 1 from Eggers test
#' @return List containing:
#' result - predicted GRADE score;
#' RoB - predicted levels downgraded due to RoB
#' Imprecision - predicted levels downgraded due to imprecision
#' Inconsistency - predicted levels downgraded due to inconsistency
#' Pubbias - predicted levels downgraded due to publication bias

PredictedScore <- function(
    data,
    overall_rob,
    random_selection,
    allocation_selection,
    performance,
    detection,
    attrition,
    reporting,
    other,
    weights,
    event_cols,
    CI_lb_col,
    CI_ub_col,
    estimates,
    variances,
    rob_tool,
    outcome,
    null_effect = 0,
    discrepancies = FALSE,
    industry = FALSE,
    search = FALSE,
    RoB_1_threshold = 1.5,
    RoB_2_threshold = NULL,
    events_1_threshold = 300,
    events_2_threshold = 100,
    CI_threshold = NULL,
    Jaccard_threshold = 0.4,
    variation_threshold = 0.8,
    Eggers_threshold = 0.9
) {
  
}

