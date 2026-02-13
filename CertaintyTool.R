#------------------------------------------------------------------------------------------------#
# Tool to ascertain if new studies for a meta-analysis will change the certainty of the evidence #
#------------------------------------------------------------------------------------------------#
#-------------------------------#
# Clareece Nevill November 2025 #
#-------------------------------#

library(dplyr)
library(metafor)
source("../3. Create tool/RoB_downgrades.R")
source("../3. Create tool/Inconsistency_downgrades.R")
source("../3. Create tool/Imprecision_downgrades.R")
source("../3. Create tool/Publication_bias_downgrades.R")

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
#' @param events_trt Column name that refers to the number of events in the treatment arm
#' @param events_ctrl Column name that refers to the number of events in the control arm
#' @param n_trt Column name that refers to the total number of people in the treatment arm
#' @param n_ctrl Column name that refers to the total number of people in the control arm
#' @param rob_tool What type of risk of bias tool was used (1 or 2)
#' @param outcome Outcome measure
#' @param model Meta-analysis model (as per metafor options (EE, MH, or DL))
#' @param ma rma object containing meta-analysis results (optional)
#' @param pubbias_min_studies Minimum number of studies needed to calculate the funnel asymetry statistic
#' @param null_effect Value for which to assess the point estimates against
#' @param industry Logical argument for whether there exists industry influence
#' @param search Logical argument for whether there are concerns regarding the search integrity
#' @param RoB_1_threshold Threshold for rating down by 1 from RoB
#' @param RoB_2_threshold Threshold for rating down by 2 from RoB
#' @param events_1_threshold Threshold for rating down by 1 from total number of events
#' @param events_2_threshold Threshold for rating down by 2 from total number of events
#' @param CI_threshold_pos Threshold for meaningful +ve effect to then downgrade for imprecision
#' @param CI_threshold_neg Threshold for meaningful -ve effect to then downgrade for imprecision
#' @param Jaccard_threshold Threshold for downgrading by 1 from Jaccard index
#' @param variation_threshold Threshold for downgrading by 1 from variation in estimates
#' @param Eggers_threshold Threshold for downgrading by 1 from Eggers test
#' @param indirectness Number of levels to downgrade due to indirectness
#' @return List containing:
#' result - predicted GRADE score;
#' RoB - predicted levels downgraded due to RoB
#' Imprecision - predicted levels downgraded due to imprecision
#' Inconsistency - predicted levels downgraded due to inconsistency
#' Pubbias - predicted levels downgraded due to publication bias
#' Indirectness - predicted levels downgraded due to indirectness

PredictedGRADEdomains <- function(
    data,
    overall_rob,
    random_selection = NA,
    allocation_selection = NA,
    performance = NA,
    detection = NA,
    attrition = NA,
    reporting = NA,
    other = NA,
    weights,
    event_cols,
    CI_lb_col,
    CI_ub_col,
    estimates,
    variances,
    events_trt,
    events_ctrl,
    n_trt,
    n_ctrl,
    rob_tool,
    outcome,
    model,
    ma = NULL,
    pubbias_min_studies = 10,
    null_effect = 0,
    industry = FALSE,
    search = FALSE,
    RoB_1_threshold = 1.5,
    RoB_2_threshold = 2.5,
    events_1_threshold = 300,
    events_2_threshold = 100,
    CI_threshold_pos = NULL,
    CI_threshold_neg = NULL,
    Jaccard_threshold = 0.4,
    variation_threshold = 0.8,
    Eggers_threshold = 0.1,
    indirectness = 0
) {
  
  #--------------------#
  # Predict RoB levels #
  #--------------------#
  
  # Calculate RoB_avg
  RoB_avg <- weighted_RoB(
    data = data, 
    rob_tool = rob_tool,
    overall_rob = overall_rob,
    random_selection = random_selection,
    allocation_selection = allocation_selection,
    performance = performance,
    detection = detection,
    attrition = attrition,
    reporting = reporting,
    other = other,
    weights = weights
  )
  
  # Set levels
  RoB_levels <- RoB_downgrades(
    RoB_avg = RoB_avg,
    one_level_threshold = RoB_1_threshold,
    two_level_threshold = RoB_2_threshold
  )
  
  #----------------------------#
  # Predict Imprecision levels #
  #----------------------------#
  
  # Calculate total number of events
  if (outcome == "RR") {
    total_events <- sum_events(
      data = data,
      event_cols = event_cols
    )
  } else {
    total_events <- NULL
  }
  
  # Conduct meta-analysis (if not already given)
  if (is.null(ma)) {
    if (model == "MH") {
      meta <- metafor::rma.mh(ai = data[[events_trt]], ci = data[[events_ctrl]], n1i = data[[n_trt]], n2i = data[[n_ctrl]], measure = outcome,
                              drop00 = c(TRUE, TRUE), add = c(0.5, 0.5), to = c("only0", "only0"))
    } else {
      meta <- metafor::rma(yi = data[[estimates]], vi = data[[variances]], measure = outcome, method = model)
    }
  } else {
    meta <- ma
  }
  
  # Set levels
  imprecision_levels <- imprecision_downgrades(
    outcome = outcome,
    total_events = total_events,
    event_threshold_2 = events_2_threshold,
    event_threshold_1 = events_1_threshold,
    CI_lb = meta$ci.lb,
    CI_ub = meta$ci.ub,
    CI_threshold_pos = CI_threshold_pos,
    CI_threshold_neg = CI_threshold_neg
  )
  
  #------------------------------#
  # Predict Inconsistency levels #
  #------------------------------#
  
  # Set levels
  inconsistency_levels <- inconsistency_downgrades(
    data = data,
    CI_lb_col = CI_lb_col,
    CI_ub_col = CI_ub_col,
    estimates = estimates,
    null_effect = null_effect,
    Jaccard_threshold = Jaccard_threshold,
    variation_threshold = variation_threshold
  )
  
  #---------------------------------#
  # Predict Publication bias levels #
  #---------------------------------#
  
  pubbias_levels <- pubbias_downgrades(
    data = data,
    estimates = estimates,
    variances = variances,
    events_trt = events_trt,
    events_ctrl = events_ctrl,
    n_trt = n_trt,
    n_ctrl = n_ctrl,
    model = model,
    min_studies = pubbias_min_studies,
    threshold = Eggers_threshold,
    industry = industry,
    search = search
  )
  
  #----------------------#
  # Predict GRADE rating #
  #----------------------#
  
  rating <- max(1, 
                4 - RoB_levels - imprecision_levels - inconsistency_levels - pubbias_levels - indirectness)
  rating <- factor(rating,
                   levels = c(1,2,3,4),
                   labels = c("Very low", "Low", "Moderate", "High"))
  
  # Return #
  return(list(
    result = as.character(rating),
    RoB = RoB_levels,
    Imprecision = imprecision_levels,
    Inconsistency = inconsistency_levels,
    Pubbias = pubbias_levels,
    Indirectness = indirectness
  ))
  
}

