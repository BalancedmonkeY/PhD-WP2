#-----------------------------------------------------------#
# Tools for estimating any downgrading due to inconsistency #
#-----------------------------------------------------------#

library(dplyr)

#' @param CI_lb_1 Lower bound of first CI
#' @param CI_ub_1 Upper bound of first CI
#' @param CI_lb_2 Lower bound of second CI
#' @param CI_ub_2 Upper bound of second CI
#' @return Jaccard = Jaccard index for the two CIs

Jaccard_index <- function(
    CI_lb_1,
    CI_ub_1,
    CI_lb_2,
    CI_ub_2
) {
  
  intersection <- max(0, min(CI_ub_1, CI_ub_2) - max(CI_lb_1, CI_lb_2))
  union <- max(CI_ub_1, CI_ub_2) - min(CI_lb_1, CI_lb_2)
  Jaccard <- intersection / union
  
  return(Jaccard)
  
}

#' @param CI_lb_col Column of data that refer to the CI lower bound (in log-form for ratios)
#' @param CI_ub_col Column of data that refers to the CI upper bound (in log forms for ratios)
#' @return index = average pairwise Jaccard Index of the data

pairwise_Jaccard <- function(
    CI_lb_col,
    CI_ub_col
) {
  
  # Transform for positive numbers
  smallest <- min(CI_lb_col)
  if (smallest < 0) {
    CI_lb_col <- CI_lb_col - smallest
    CI_ub_col <- CI_ub_col - smallest
  }
  
  # Calculate all pairwise Jaccard indexes and sum them
  n <- length(CI_lb_col)
  if (n > 1) {
    all_Jaccard <- sum(
      combn(n, 2, function(idx) {
        i <- idx[1]
        j <- idx[2]
        Jaccard_index(CI_lb_1 = CI_lb_col[i], CI_ub_1 = CI_ub_col[i], 
                      CI_lb_2 = CI_lb_col[j], CI_ub_2 = CI_ub_col[j])
      })
    )
    
    # Calculate average value
    index <- (2*all_Jaccard)/(n*(n-1))
    
  } else {
    
    print("There needs to be at least two rows of data")
    index <- NA
    
  }
  
  return(index)
  
}

#' @param threshold Value for which to assess the point estimates against
#' @param estimates Column of data that refers to the point estimates of all studies (in the same scale as the threshold)
#' @return variation = proportion of studies where the point estimate was one side of the threshold (between 0.5 and 1)

estimate_variation <- function(
    threshold = 0,
    estimates
) {
  
  variation <- sum(estimates > threshold) / length(estimates)
  if (variation < 0.5) {variation <- 1 - variation} # flip if counted on the 'other' side
  
  return(variation)
  
}

#' @param data Dataset containing information
#' @param CI_lb_col Column name that refer to the CI lower bound (in log-form for ratios)
#' @param CI_ub_col Column name that refers to the CI upper bound (in log forms for ratios)
#' @param estimates Column name that refers to the point estimates of all studies (in log forms for ratios)
#' @param null_effect Value for which to assess the point estimates against
#' @param Jaccard_threshold Threshold for the Jaccard index to suggest inconsistency
#' @param variation_threshold Threshold for the variation in estimates to suggest inconsistency
#' @param auto_adjust Logical parameter to indicate whether to automatically adjust the thresholds based on previous versions of the review
#' @param prev_levels The number of levels the evidence was downgraded due to inconsistency in the previous version
#' @param prev_Jaccard Jaccard index for previous update
#' @param prev_est_var Variation in estimates score for previous update
#' @return list containing: levels = Number of levels the evidence is likely to be downgraded due to inconsistency (0 or 1)
#'                          Jaccard_threshold = threshold value used (may differ if autoadjust was used)
#'                          variation_threshold = threshold value used (may differ if autoadjust was used)

inconsistency_downgrades <- function(
    data,
    CI_lb_col,
    CI_ub_col,
    estimates,
    null_effect = 0,
    Jaccard_threshold = 0.4,
    variation_threshold = 0.8,
    autoadjust = FALSE,
    prev_levels,
    prev_Jaccard,
    prev_est_var
) {
  
  # Autoadjustment if required
  if (auto_adjust) {
    if (prev_levels == 1) {
      if (prev_est_var >= variation_threshold) {
        variation_threshold <- ceiling((prev_est_var+0.01)*20)/20
        if (variation_threshold >= 1) { # may happen if the previous est_var equaled 1 exactly
          variation_threshold = 1
        }
        if (prev_est_var >= 0.95) { # i.e., threshold will become 1 which will never lead to a downgrading
          if (prev_Jaccard >= Jaccard_threshold) {
            Jaccard_threshold <- ceiling((prev_Jaccard+0.01)*20)/20
          }
        }
      } else {
        if (prev_Jaccard >= Jaccard_threshold) {
          Jaccard_threshold <- ceiling((prev_Jaccard+0.01)*20)/20
        }
      }
    } else if (prev_levels == 0) {
      if (prev_est_var < variation_threshold) {
        variation_threshold <- floor((prev_est_var)*20)/20
        if (prev_est_var < 0.55) { # i.e., threshold will become 0.5 which is is the lowest possible and indicative of checking the Jaccard
          if (prev_Jaccard < Jaccard_threshold) {
            Jaccard_threshold <- floor((prev_Jaccard)*20)/20
          }
        } 
      }
    }
  }
  
  # Obtain Jaccard index and variation in estimates
  Jaccard <- pairwise_Jaccard(
    CI_lb_col = data[[CI_lb_col]],
    CI_ub_col = data[[CI_ub_col]]
  )
  var_est <- estimate_variation(
    estimates = data[[estimates]]
  )
  
  # Assess whether to downgrade
  if (var_est >= variation_threshold) {
    levels = 0
  } else if (Jaccard >= Jaccard_threshold) {
    levels = 0
  } else {
    levels = 1
  }
  
  return(list(levels = levels,
              Jaccard_threshold = Jaccard_threshold,
              variation_threshold = variation_threshold))
  
}