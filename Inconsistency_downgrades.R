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
#' @return levels = Number of levels the evidence is likely to be downgraded due to inconsistency (0 or 1)

inconsistency_downgrades <- function(
    data,
    CI_lb_col,
    CI_ub_col,
    estimates,
    null_effect = 0,
    Jaccard_threshold = 0.4,
    variation_threshold = 0.8
) {
  
  # Assign zero if only one study
  data <- data %>% filter(!is.na(.data[[CI_lb_col]])) # remove rows where there may be NAs (often the case if zero events)
  if (nrow(data) == 1) {
    levels = 0
  } else {
  
    # Obtain Jaccard index and variation in estimates
    Jaccard <- pairwise_Jaccard(
      CI_lb_col = data[[CI_lb_col]],
      CI_ub_col = data[[CI_ub_col]]
    )
    var_est <- estimate_variation(
      estimates = data[[estimates]],
      threshold = null_effect
    )
  
    # Assess whether to downgrade
    if (var_est >= variation_threshold) {
      levels = 0
    } else if (Jaccard >= Jaccard_threshold) {
      levels = 0
    } else {
      levels = 1
    }
  }
  
  return(levels = levels)
  
}
