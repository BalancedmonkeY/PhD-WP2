#--------------------------------------------------------------#
# Tools for estimating any downgrading due to publication bias #
#--------------------------------------------------------------#

library(metafor)

#' @param data Dataset containing information (including variance/standard error, and on log scale for ratios) as calculated by metafor escalc() command
#' @param estimates Column name for study estimates
#' @param variances Column name for study variances
#' @param threshold Threshold value of funnel plot stat that means downgrading by one level
#' @param industry Logical argument for whether there exists industry influence
#' @param search Logical argument for whether there are concerns regarding the search integrity
#' @return levels = Number of levels the evidence is likely to be downgraded due to publication bias (0 or 1)

pubbias_downgrades <- function(
    data,
    estimates,
    variances,
    threshold = 0.9,
    industry = FALSE,
    search = FALSE,
) {
  
  # Check that there are at least 5 studies
  n <- nrow(data)
  
  if (n < 5) {
    print("The funnel plot statistic is not suitable for datasets with less than 5 studies")
    stat <- NA
  } else {
    
    # Calculate funnel plot statistic (Egger's)
    res <- metafor::rma(data = data, yi = estimates, vi = variances)
    test_info <- metafor::regtest(res, model = "lm", predictor = "sei")
    stat <- test_info$pval
  }
  
  # Downgrade
  if (stat >= threshold | industry | search) {
    levels = 1
  } else {
    levels = 0
  }
  
  return(levels = levels)
  
}