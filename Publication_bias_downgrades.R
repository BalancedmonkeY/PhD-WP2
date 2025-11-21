#--------------------------------------------------------------#
# Tools for estimating any downgrading due to publication bias #
#--------------------------------------------------------------#

library(metafor)

#' @param data Dataset containing information (including variance/standard error, and on log scale for ratios) as calculated by metafor escalc() command
#' @param threshold Threshold value of funnel plot stat that means downgrading by one level
#' @param descrepancies Logical argument for whether there exist discrepancies between published and unpublished data
#' @param industry Logical argument for whether there exists industry influence
#' @param search Logical argument for whether there are concerns regarding the search integrity
#' @param auto_adjust Logical parameter to indicate whether to automatically adjust the thresholds based on previous versions of the review
#' @param prev_stat Funnel plot statistics from previous update
#' @param prev_levels The number of levels the evidence was downgraded due to publication bias in the previous version
#' @return list containing: levels = Number of levels the evidence is likely to be downgraded due to publication bias (0 or 1)
#'                          threshold = threshold value used (may differ if autoadjust was used)

pubbias_downgrades <- function(
    data,
    threshold = 0.9,
    discrepancies = FALSE,
    industry = FALSE,
    search = FALSE,
    auto_adjust = FALSE,
    prev_stat,
    prev_levels
) {
  
  # Check that there are at least 5 studies
  n <- nrow(data)
  
  if (n < 5) {
    print("The funnel plot statistic is not suitable for datasets with less than 5 studies")
    stat <- NA
  } else {
    
    # Auto adjustment if required
    if (autoadjust) {
      if (prev_levels == 0 & prev_stat >= threshold) {
        threshold <- ceiling((prev_stat+0.01)*20)/20
      } else if (prev_levels == 1 & prev_stat < threshold) {
        threshold <- floor((prev_stat)*20)/20
      }
    }
    
    # Calculate funnel plot statistic (Egger's)
    res <- metafor::rma(data = data, yi = yi, vi = vi)
    test_info <- metafor::regtest(res, model = "lm", predictor = "sei")
    stat <- test_info$pval
  }
  
  # Downgrade
  if (stat >= threshold | discrepancies | industry | search) {
    levels = 1
  } else {
    levels = 0
  }
  
  return(list(levels = levels,
              threshold = threshold))
  
}