#--------------------------------------------------------------#
# Tools for estimating any downgrading due to publication bias #
#--------------------------------------------------------------#

library(metafor)

#' @param data Dataset containing information (including variance/standard error, and on log scale for ratios) as calculated by metafor escalc() command
#' @param estimates Column name for study estimates
#' @param variances Column name for study variances
#' @param events_trt Column name that refers to the number of events in the treatment arm
#' @param events_ctrl Column name that refers to the number of events in the control arm
#' @param n_trt Column name that refers to the total number of people in the treatment arm
#' @param n_ctrl Column name that refers to the total number of people in the control arm
#' @param model Meta-analysis model (as per metafor options (EE, MH, or DL))
#' @param min_studies Minimum number of studies needed to calculate the funnel asymetry statistic
#' @param threshold Threshold value of funnel plot stat that means downgrading by one level
#' @param industry Logical argument for whether there exists industry influence
#' @param search Logical argument for whether there are concerns regarding the search integrity
#' @return levels = Number of levels the evidence is likely to be downgraded due to publication bias (0 or 1)

pubbias_downgrades <- function(
    data,
    estimates,
    variances,
    events_trt,
    events_ctrl,
    n_trt,
    n_ctrl,
    model,
    min_studies = 10,
    threshold = 0.1,
    industry = FALSE,
    search = FALSE
) {
  
  # Check that there are at least 10 studies
  n <- nrow(data)
  
  if (n < min_studies) {
    paste0("The funnel plot statistic is not suitable for datasets with less than ", min_studies, " studies")
    stat <- 1 # this ensures that no downgrading will happen due to asymmetry
  } else {
    
    # Calculate funnel plot statistic (Egger's)
    if (model == "MH") {
      res <- metafor::rma.mh(ai = data[[events_trt]], ci = data[[events_ctrl]], n1i = data[[n_trt]], n2i = data[[n_ctrl]],
                             drop00 = c(TRUE, TRUE), add = c(0.5, 0.5), to = c("only0", "only0"))
    } else {
      res <- metafor::rma(yi = data[[estimates]], vi = data[[variances]])
    }
    test_info <- metafor::regtest(res, model = "lm", predictor = "sei")
    stat <- test_info$pval
  }
  
  # Downgrade
  if (stat < threshold | industry | search) {
    levels = 1
  } else {
    levels = 0
  }
  
  return(levels = levels)
  
}