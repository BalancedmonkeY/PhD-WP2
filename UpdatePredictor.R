#------------------------------------------------------------------------#
# Tool to ascertain if new studies for a meta-analysis suggest an update #
#------------------------------------------------------------------------#

#-------------------------------#
# Clareece Nevill December 2025 #
#-------------------------------#

source("../3. Create tool/CertaintyTool.R")
source("../3. Create tool/InterpretationTool.R")
source("../3. Create tool/ThresholdDescribers.R")

#' @param data Data frame where each row contains a study, and the columns include: (i) RoB ratings, (ii) estimates & CI bounds & vi (in log form if ratio), (iii) Weight, (iv) Number of events (if dichotomous), (v) Search data study was found from
#' @param search_col Column name indicating the search dates for which the respective study was identified
#' @param last_search_date Date that indicates the latest search date that was included in the most recent published update of the review (in year-month-date format)
#' @param overall_rob Overall rating of risk of bias
#' @param random_selection Random sequence generation (selection bias) rating (RoB 1 only)
#' @param allocation_selection Allocation concealment (selection bias) rating (RoB 1 only)
#' @param performance Performance bias rating (RoB 1 only)
#' @param detection Detection bias rating (RoB 1 only)
#' @param attrition Attrition bias rating (RoB 1 only)
#' @param reporting Reporting bias rating (RoB 1 only)
#' @param other Other bias rating (RoB 1 only)
#' @param event_cols Column names for events (if dichotomous)
#' @param CI_lb_col Column name that refer to the CI lower bound (keep in original units)
#' @param CI_ub_col Column name that refers to the CI upper bound (keep in original units)
#' @param estimates Column name that refers to the point estimates of all studies (in log forms for ratios)
#' @param variances Column name that refers to the variances of each study
#' @param events_trt_name Column name that refers to the number of events in the treatment arm
#' @param events_ctrl_name Column name that refers to the number of events in the control arm
#' @param n_trt_name Column name that refers to the total number of people in the treatment arm
#' @param n_ctrl_name Column name that refers to the total number of people in the control arm
#' @param rob_tool What type of risk of bias tool was used (1 or 2)
#' @param outcome Outcome measure
#' @param model Meta-analysis model (as per metafor options)
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
#' @param prev_RoB Previous levels downgraded for RoB
#' @param prev_imprecision Previous levels downgraded for imprecision
#' @param prev_inconsistency Previous levels downgraded for inconsistency
#' @param prev_pubbias Previous levels downgraded for publication bias
#' @param prev_indirectness Previous levels downgraded for indirectness
#' @param prev_est Previous pooled estimate (optional - suitable for when studies might be excluded from 'current_data')
#' @param prev_lb Lower CI bound for previous pooled estimate
#' @param prev_ub Upper CI bound for previous pooled estimate
#' @param prev_p_value P-value from previous meta-analysis
#' @param sig_level significance level
#' @param ylim limits of the y axis	in form c(y1, y2) 
#' @param xlim limits of the x axis in form c(x1, x2)
#' @param contour_points number of points for creating contours with a random-effects model - more means a smoother contour but takes longer to compute	
#' @param draw_plot TRUE/FALSE for drawing an extended funnel plot of the contour thresholds used
#' @param legend - TRUE/FALSE for displaying key/legend	
#' @param expxticks custom ticks for the x axis on a exponential scale (assumes data is already on log scale)
#' @param xticks custom ticks for the x axis
#' @param yticks custom ticks for the y axis
#' @param effect_zero value for the null effect (if interesting in statistical significant regarding the pooled effect)
#' @param effect_lower value for which the user wants the lower CI bound to be as high as (regarding pooled effect)
#' @param effect_upper value for which the user wants the upper CI bound to be as low as (regarding pooled effect)
#' @param effect_est_pos Value for which the user wants to specify the estimate to be as high as
#' @param effect_est_neg Value for which the user wants to specify the estimate to be as low as	
#' @param xlab label for the x axis		
#' @param ylab label for the y axis
#' @param plot_threshold TRUE/FALSE plot the threshold value vertical line as defined by the argument 'zero', 'lower', or 'upper'	
#' @param plot_summ_current TRUE/FALSE plot the updated pooled effect vertical line of current studies
#' @param plot_summ_updated TRUE/FALSE plot the updated pooled effect vertical line of current studies
#' @param legendpos position of legend as per ggplot styling	
#' @param summ_current TRUE/FALSE plot current summary diamond including pooled effect and confidence interval (significance level as defined by sig.level)
#' @param summ_updated TRUE/FALSE plot updated summary diamond including pooled effect and confidence interval (significance level as defined by sig.level)
#' @param summ_pos adjustment of position of summary diamond
#' @param new_points TRUE/FALSE add points of new study(ies) to plot
#' @param points - TRUE/FALSE whether the study points should be displayed at all (TRUE default)
#' @param pred_interval TRUE/FALSE display predictive interval along with the summary diamond
#' @param tau2 TRUE/FALSE display tau2 values (TRUE default)
#' @param rand.load TRUE/FALSE show percentage of computations complete when the random effects contours are calculated
#' returns a list of the following:
#' text_result = text description of whether or not to update
#' GRADE_results = output from certainty tool
#' pooled_results = output from interpretation tool

UpdatePredictor <- function(
    data,
    search_col,
    last_search_date,
    overall_rob,
    random_selection = NA,
    allocation_selection = NA,
    performance = NA,
    detection = NA,
    attrition = NA,
    reporting = NA,
    other = NA,
    event_cols,
    CI_lb_col,
    CI_ub_col,
    estimates = NULL,
    variances = NULL,
    events_trt_name = NULL,
    events_ctrl_name = NULL,
    n_trt_name = NULL,
    n_ctrl_name = NULL,
    rob_tool,
    outcome,
    model,
    null_effect = 0,
    pubbias_min_studies = 10,
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
    Eggers_threshold = 0.025,
    indirectness = 0,
    prev_RoB,
    prev_imprecision,
    prev_inconsistency,
    prev_pubbias,
    prev_indirectness,
    prev_est = NULL,
    prev_lb = NULL,
    prev_ub = NULL,
    prev_p_value = NULL,
    sig_level = 0.05,
    ylim = NULL,
    xlim = NULL,
    contour_points = 200,
    draw_plot = FALSE,
    legend = TRUE,
    expxticks = NULL,
    xticks = NULL,
    yticks = NULL,
    effect_zero = NA,
    effect_lower = NA,
    effect_upper = NA,
    effect_est_pos = NA,
    effect_est_neg = NA,
    xlab = paste0 ("Effect (", outcome, ")"),
    ylab = "Standard Error",
    plot_threshold = TRUE,
    plot_summ_current = TRUE,
    plot_summ_updated = TRUE,
    legendpos = NULL,
    summ_current = TRUE,
    summ_updated = TRUE,
    summ_pos = 0,
    new_points = TRUE,
    points = TRUE,
    pred_interval = FALSE,
    tau2 = TRUE,
    rand_load = 100
  )
{
  
  #--------------------------------------#
  # Split data into 'original' and 'new' #
  #--------------------------------------#
  
  prev_df <- data %>% filter(data[[search_col]] <= last_search_date)
  new_df <- data %>% filter(data[[search_col]] > last_search_date)
  
  #-----------------------------------#
  # Predict change in pooled estimate #
  #-----------------------------------#
  
  pooled_results <- InterpretationThreshold(
    SS <- if (!is.null(estimates)) prev_df[[estimates]] else NULL,
    seSS = if (!is.null(variances))  sqrt(prev_df[[variances]]) else NULL,
    SSnew = if (!is.null(estimates))  new_df[[estimates]] else NULL,
    seSSnew = if (!is.null(variances))  sqrt(new_df[[variances]]) else NULL,
    events_trt = if (!is.null(events_trt_name))  prev_df[[events_trt_name]] else NULL,
    events_ctrl = if (!is.null(events_ctrl_name))  prev_df[[events_ctrl_name]] else NULL,
    n_trt = if (!is.null(n_trt_name))  prev_df[[n_trt_name]] else NULL,
    n_ctrl = if (!is.null(n_ctrl_name))  prev_df[[n_ctrl_name]] else NULL,
    events_trt_new = if (!is.null(events_trt_name))  new_df[[events_trt_name]] else NULL,
    events_ctrl_new = if (!is.null(events_ctrl_name))  new_df[[events_ctrl_name]] else NULL,
    n_trt_new = if (!is.null(n_trt_name))  new_df[[n_trt_name]] else NULL,
    n_ctrl_new = if (!is.null(n_ctrl_name))  new_df[[n_ctrl_name]] else NULL,
    sig_level = sig_level,
    method = model,
    outcome = outcome,
    ylim = ylim,
    xlim = xlim,
    contour_points = contour_points,
    draw_plot = draw_plot,
    legend = legend,
    expxticks = expxticks,
    xticks = xticks,
    yticks = yticks,
    zero = effect_zero,
    lower = effect_lower,
    upper = effect_upper,
    est_pos = effect_est_pos,
    est_neg = effect_est_neg,
    xlab = xlab,
    ylab = ylab,
    plot_threshold = plot_threshold,
    plot_summ_current = plot_summ_current,
    plot_summ_updated = plot_summ_updated,
    legendpos = legendpos,
    summ_current = summ_current,
    summ_updated = summ_updated,
    summ_pos = summ_pos,
    new_points = new_points,
    points = points,
    pred_interval = pred_interval,
    tau2 = tau2,
    rand_load = rand_load
  )
  
  # Threshold result for the original meta-analysis (taken externally due to exclusion set-up)
  if (outcome %in% c('OR', 'RR')) {
    if (!is.na(effect_zero)) {
      effect_zero <- log(effect_zero)
    } else if (!is.na(effect_lower)) {
      effect_lower <- log(effect_lower)
    } else if (!is.na(effect_upper)) {
      effect_upper <- log(effect_upper)
    } else if (!is.na(effect_est_pos) | !is.na(effect_est_neg)) {
      if (!is.na(effect_est_pos)) {effect_est_pos <- log(effect_est_pos)}
      if (!is.na(effect_est_neg)) {effect_est_neg <- log(effect_est_neg)}
    }
    
    if (!is.null(prev_est)) {
      prev_est <- log(prev_est)
    }
    if (!is.null(prev_lb)) {
      prev_lb <- log(prev_lb)
    }
    if (!is.null(prev_ub)) {
      prev_ub <- log(prev_ub)
    }
    null_effect <- log(null_effect)
  }
  og_threshold_result <- Threshold_Description(est = prev_est, ci_lb = prev_lb, ci_ub = prev_ub, pvalue = prev_p_value, outcome = outcome,
                                               sig_level=sig_level, zero=effect_zero, lower=effect_lower, 
                                               upper=effect_upper, est_pos=effect_est_pos, est_neg=effect_est_neg,
                                               new_or_og = "og")
  
  #---------------------------------------------------#
  # Predict new GRADE rating if including all studies #
  #---------------------------------------------------#
  
  data$weights <- 0
  data$weights[!is.na(data$Mean)] <- weights(pooled_results$new_ma)
  
  
  newGRADE <- PredictedGRADEdomains(
    data = data,
    overall_rob = overall_rob,
    random_selection = random_selection,
    allocation_selection = allocation_selection,
    performance = performance,
    detection = detection,
    attrition = attrition,
    reporting = reporting,
    other = other,
    weights = "weights",
    event_cols = event_cols,
    CI_lb_col = CI_lb_col,
    CI_ub_col = CI_ub_col,
    estimates = if (!is.null(estimates)) estimates else NULL,
    variances = if (!is.null(variances)) variances else NULL,
    events_trt_name = if (!is.null(events_trt_name)) events_trt_name else NULL,
    events_ctrl_name = if (!is.null(events_ctrl_name)) events_ctrl_name else NULL,
    n_trt_name = if (!is.null(n_trt_name)) n_trt_name else NULL,
    n_ctrl_name = if (!is.null(n_ctrl_name)) n_ctrl_name else NULL,
    rob_tool = rob_tool,
    outcome = outcome,
    model = model,
    ma = pooled_results$new_ma,
    pubbias_min_studies = pubbias_min_studies,
    null_effect = null_effect,
    industry = industry,
    search = search,
    RoB_1_threshold = RoB_1_threshold,
    RoB_2_threshold = RoB_2_threshold,
    events_1_threshold = events_1_threshold,
    events_2_threshold = events_2_threshold,
    CI_threshold_pos = CI_threshold_pos,
    CI_threshold_neg = CI_threshold_neg,
    Jaccard_threshold = Jaccard_threshold,
    variation_threshold = variation_threshold,
    Eggers_threshold = Eggers_threshold,
    indirectness = indirectness
  ) 
  
  #---------------------------------#
  # Assess if there are any changes #
  #---------------------------------#
  
  
  # Any changes in domains #
  if (prev_RoB != newGRADE$RoB) {
    RoB_text <- paste0("Currently, the evidence has been downgraded ", prev_RoB, 
                       " levels due to risk of bias. LSRUpdateR predicts the addition of new studies will downgrade the evidence by ", 
                       newGRADE$RoB, " level(s) instead.")
  } else {RoB_text <- ""}
  
  if (prev_imprecision != newGRADE$Imprecision) {
    imprecision_text <- paste0("Currently, the evidence has been downgraded ", prev_imprecision, 
                               " levels due to imprecision. LSRUpdateR predicts the addition of new studies will downgrade the evidence by ", 
                               newGRADE$Imprecision, " level(s) instead.")
  } else {imprecision_text <- ""}
  
  if (prev_inconsistency != newGRADE$Inconsistency) {
    inconsistency_text <- paste0("Currently, the evidence has been downgraded ", prev_inconsistency, 
                               " levels due to inconsistency. LSRUpdateR predicts the addition of new studies will downgrade the evidence by ", 
                               newGRADE$Inconsistency, " level(s) instead.")
  } else {inconsistency_text <- ""}
  
  if (prev_pubbias != newGRADE$Pubbias) {
    pubbias_text <- paste0("Currently, the evidence has been downgraded ", prev_pubbias, 
                                 " levels due to publication bias. LSRUpdateR predicts the addition of new studies will downgrade the evidence by ", 
                                 newGRADE$Pubbias, " level(s) instead.")
  } else {pubbias_text <- ""}
  
  if (prev_indirectness != newGRADE$Indirectness) {
    indirectness_text <- paste0("Currently, the evidence has been downgraded ", prev_indirectness, 
                           " levels due to indirectness. LSRUpdateR predicts the addition of new studies will downgrade the evidence by ", 
                           newGRADE$Indirectness, " level(s) instead.")
  } else {indirectness_text <- ""}
  
  # Calculate whether change of overall rating & add respective text
  old_GRADE_rating <- as.character(factor(max(4-prev_RoB-prev_imprecision-prev_inconsistency-prev_pubbias-prev_indirectness,1), levels = c(1,2,3,4), labels = c("Very low", "Low", "Moderate", "High")))
  if (old_GRADE_rating != newGRADE$result) {
    GRADE_domains_text <- paste0(RoB_text, imprecision_text, inconsistency_text, pubbias_text, indirectness_text)
    Certainty_text <- paste0("LSRUpdateR predicts that the addition of new studies will change the GRADE rating of evidence. Currently, the evidence is graded at ",
                             old_GRADE_rating, ", but LSRUpdateR predicts it will become ", newGRADE$result, " after the inclusion of new studies. ",
    GRADE_domains_text)
  } else {
    Certainty_text <- paste0("Currently, the evidence is graded at ", old_GRADE_rating,
                             ". LSRUpdateR doesn't predict that this will change with the inclusion of new studies.")
  }
  
  # Any changes in pooled effect #
  # Indicators of whether current and updated meta-analysis hit any thresholds
  og_hit <- grepl("do give", og_threshold_result)
  new_hit <- grepl("will give", pooled_results$threshold_result)
  # see if changed
  if (og_hit != new_hit) {
    Interpretation_text <- paste0("LSRUpdateR predicts that the addition of new studies will change the interpretation of the results. ",
                                  og_threshold_result, ". ",
                                  pooled_results$threshold_result, ".")
  } else {
    Interpretation_text <- paste0("LSRUpdateR doesn't predict any change in interpretation: ", og_threshold_result, ". ", pooled_results$threshold_result, ".")
  }
  
  # Any changes for either #
  if (grepl("doesn't predict", Interpretation_text) & grepl("doesn't predict", Certainty_text)) {
    update_text <- paste("LSRUpdateR doesn't predict the need to conduct an update.",
                          Interpretation_text, Certainty_text)
  } else {
    update_text <- paste("LSRUpdateR predicts that the review should be updated for the following reasons.",
                  Interpretation_text, Certainty_text)
  }
  
  return(list(
    text_result = update_text,
    GRADE_results = newGRADE,
    pooled_results = pooled_results
    ))
  
  
}