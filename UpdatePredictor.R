#------------------------------------------------------------------------#
# Tool to ascertain if new studies for a meta-analysis suggest an update #
#------------------------------------------------------------------------#

#-------------------------------#
# Clareece Nevill December 2025 #
#-------------------------------#

source("CertaintyTool.R")
source("InterpretationTool.R")

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
#' @param CI_lb_col Column name that refer to the CI lower bound (in log-form for ratios)
#' @param CI_ub_col Column name that refers to the CI upper bound (in log forms for ratios)
#' @param estimates Column name that refers to the point estimates of all studies (in log forms for ratios)
#' @param variances Column name that refers to the variances of each study
#' @param rob_tool What type of risk of bias tool was used (1 or 2)
#' @param outcome Outcome measure
#' @param model Meta-analysis model (as per metafor options)
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
    estimates,
    variances,
    rob_tool,
    outcome,
    model,
    null_effect = 0,
    industry = FALSE,
    search = FALSE,
    RoB_1_threshold = 1.5,
    RoB_2_threshold = NULL,
    events_1_threshold = 300,
    events_2_threshold = 100,
    CI_threshold_pos = NULL,
    CI_threshold_neg = NULL,
    Jaccard_threshold = 0.4,
    variation_threshold = 0.8,
    Eggers_threshold = 0.1,
    indirectness = 0,
    prev_RoB,
    prev_imprecision,
    prev_inconsistency,
    prev_pubbias,
    prev_indirectness,
    sig_level = 0.05,
    ylim = NULL,
    xlim = NULL,
    contour_points = 200,
    draw_plot = TRUE,
    legend = TRUE,
    expxticks = NULL,
    xticks = NULL,
    yticks = NULL,
    effect_zero = NA,
    effect_lower = NA,
    effect_upper = NA,
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
  
  prev_data <- data %>% filter(.data[[search_col]] <= last_search_date)
  new_data <- data %>% filter(.data[[search_col]] > last_search_date)
  
  #-----------------------------------#
  # Predict change in pooled estimate #
  #-----------------------------------#
  
  pooled_results <- InterpretationThreshold(
    SS = prev_data[[estimates]],
    seSS = sqrt(prev_data[[variances]]),
    SSnew = new_data[[estimates]],
    seSSnew = sqrt(new_data[[variances]]),
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
  
  #---------------------------------------------------#
  # Predict new GRADE rating if including all studies #
  #---------------------------------------------------#
  
  data$weights <- weights(pooled_results$new_ma)
  
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
    estimates = estimates,
    variances = variances,
    rob_tool = rob_tool,
    outcome = outcome,
    model = model,
    ma = pooled_results$new_ma,
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
  
  # Combine domains #
  GRADE_domains_text <- paste0(RoB_text, imprecision_text, inconsistency_text, pubbias_text, indirectness_text)
  if (GRADE_domains_text != "") {
    Certainty_text <- paste0("LSRUpdateR predicts that the addition of new studies will change the GRADE rating of evidence. Currently, the evidence is graded at ",
                             as.character(factor(max(4-prev_RoB-prev_imprecision-prev_inconsistency-prev_pubbias-prev_indirectness,1), levels = c(1,2,3,4), labels = c("Very low", "Low", "Moderate", "High"))),
    ", but LSRUpdateR predicts it will become ", newGRADE$result, " after the inclusion of new studies. ",
    GRADE_domains_text)
  } else {
    Certainty_text <- paste0("Currently, the evidence is graded at ",
                             as.character(factor(max(4-prev_RoB-prev_imprecision-prev_inconsistency-prev_pubbias-prev_indirectness,1), levels = c(1,2,3,4), labels = c("Very low", "Low", "Moderate", "High"))),
                             ". LSRUpdateR doesn't predict that this will change with the inclusion of new studies.")
  }
  
  # Any changes in pooled effect #
  if (grepl("will give", pooled_results$threshold_result)) {
    Interpretation_text <- paste0("LSRUpdateR predicts that the addition of new studies will change the interpretation of the results. ",
                                  pooled_results$threshold_result, ".")
  } else {
    Interpretation_text <- paste0("LSRUpdateR doesn't predict any change in interpretation: ", pooled_results$threshold_result, ".")
  }
  
  # Any changes for either #
  if (grepl("doesn't predict", Interpretation_text) & grepl("doesn't predict", Certainty_text)) {
    update_text <- paste("LSRUpdateR doesn't predict the need to conduct an update.",
                          Interpretation_text, Certainty_text)
  } else {
    update_text <- paste("LSRUpdateR predicts that the review should be updated for the following reasons.",
                  Interpretation_text, Certainty_text)
  }
  
  return(update_text)
  
  
}