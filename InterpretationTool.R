#----------------------------------------------------------------------------------------------------------#
# Tool to ascertain if new studies for a meta-analysis will change the interpretation of the pooled effect #
#----------------------------------------------------------------------------------------------------------#
#----------------------------------------#
# Developed using code from Langan et al #
#----------------------------------------#
#-----------------------------#
# Clareece Nevill August 2025 #
#-----------------------------#

library(ggplot2)
library(metafor)
library(progressr)
library(rmeta)

#' @param SS effect estimates of the current studies (log scale if ratio)		
#' @param seSS standard errors of the current studies (log scale if ratio)
#' @param SSnew effect estimates of the new study(ies)
#' @param seSSnew standard errors of the new study(ies)
#' @param events_trt Number of events in treatment arm of current studies
#' @param events_ctrl Number of events in control arm of current studies
#' @param n_trt Total sample size of treatment arm of current studies
#' @param n_ctrl Total sample size of control arm of current studies
#' @param events_trt_new Number of events in treatment arm of new studies
#' @param events_ctrl_new Number of events in control arm of new studies
#' @param n_trt_new Total sample size of treatment arm of new studies
#' @param n_ctrl_new Total sample size of control arm of new studies
#' @param sig_level significance level
#' @param method type of meta-analysis - "EE", "MH" or "DL"
#' @param outcome type of outcome measure - "OR", "RR", "RD", "MD", "SMD"
#' @param ylim limits of the y axis	in form c(y1, y2) 
#' @param xlim limits of the x axis in form c(x1, x2)
#' @param contour_points number of points for creating contours with a random-effects model - more means a smoother contour but takes longer to compute	
#' @param draw_plot TRUE/FALSE for drawing an extended funnel plot of the contour thresholds used
#' @param legend - TRUE/FALSE for displaying key/legend	
#' @param expxticks custom ticks for the x axis on a exponential scale (assumes data is already on log scale)
#' @param xticks custom ticks for the x axis
#' @param yticks custom ticks for the y axis
#' @param zero value for the null effect (usually 0, even when expxticks is used for odds ratios)
#' @param lower value for which the user wants the lower CI bound to be as high as
#' @param upper value for which the user wants the upper CI bound to be as low as	
#' @param est_pos Value for which the user wants to specify the estimate to be as high as
#' @param est_neg Value for which the user wants to specify the estimate to be as low as
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
#' @return List containing the following: 
#' threshold_result - a string containing a description of whether the threshold of the new analysis was met;
#' threshold_plot - the extended funnel plot that visualises the threshold and whether it was met
#' original_ma - rma object containing meta-analysis details of original dataset
#' new_ma - ma object containing meta-analysis details of new dataset              
InterpretationThreshold <- function(
    SS = NULL,
    seSS = NULL,
    SSnew = NULL,
    seSSnew = NULL,
    events_trt = NULL,
    events_ctrl = NULL,
    n_trt = NULL,
    n_ctrl = NULL,
    events_trt_new = NULL,
    events_ctrl_new = NULL,
    n_trt_new = NULL,
    n_ctrl_new = NULL,
    sig_level = 0.05,
    method = 'EE',
    outcome = 'RR',
    ylim = NULL,
    xlim = NULL,
    contour_points = 200,
    draw_plot = TRUE,
    legend = TRUE,
    expxticks = NULL,
    xticks = NULL,
    yticks = NULL,
    zero = NA,
    lower = NA,
    upper = NA,
    est_pos = NA,
    est_neg = NA,
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
) {
  

  #----------------------------------------#
  # Calculate initial arguments/parameters #
  #----------------------------------------#
  
  #Convert the significance level into 'z' from the normal distribution (i.e. 0.05 -> 1.96)
  ci <- qnorm(1 - sig_level / 2)
  
  #Calculate updated meta-analysis (using rma from {metafor})
  if (method == "MH") {
    updated_meta <- metafor::rma.mh(ai = c(events_trt, events_trt_new), ci = c(events_ctrl,events_ctrl_new), 
                                    n1i = c(n_trt, n_trt_new), n2i = c(n_ctrl, n_ctrl_new), level = (1 - sig_level), measure = outcome,
                                    drop00 = c(TRUE, TRUE), add = c(0.5, 0.5), to = c("only0", "only0"))
  } else {
   updated_meta <- metafor::rma(yi = c(SS, SSnew), sei = c(seSS, seSSnew), method = method, level = (1 - sig_level), measure = outcome) 
  }
  
  updated_tau2 <- updated_meta$tau2
  
  # Calculate SS and seSS when MH is chosen
  if (method == "MH") {
    SS <- as.vector(na.omit(updated_meta$yi.f[1:length(events_trt)]))
    SSnew <- as.vector(na.omit(updated_meta$yi.f[(length(events_trt) + 1):(length(events_trt) + length(events_trt_new))]))
    seSS <- sqrt(as.vector(na.omit(updated_meta$vi.f[1:length(events_trt)])))
    seSSnew <- sqrt(as.vector(na.omit(updated_meta$vi.f[(length(events_trt) + 1):(length(events_trt) + length(events_trt_new))])))
  }
  
  #Calculate current meta-analysis (using rma from {metafor}) (if there is suitable data)
  if (length(na.omit(SS)) != 0) {
    if (method == "MH") {
      current_meta <- metafor::rma.mh(ai = events_trt, ci = events_ctrl, n1i = n_trt, n2i = n_ctrl, level = (1 - sig_level), measure = outcome,
                                      drop00 = c(TRUE, TRUE), add = c(0.5, 0.5), to = c("only0", "only0"))
    } else {
      current_meta <- metafor::rma(yi = SS, sei = seSS, method = method, level = (1 - sig_level), measure = outcome)
    }
    current_tau2 <- current_meta$tau2
  } else {
    current_meta <- "The 'current' meta-analysis was not suitable as there was no estimable data"
  }
  
  
  #---------------------------#
  # NEW WEIGHTINGS of studies #
  #---------------------------#
  
  # Weights are calculated the same whether 'current' or 'new'
  if (method == "MH") { # weight estimate from Greenland and Robins
    # zero adjustment
    zero_row <- xor(events_ctrl == 0, events_trt == 0)
    size_current <- ifelse(zero_row,
                           ((events_ctrl + 0.5)*(n_trt + 1))/((n_trt + 1) + (n_ctrl + 1)),
                           (events_ctrl*n_trt)/(n_trt + n_ctrl))
    zero_row_new <- xor(events_ctrl_new == 0, events_trt_new == 0)
    size_new <- ifelse(zero_row_new,
                       ((events_ctrl_new + 0.5)*(n_trt_new + 1))/((n_trt_new + 1) + (n_ctrl_new + 1)),
                       (events_ctrl_new*n_trt_new)/(n_trt_new + n_ctrl_new))
  } else if (method == "DL") {
    size_current <- 1 / ((seSS^2) + updated_tau2)  # standard inverse-variance weighting
    size_new <- 1 / ((seSSnew^2) + updated_tau2)
  } else {
    size_current <- 1 / (seSS^2)
    size_new <- 1 / (seSSnew^2)
  }
  
  # New studies 'average' (only when multiple)
  new_avg_est <- sum(size_new * SSnew) / sum(size_new)
  new_avg_se <- sqrt(1 / sum(size_new))

  
  #-----------------------------------#
  # Intelligent y-axis default limits #
  #-----------------------------------#
  
  if (length(na.omit(SS)) != 0) {
  
  sediff <- max(c(seSS, seSSnew), na.rm = TRUE) - min(c(seSS, seSSnew), na.rm = TRUE)
  
  if (!is.null(ylim) && ylim[1] < ylim[2]) {
    ylim <- rev(ylim)  # for when the user has already defined y limits
  }
  
  if (is.null(ylim)) {
    ylim <- c(max(c(seSS, seSSnew), na.rm = TRUE) + 0.20 * sediff, min(c(seSS, seSSnew), na.rm = TRUE) - 0.25 * sediff)
    if (ylim[2] < 0) {
      ylim[2] <- 0
    }
  }
  
  axisdiff <- ylim[2] - ylim[1]
  
  }
  
  #-----------------------------------#
  # Intelligent x-axis default limits #
  #-----------------------------------#
  
  if (length(na.omit(SS)) != 0) {
  
  # log if user specified and outcome is a ratio
  if (!is.null(xlim) & outcome %in% c('OR', 'RR')) {
    xlim <- log(xlim)
  }
  
  # Default limits
  SSdiff <- max(c(SS, SSnew), na.rm = TRUE) - min(c(SS, SSnew), na.rm = TRUE)
  
  if (is.null(xlim)) {
    xlim <- c(min(c(SS, SSnew), na.rm = TRUE) - 0.2 * SSdiff, max(c(SS, SSnew), na.rm = TRUE) + 0.2 * SSdiff)
  }
  
  }

  
  #-------------------------------------------------#
  # Vector of weights/sizes to define contours upon #
  #-------------------------------------------------#
  
  if (length(na.omit(SS)) != 0) {
  
   cSS <- seq(xlim[1], xlim[2], length.out = contour_points)  # granulated vector for effect size (x-axis)
   csize <- seq(ylim[1], ylim[2], length.out = contour_points)  # granulated vector for standard error (y-axis)
   csize[csize <= 0] <- 0.0000001 * min(c(seSS, seSSnew))
   for (k in 2:length(csize)) {
     if (csize[k] == 0 & csize[k-1] == 0) {
       csize[k] <- NA
     }
   }
   csize <- csize[!is.na(csize)]   # remove unnecessary data points
   
  }
  
  #------------------------------#
  # Create significance contours #
  #------------------------------#
   
  # Transform user-specified thresholds if outcome is a ratio
  if (outcome %in% c('OR', 'RR')) {
   if (!is.na(zero)) {
     zero <- log(zero)
   } else if (!is.na(lower)) {
      lower <- log(lower)
   } else if (!is.na(upper)) {
      upper <- log(upper)
   } else if (!is.na(est_pos) | !is.na(est_neg)) {
      if (!is.na(est_pos)) {est_pos <- log(est_pos)}
      if (!is.na(est_neg)) {est_neg <- log(est_neg)}
   }
  }
  
  if (length(na.omit(SS)) != 0) {
  
  # fixed-effect model #
  if (method %in% c("EE", "MH"))  {
    vwt <- 1 / (csize^2)     # weight for each point on plot (inverse of variance) (i.e. weight of new study)
    if (!is.na(zero)) {
      c1SS <- (1 / vwt) * (zero * (sum(size_current) + vwt) - sum(size_current * SS) +  ci * (sum(size_current) + vwt)^0.5)   # formulae based on CI boundaries of new MA meeting no effect
      c2SS <- (1 / vwt) * (zero * (sum(size_current) + vwt) - sum(size_current * SS) -  ci * (sum(size_current) + vwt)^0.5)
    } else if (!is.na(lower)) {
      c1SS <- (1 / vwt) * (lower * (sum(size_current) + vwt) - sum(size_current * SS) +  ci * (sum(size_current) + vwt)^0.5)
    } else if (!is.na(upper)) {
      c2SS <- (1 / vwt) * (upper * (sum(size_current) + vwt) - sum(size_current * SS) -  ci * (sum(size_current) + vwt)^0.5)
    } else if (!is.na(est_pos) | !is.na(est_neg)) {
      if (!is.na(est_pos)) {
        c1SS <- (1 / vwt) * (est_pos * (sum(size_current) + vwt) - sum(size_current * SS))
      }
      if (!is.na(est_neg)) {
        c2SS <- (1 / vwt) * (est_neg * (sum(size_current) + vwt) - sum(size_current * SS))
      }
    }
  }
    
  # random-effects model #
  if (method == "DL")  {
    
    # calculate tau2 adjustment factor K
    # Create V, X, Y, and Z
    size_current_fixed <- 1 / (seSS^2)
    V <- sum(size_current_fixed^2)
    X <- sum(size_current_fixed * SS)
    Y <- sum(size_current_fixed)
    Z <- sum(size_current_fixed * SS^2)
    K <- updated_tau2 / ((Z + (new_avg_est^2/new_avg_se^2) - ((X + new_avg_est/new_avg_se^2)^2/(Y + 1/new_avg_se^2)) - (length(SS) + length(SSnew) - 1))/(Y + 1/new_avg_se^2 - ((V + 1/new_avg_se^4)/(Y + 1/new_avg_se^2))))
    
    # optimised grid of all contour points
    contour_tiles <- expand.grid(i = seq_along(cSS), j = seq_along(csize))
    
    progressr::with_progress({
      p <- progressr::progressor(steps = nrow(contour_tiles))  # progress bar
    
      # calculate results for each grid point
      results <- purrr::pmap(contour_tiles, function(i, j) {
        
        if (length(SSnew) == 1) {
          # using rmeta::meta.summaries as it's much faster than metafor::rma
          metacont <- rmeta::meta.summaries(d = c(SS, cSS[i]), se = c(seSS, csize[j]), method = "random", conf.level = (1-sig_level))
          est <- metacont$summary
          lc <- metacont$summary - ci*metacont$se.summary
          uc <- metacont$summary + ci*metacont$se.summary
        # When there are multiple new studies, we are estimating tau2 using method of moments and an adjustment factor  
        } else if (length(SSnew) > 1) {
          tau2_est <- ((Z + (cSS[i]^2/csize[j]^2) - ((X + cSS[i]/csize[j]^2)^2/(Y + 1/csize[j]^2)) - (length(SS) + length(SSnew) - 1))/(Y + 1/csize[j]^2 - ((V + 1/csize[j]^4)/(Y + 1/csize[j]^2)))) * K
          if (tau2_est < 0) {tau2_est <- 0} # ensures no errors
          size_current_pixel <- 1 / ((seSS^2) + tau2_est)
          sc_sum <- sum(size_current_pixel)
          est <- (sum(size_current_pixel * SS) + cSS[i]/(csize[j]^2)) / (sc_sum + 1/(csize[j]^2))
          se <- sqrt(1/(sc_sum + 1/(csize[j]^2)))
          lc <- est - ci*se
          uc <- est + ci*se
        }
        
        p() #update progress
        
        # code according to threshold
        if (!is.na(zero)) {
          if (lc < zero & uc < zero) return("sigless_col")   # sig < 0
          if (lc < zero & uc > zero) return("nosig_col")   # not sig
          if (lc > zero & uc > zero) return("sigmore_col")   # sig > 0
        } else if (!is.na(lower)) {
          if (lc > lower) return("clinsig_col")
          if (lc <= lower) return("noclinsig_col")
        } else if (!is.na(upper)) {
          if (uc < upper) return("clinsig_col")
          if (uc >= upper) return("noclinsig_col")
        } else if (!is.na(est_pos) & !is.na(est_neg)) {
          if (est < est_neg) return("sigless_col")  
          if (est >= est_neg & est <= est_pos) return("nosig_col")
          if (est > est_pos) return("sigmore_col")
        } else if (!is.na(est_pos)) {
          if (est <= est_pos) return("nosig_col")
          if (est > est_pos) return("sigmore_col")
        } else if (!is.na(est_neg)) {
          if (est < est_neg) return("sigless_col")  
          if (est >= est_neg) return("nosig_col")
        }
        NA_character_
      })
      
      # Assign values to contour_tiles
      contour_tiles$code <- unlist(results)
      
    })
    
    contour_tiles <- cbind(contour_tiles, expand.grid(cSS = cSS, csize = csize))
  }
    
  }
   
  #------------------------------#
  # Assess whether threshold met #
  #------------------------------#  
   
   source("../3. Create tool/ThresholdDescribers.R")
   
   # Threshold result for new studies
   threshold_result <- Threshold_Description(meta = updated_meta, sig_level=sig_level, zero=zero, lower=lower, outcome = outcome, 
                                             upper=upper, est_pos=est_pos, est_neg=est_neg, SSnew=SSnew,
                                             new_or_og = "new")
  
  
  
  #------------------#
  # Summary diamonds #
  #------------------#
   
   if (length(na.omit(SS)) != 0) {
  
  if (summ_current & length(na.omit(SS)) != 0) {
    summary_diamond_current <- data.frame(
      xsumm = c(current_meta$ci.lb, current_meta$b, current_meta$ci.ub, current_meta$b),
      ysumm = c(ylim[2] - 0.10 * axisdiff + summ_pos, ylim[2] - 0.07 * axisdiff + summ_pos, 
                ylim[2] - 0.10 * axisdiff + summ_pos, ylim[2] - 0.13 * axisdiff + summ_pos)
    )
    
    if (pred_interval & length(na.omit(SS)) != 0) {	
      if (method == 'DL') {
        predint1_current <- predict(current_meta)$pi.lb
        predint2_current <- predict(current_meta)$pi.ub
        # update x-axis limits if predictive interval is wider
        xlim <- c(min(predint1_current - 0.2 * SSdiff, xlim[1]), max(predint2_current + 0.2 * SSdiff, xlim[2]))
      } else {
        print("For fixed-effects models, tau-squared is equal to 0 and therefore the PI becomes equivalent to the CI")
      }
    }
  }
   
   if (summ_updated) {
     summary_diamond_updated <- data.frame(
       xsumm = c(updated_meta$ci.lb, updated_meta$b, updated_meta$ci.ub, updated_meta$b),
       ysumm = c(ylim[2] - 0.10 * axisdiff + summ_pos, ylim[2] - 0.07 * axisdiff + summ_pos, 
                 ylim[2] - 0.10 * axisdiff + summ_pos, ylim[2] - 0.13 * axisdiff + summ_pos)
     )
     
     if (pred_interval) {	
       if (method == 'DL') {
         predint1_updated <- predict(updated_meta)$pi.lb
         predint2_updated <- predict(updated_meta)$pi.ub
         # update x-axis limits if predictive interval is wider
         xlim <- c(min(predint1_updated - 0.2 * SSdiff, xlim[1]), max(predint2_updated + 0.2 * SSdiff, xlim[2]))
       } else {
         print("For fixed-effects models, tau-squared is equal to 0 and therefore the PI becomes equivalent to the CI")
       }
     }
   } 
  
  #---------------#
  # Design legend #
  #---------------#
  
  # empty data frame ready for filling (one for each type of legend)
  legendmat.col.values <- NULL
  legendmat.col <- data.frame(labels = rep(NA, 7), linetype = rep(NA, 7), shape = rep(NA, 7), color = rep(NA, 7), fill = rep(NA, 7))
  legendmat.fill.values <- NULL
  legendmat.fill.labels <- NULL
  
  if (points) {
    legendmat.col.values <- c(legendmat.col.values, "point_col" = "black")
    legendmat.col$labels[1] <- "Current studies"
    legendmat.col$linetype[1] <- "blank"
    legendmat.col$shape[1] <- 19
    legendmat.col$color[1] <- "black"
  }
  
  if (new_points) {
    legendmat.col.values <- c(legendmat.col.values, "new_point_col" = "black")
    legendmat.col$labels[2] <- ifelse(length(SSnew) > 1, "New studies", "New study")
    legendmat.col$linetype[2] <- "blank"
    legendmat.col$shape[2] <- 21
    legendmat.col$color[2] <- "black"
    legendmat.col$fill[2] <- "lightblue"
  }
  
  if (new_points & length(SSnew) > 1) {
    legendmat.col.values <- c(legendmat.col.values, "new_point_avg_col" = "black")
    legendmat.col$labels[3] <- "New studies 'average'"
    legendmat.col$linetype[3] <- "blank"
    legendmat.col$shape[3] <- 24
    legendmat.col$color[3] <- "black"
    legendmat.col$fill[3] <- "lightblue"
  }
  
  if (pred_interval) {
    legendmat.col.values <- c(legendmat.col.values, "pred_col" = "black")
    legendmat.col$labels[4] <- "95% Predictive Interval"
    legendmat.col$linetype[4] <- "solid"
    legendmat.col$color[4] <- "black"
  }
  
  if (plot_summ_current) {
    legendmat.col.values <- c(legendmat.col.values, "summ_current_col" = "slategrey")
    legendmat.col$labels[5] <- "Current Pooled Effect"
    legendmat.col$linetype[5] <- "solid"
    legendmat.col$color[5] <- "slategrey"
  }
  
  if (plot_summ_updated) {
    legendmat.col.values <- c(legendmat.col.values, "summ_updated_col" = "cadetblue4")
    legendmat.col$labels[6] <- "Updated Pooled Effect"
    legendmat.col$linetype[6] <- "solid"
    legendmat.col$color[6] <- "cadetblue4"
  }
  
  if (plot_threshold) {
    legendmat.col.values <- c(legendmat.col.values, "threshold_col" = "lightgrey")
    legendmat.col$labels[7] <- "Threshold value"
    legendmat.col$linetype[7] <- "solid"
    legendmat.col$color[7] <- "lightgray"
  }
  
  if (summ_current) {
    legendmat.fill.values <- c(legendmat.fill.values, "diamond_fill_current" = "lavenderblush4")
    legendmat.fill.labels <- c(legendmat.fill.labels, paste0("Current Pooled Result (diamond - ", round((1-sig_level)*100, 1), "% CI)"))
  }
  
  if (summ_updated) {
    legendmat.fill.values <- c(legendmat.fill.values, "diamond_fill_updated" = "cornflowerblue")
    legendmat.fill.labels <- c(legendmat.fill.labels, paste0("Updated Pooled Result (diamond - ", round((1-sig_level)*100, 1), "% CI)"))
  }
  
  if (!is.na(zero)) {
    legendmat.fill.values <- c(legendmat.fill.values, "nosig_col" = "white", "sigless_col" = "gray91", "sigmore_col" = "gray72")
    legendmat.fill.labels <- c(legendmat.fill.labels, paste0("Non Sig Effect (", sig_level*100, "% level)"), paste0("Sig Effect < NULL (", sig_level*100, "% level)"), paste0("Sig Effect > NULL (", sig_level*100, "% level)"))
  } else if (!is.na(lower) | !is.na(upper)) {
    legendmat.fill.values <- c(legendmat.fill.values, "noclinsig_col" = "white", "clinsig_col" = "gray72")
    if (!is.na(lower)) {
      legendmat.fill.labels <- c(legendmat.fill.labels, paste0("Non Sig Effect (", round((1-sig_level)*100, 1), "%CI crosses ", ifelse(outcome %in% c('OR', 'RR'), exp(lower), lower), ")"), paste0("Sig Effect (", round((1-sig_level)*100, 1), "%CI is higher than ", ifelse(outcome %in% c('OR', 'RR'), exp(lower), lower), ")"))
    } else {
      legendmat.fill.labels <- c(legendmat.fill.labels, paste0("Non Sig Effect (", round((1-sig_level)*100, 1), "%CI crosses ", ifelse(outcome %in% c('OR', 'RR'), exp(upper), upper), ")"), paste0("Sig Effect (", round((1-sig_level)*100, 1), "%CI is lower than ", ifelse(outcome %in% c('OR', 'RR'), exp(upper), upper), ")"))
    }
  } else if (!is.na(est_pos) & !is.na(est_neg)) {
    legendmat.fill.values <- c(legendmat.fill.values, "nosig_col" = "white", "sigless_col" = "gray91", "sigmore_col" = "gray72")
    legendmat.fill.labels <- c(legendmat.fill.labels, "Non clinically significant result", paste0("Clinically negative effect (less than ", ifelse(outcome %in% c('OR', 'RR'), exp(est_neg), est_neg), " )"), paste0("Clinically positive effect (more than ", ifelse(outcome %in% c('OR', 'RR'), exp(est_pos), est_pos), " )"))
  } else if (!is.na(est_pos)) {
    legendmat.fill.values <- c(legendmat.fill.values, "nosig_col" = "white", "sigmore_col" = "gray72")
    legendmat.fill.labels <- c(legendmat.fill.labels, "Non clinically significant result", paste0("Clinically positive effect (more than ", ifelse(outcome %in% c('OR', 'RR'), exp(est_pos), est_pos), " )"))
  } else if (!is.na(est_neg)) {
    legendmat.fill.values <- c(legendmat.fill.values, "nosig_col" = "white", "sigless_col" = "gray91")
    legendmat.fill.labels <- c(legendmat.fill.labels, "Non clinically significant result", paste0("Clinically negative effect (less than ", ifelse(outcome %in% c('OR', 'RR'), exp(est_neg), est_neg), " )"))
  }
  
  # drop rows that are not included (based on inputs)
  legendmat.col <- legendmat.col[!is.na(legendmat.col$labels), ]
  
   }
  #-------------------#
  # Put together plot #
  #-------------------#
  
  if (draw_plot & length(na.omit(SS)) != 0) {
  
    # empty frame
    plot <- ggplot(data = data.frame(x = SS, y = seSS), aes(x = x, y = y)) +
      labs(x = xlab, y = ylab) +
      theme_classic() + theme(aspect.ratio = 1, panel.background = element_rect(colour = "black")) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_reverse(expand = c(0, 0)) +
      coord_cartesian(xlim = xlim, ylim = ylim)
    
    # Exponential x axis ticks if outcome is a ratio (and no ticks given)
    if (outcome %in% c('OR', 'RR') & is.null(expxticks)) {
      plot <- plot +
        scale_x_continuous(breaks = log(c(0.25, 0.5, 1, 2, 4)), labels = c(0.25, 0.5, 1, 2, 4), expand = c(0,0))
    }
  
    # Specify axis ticks if specified
    # x axis ticks for exponential effects
    if (!is.null(expxticks)) {
      plot <- plot +
        scale_x_continuous(breaks = log(expxticks), labels = expxticks, expand = c(0, 0))
    }
    # x axis ticks (non exp)
    if (!is.null(xticks)) {
    plot <- plot +
        scale_x_continuous(breaks = xticks, labels = xticks, expand = c(0, 0))
    }
    # y axis ticks
    if (!is.null(yticks)) {
      plot <- plot +
        scale_y_reverse(breaks = yticks, labels = yticks, expand = c(0, 0))
    }
  
    # contours
  if (method %in% c('EE', 'MH')) {
    if (!is.na(zero)) {
        plot <- plot +
          geom_polygon(
            data = data.frame(x = c(c1SS, rev(c2SS)), y = c(csize, rev(csize))),
            aes(x = x, y = y, fill = "nosig_col"),
            color = "white"
          ) +
          geom_polygon(
            data = data.frame(x = c(c2SS, xlim[1], xlim[1]), y = c(csize, ylim[2], ylim[1])),
            aes(x = x, y = y, fill = "sigless_col"),
            color = "gray91"
          ) +
          geom_polygon(
            data = data.frame(x = c(c1SS, xlim[2], xlim[2]), y = c(csize, ylim[2], ylim[1])),
            aes(x = x, y = y, fill = "sigmore_col"),
            color = "gray72"
          )
    } else if (!is.na(lower)) {
      plot <- plot +
        geom_polygon(
          data = data.frame(x = c(xlim[1], c1SS, xlim[2], xlim[1]), y = c(ylim[1], csize, ylim[2], ylim[2])),
          aes(x = x, y = y, fill = "noclinsig_col"),
          color = "white"
        ) +
        geom_polygon(
          data = data.frame(x = c(c1SS, xlim[2], xlim[2]), y = c(csize, ylim[2], ylim[1])),
          aes(x = x, y = y, fill = "clinsig_col"),
          color = "gray72"
        )
    } else if (!is.na(upper)) {
      plot <- plot +
        geom_polygon(
          data = data.frame(x = c(xlim[2], c2SS, xlim[1], xlim[2]), y = c(ylim[1], csize, ylim[2], ylim[2])),
          aes(x = x, y = y, fill = "noclinsig_col"),
          color = "white"
        ) +
        geom_polygon(
          data = data.frame(x = c(c2SS, xlim[1], xlim[1]), y = c(csize, ylim[2], ylim[1])),
          aes(x = x, y = y, fill = "clinsig_col"),
          color = "gray72"
        )
    } else if (!is.na(est_pos) | !is.na(est_neg)) {
      if (!is.na(est_pos)) {
        plot <- plot +
          geom_polygon(
            data = data.frame(x = c(c1SS, xlim[2], xlim[2]), y = c(csize, ylim[2], ylim[1])),
            aes(x = x, y = y, fill = "sigmore_col"),
            color = "gray72"
          )
      }
      if (!is.na(est_neg)) {
        plot <- plot +
          geom_polygon(
            data = data.frame(x = c(c2SS, xlim[1], xlim[1]), y = c(csize, ylim[2], ylim[1])),
            aes(x = x, y = y, fill = "sigless_col"),
            color = "gray91"
          )
      }
      if (is.na(est_neg)) {
        plot <- plot +
          geom_polygon(
            data = data.frame(x = c(xlim[1], c1SS, xlim[2], xlim[1]), y = c(ylim[1], csize, ylim[2], ylim[2])),
            aes(x = x, y = y, fill = "nosig_col"),
            color = "white"
          )
      } else if (is.na(est_pos)) {
        plot <- plot +
          geom_polygon(
            data = data.frame(x = c(xlim[2], c2SS, xlim[1], xlim[2]), y = c(ylim[1], csize, ylim[2], ylim[2])),
            aes(x = x, y = y, fill = "nosig_col"),
            color = "white"
          )
      } else {
        plot <- plot +
          geom_polygon(
            data = data.frame(x = c(c1SS, rev(c2SS)), y = c(csize, rev(csize))),
            aes(x = x, y = y, fill = "nosig_col"),
            color = "white"
          )
      }
    } 
  }
    if (method == 'DL') {
      plot <- plot +
        geom_raster(data = contour_tiles, 
          aes(x = cSS, y = csize, fill = code))
    }
  
    # Pooled effect lines
    if (plot_summ_current & length(na.omit(SS)) != 0) {
      plot <- plot +
        geom_vline(aes(xintercept = current_meta$b, color = "summ_current_col"))
    }
    
    if (plot_summ_updated) {
      plot <- plot +
        geom_vline(aes(xintercept = updated_meta$b, color = "summ_updated_col"))
    }
  
    # summary diamonds
    if (summ_current) {
      if (pred_interval) {
        plot <- plot +
          geom_segment(
            aes(
              x = predint1_current,
              y = ylim[2] - 0.10 * axisdiff + summ_pos,
              xend = predint2_current,
              yend = ylim[2] - 0.10 * axisdiff + summ_pos,
              color = "pred_col"
            ),
            show_legend = ifelse(plot_summ_current | plot_zero, FALSE, TRUE)
          )   # needed to avoid crosshairs in legend for when there are vlines also present
      }
      plot <- plot +
        geom_polygon(data = summary_diamond_current, 
                     aes(x = xsumm, y = ysumm, fill = "diamond_fill_current"), 
                     color = "black", alpha = 0.8)
    }
    
    if (summ_updated) {
      if (pred_interval) {
        plot <- plot +
          geom_segment(
            aes(
              x = predint1_updated,
              y = ylim[2] - 0.10 * axisdiff + summ_pos,
              xend = predint2_updated,
              yend = ylim[2] - 0.10 * axisdiff + summ_pos,
              color = "pred_col"
            ),
            show_legend = ifelse(plot_summ_updated | plot_zero, FALSE, TRUE)
          )   # needed to avoid crosshairs in legend for when there are vlines also present
      }
      plot <- plot +
        geom_polygon(data = summary_diamond_updated, 
                     aes(x = xsumm, y = ysumm, fill = "diamond_fill_updated"), 
                     color = "black", alpha = ifelse(summ_current, 0.4, 0.8))
    }
  
    # Null vertical line
    if (plot_threshold) {
      plot <- plot +
        geom_vline(aes(xintercept = ifelse(!is.na(zero), zero, ifelse(!is.na(lower), lower, ifelse(!is.na(upper), upper, ifelse(!is.na(est_pos) & !is.na(est_neg), c(est_pos, est_neg), ifelse(!is.na(est_pos), est_pos, est_neg))))), 
                       color = "threshold_col"))
    }
  
    # Study points
    if (points) {
      plot <- plot +
        geom_point(aes(color = "point_col"), shape = 19)
    }
    
    if (new_points & length(SSnew) > 1) {
      # Dashed lines connecting new study to 'average' point
      lines_df <- data.frame(
        x = SSnew,
        xend = rep(new_avg_est, length(SSnew)),
        y = seSSnew,
        yend = rep(new_avg_se, length(seSSnew))
      )
      plot <- plot +
        geom_segment(data = lines_df,
                     aes(x = x, xend = xend, y = y, yend = yend),
                     linetype = "dashed", color = "black", linewidth = 0.5)
      # Average point
      plot <- plot +
        geom_point(data = data.frame(x = new_avg_est, y = new_avg_se),
          aes(x = x, y = y, color = "new_point_avg_col"), 
          shape = 24, fill = "lightblue", size = 3)
    }
    
    if (new_points) {
      plot <- plot +
        geom_point(data = data.frame(x = SSnew, y = seSSnew), 
                   aes(x = x, y = y, color = "new_point_col"), 
                   shape = 21, fill = "lightblue")
    }
    
    #Tau2 display
    if (tau2 & length(na.omit(SS)) != 0) {
      if (method == 'DL') {
        tau2_label <- paste0(
          "atop(", # atop() stacks the new lines
          "Current~tau^2==", round(current_tau2, 3), ",",
          "Updated~tau^2==", round(updated_tau2, 3),
          ")"
        )
      
        plot <- plot +
          annotate(
            geom = "label",
            x = xlim[1] + 0.05 * (xlim[2] - xlim[1]),
            y = ylim[2] - 0.85 * (ylim[2] - ylim[1]),
            label = tau2_label,
            parse = TRUE,
            hjust = 0,
            vjust = 1,
            size = 3
          )
      }
    }
  
    # Add legend
    if (!is.null(legendmat.col.values)) {
      plot <- plot +
        scale_colour_manual(name = "",
                            values = legendmat.col.values,
                            labels = legendmat.col$labels,
                            guide = guide_legend(override.aes = list(linetype = legendmat.col$linetype,
                                                                     shape = legendmat.col$shape,
                                                                     color = legendmat.col$color,
                                                                     fill = legendmat.col$fill)))
    }
    if (!is.null(legendmat.fill.values)) {
      plot <- plot +
        scale_fill_manual(name = "",
                          values = legendmat.fill.values,
                          labels = legendmat.fill.labels,
                          breaks = names(legendmat.fill.values))
    }
    if (!legend) {
      plot <- plot +
        theme(legend.position = "none")   # turns off legend
    }
    if (legend & !is.null(legendpos)) {
      plot <- plot +
        theme(legend.position = legendpos)   # if user wants legend and specified a position
    }
  
    # return final plot
    threshold_plot <- plot
    
  } else {
    threshold_plot <- "'draw_plot' was not selected to be TRUE or there was no estimable data in the 'current' version"  # message for when users want the plot but didn't select for it
  }
  
  return(list(threshold_result = threshold_result,
              threshold_plot = threshold_plot,
              original_ma = current_meta,
              new_ma = updated_meta))
  
}

#-----------------------------------------#
#             End of function             #
#-----------------------------------------#
