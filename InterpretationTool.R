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


### Test Data antibiotics vs. control for the common cold to alleviate symptoms by 7 days ###
raw_data_bin <- data.frame(StudyID = c(1, 2, 3, 4, 5, 6), Study = c("Herne_1980", "Hoaglund_1950", "Kaiser_1996", "Lexomboon_1971", "McKerrow_1961", "Taylor_1977"),
                   R.1 = c(7, 39, 97, 8, 5, 12), N.1 = c(7+39, 39+115, 97+49, 8+166, 5+10, 12+117), T.1 = rep("Treatment", 6),
                   R.2 = c(10, 51, 94, 4, 8, 3), N.2 = c(10+12, 51+104, 94+48, 4+83, 8+10, 3+56), T.2 = rep("Control", 6))
raw_data_new <- data.frame(StudyID = c(7), Study = c("Fake_1"),
                           R.1 = c(63), N.1 = c(63+252), T.1 = rep("Treatment", 1),
                           R.2 = c(77), N.2 = c(77+105), T.2 = rep("Control", 1))
raw_data_new_multiple <- data.frame(StudyID = c(7, 8, 9), Study = c("Fake_1", "Fake_2", "Fake_3"),
                           R.1 = c(63, 12, 45), N.1 = c(63+252, 12+34, 45+60), T.1 = rep("Treatment", 3),
                           R.2 = c(77, 21, 58), N.2 = c(77+105, 21+25, 58+42), T.2 = rep("Control", 1))
# raw_data_new <- data.frame(StudyID = c(7, 8, 9), Study = c("Fake_1", "Fake_2", "Fake_3"),
#                            R.1 = c(9, 42, 90), N.1 = c(9+36, 42+105, 90+51), T.1 = rep("Treatment", 3),
#                            R.2 = c(11, 48, 99), N.2 = c(11+15, 48+88, 99+62), T.2 = rep("Control", 3))
# raw_data_con <- data.frame(StudyID = c(1, 2, 3, 4, 5, 6, 7), Study = c("Connor_2002", "Geier_2004", "Kinzler_1991", "Lehri_2004", "Halsch_2001", "Volz_1997", "Warnecke_1991"),
#                        Mean.1 = c(5.7, 12.7, 12.3, 10.6, 3, 21, 25.61), Mean.2 = c(8.5, 12.3, 3.6, 9.2, 0.6, 16.2, 7.65),
#                        SD.1 = c(7.6, 6.7, 8.7, 7.3, 7.5, 13, 12.8), SD.2 = c(4.2, 7.3, 8.4, 10, 4.6, 14.3, 15.9),
#                        N.1 = c(17, 25, 29, 34, 20, 52, 20), N.2 = c(18, 25, 29, 23, 30, 48, 20),
#                        T.1 = rep("Treatment", 7), T.2 = rep("Control", 7))
### Obtain study effects and standard errors #
MAdata_bin <- metafor::escalc(measure = 'RR', ai = R.1, bi = N.1-R.1, ci = R.2, di = N.2-R.2, data = raw_data_bin)   # gives ES (effect estimate) and seES (sampling variances) on logOR scale for binary data
MAdata_bin_new <- metafor::escalc(measure = 'RR', ai = R.1, bi = N.1-R.1, ci = R.2, di = N.2-R.2, data = raw_data_new)
MAdata_bin_new_multiple <- metafor::escalc(measure = 'RR', ai = R.1, bi = N.1-R.1, ci = R.2, di = N.2-R.2, data = raw_data_new_multiple)
# MAdata_con <- escalc(measure = "MD", m1i = Mean.1, m2i = Mean.2, sd1i = SD.1, sd2i = SD.2, n1i = N.1, n2i = N.2, data = raw_data_con)
MAdata_bin$sei <- sqrt(MAdata_bin$vi)  # Calculate standard errors
MAdata_bin_new$sei <- sqrt(MAdata_bin_new$vi)
MAdata_bin_new_multiple$sei <- sqrt(MAdata_bin_new_multiple$vi)
# MAdata_con$sei <- sqrt(MAdata_con$vi)


#' @param SS effect estimates of the current studies		
#' @param seSS standard errors of the current studies
#' @param SSnew effect estimates of the new study(ies)
#' @param seSSnew standard errors of the new study(ies)
#' @param sig_level significance level
#' @param method type of meta-analysis - "fixed" or "random"
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
#' @param rand.load TRUE/FALSE show percentage of computations complete when the random effects contours are calculated
#' @return List containing the following: 
#' threshold_result - a string containing a description of whether the threshold of the new analysis was met;
#' threshold_plot - the extended funnel plot that visualises the threshold and whether it was met               
InterpretationThreshold <- function(
    SS,
    seSS,
    SSnew,
    seSSnew,
    sig_level = 0.05,
    method = 'fixed',
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
    rand_load = 100
) {
  
  progressr::handlers(global = TRUE)
  
  
  #----------------------------------------#
  # Calculate initial arguments/parameters #
  #----------------------------------------#
  
  #Convert the significance level into 'z' from the normal distribution (i.e. 0.05 -> 1.96)
  ci <- qnorm(1 - sig_level / 2)
  
  #Calculate current meta-analysis (using rma from {metafor})
  current_meta <- metafor::rma(yi = SS, sei = seSS, method = ifelse(method == 'random', "DL", "EE"), level = (1 - sig_level), measure = outcome)
  
  #Calculate updated meta-analysis (using rma from {metafor})
  updated_meta <- metafor::rma(yi = c(SS, SSnew), sei = c(seSS, seSSnew), method = ifelse(method == 'random', "DL", "EE"), level = (1 - sig_level), measure = outcome)
  updated_tau2 <- updated_meta$tau2
  
  
  #---------------------------#
  # NEW WEIGHTINGS of studies #
  #---------------------------#
  
  # Weights are calculated the same whether 'current' or 'new'
  if (method == "random") {
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
  
  sediff <- max(c(seSS, seSSnew)) - min(c(seSS, seSSnew))
  
  if (!is.null(ylim) && ylim[1] < ylim[2]) {
    ylim <- rev(ylim)  # for when the user has already defined y limits
  }
  
  if (is.null(ylim)) {
    ylim <- c(max(c(seSS, seSSnew)) + 0.20 * sediff, min(c(seSS, seSSnew)) - 0.25 * sediff)
    if (ylim[2] < 0) {
      ylim[2] <- 0
    }
  }
  
  axisdiff <- ylim[2] - ylim[1]
  
  #-----------------------------------#
  # Intelligent x-axis default limits #
  #-----------------------------------#
  
  # log if user specified and outcome is a ratio
  if (!is.null(xlim) & outcome %in% c('OR', 'RR')) {
    xlim <- log(xlim)
  }
  
  # Default limits
  SSdiff <- max(c(SS, SSnew)) - min(c(SS, SSnew))
  
  if (is.null(xlim)) {
    xlim <- c(min(c(SS, SSnew)) - 0.2 * SSdiff, max(c(SS, SSnew)) + 0.2 * SSdiff)
  }

  
  #-------------------------------------------------#
  # Vector of weights/sizes to define contours upon #
  #-------------------------------------------------#
  
   cSS <- seq(xlim[1], xlim[2], length.out = contour_points)  # granulated vector for effect size (x-axis)
   csize <- seq(ylim[1], ylim[2], length.out = contour_points)  # granulated vector for standard error (y-axis)
   csize[csize <= 0] <- 0.0000001 * min(c(seSS, seSSnew))
   for (k in 2:length(csize)) {
     if (csize[k] == 0 & csize[k-1] == 0) {
       csize[k] <- NA
     }
   }
   csize <- csize[!is.na(csize)]   # remove unnecessary data points
  
  #------------------------------#
  # Create significance contours #
  #------------------------------#
   
  # Transform user-specified thresholds if outcome is a ratio
  if (outcome %in% c('OR', 'RR')) {
   if (!is.na(zero)) {zero <- log(zero)}
   else if (!is.na(lower)) {lower <- log(lower)}
   else if (!is.na(upper)) {upper <- log(upper)}
  }
  
  # fixed-effect model #
  if (method == "fixed")  {
    vwt <- 1 / (csize^2)     # weight for each point on plot (inverse of variance) (i.e. weight of new study)
    if (!is.na(zero)) {
      c1SS <- (1 / vwt) * (zero * (sum(size_current) + vwt) - sum(size_current * SS) +  ci * (sum(size_current) + vwt)^0.5)   # formulae based on CI boundaries of new MA meeting no effect
      c2SS <- (1 / vwt) * (zero * (sum(size_current) + vwt) - sum(size_current * SS) -  ci * (sum(size_current) + vwt)^0.5)
    } else if (!is.na(lower)) {
      c1SS <- (1 / vwt) * (lower * (sum(size_current) + vwt) - sum(size_current * SS) +  ci * (sum(size_current) + vwt)^0.5)
    } else if (!is.na(upper)) {
      c2SS <- (1 / vwt) * (upper * (sum(size_current) + vwt) - sum(size_current * SS) -  ci * (sum(size_current) + vwt)^0.5)
    }
  }
    
  # random-effects model #
  if (method == "random")  {
    
    # optimised grid of all contour points
    contour_tiles <- expand.grid(i = seq_along(cSS), j = seq_along(csize))
    
    progressr::with_progress({
      p <- progressr::progressor(steps = nrow(contour_tiles))  # progress bar
    
      # calculate results for each grid point
      results <- purrr::pmap(contour_tiles, function(i, j) {
        # using rmeta::meta.summaries as it's much faster than metafor::rma
        metacont <- rmeta::meta.summaries(d = c(SS, cSS[i]), se = c(seSS, csize[j]), method = "random", conf.level = (1-sig_level))
        lc <- metacont$summary - ci*metacont$se.summary
        uc <- metacont$summary + ci*metacont$se.summary
        
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
        }
        NA_character_
      })
      
      # Assign values to contour_tiles
      contour_tiles$code <- unlist(results)
      
    })
    
    contour_tiles <- cbind(contour_tiles, expand.grid(cSS = cSS, csize = csize))
  }
   
  #------------------------------#
  # Assess whether threshold met #
  #------------------------------#  
   
  # Acquire updated meta-analysis confidence intervals
  new_pvalue <- updated_meta$pval 
  new_lower <- updated_meta$ci.lb
  new_upper <- updated_meta$ci.ub
  
  # Where focus is on statistical significance (where alpha = 0.05)
  if (sig_level == 0.05 & !is.na(zero) & is.na(lower) & is.na(upper)) {
    if (new_pvalue < 0.05) {
      threshold_result <- paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will give a statistically significant (alpha = 0.05) pooled estimate")
    } else {
      threshold_result <- paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will not give a statistically significant (alpha = 0.05) pooled estimate")
    }
  # Where focus is on statistical significance (where alpha != 0.05)
  } else if (sig_level != 0.05 & !is.na(zero) & is.na(lower) & is.na(upper)) {
    if (new_upper < 0 | new_lower > 0) {
      threshold_result <- paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will give a statistically significant (alpha = ", sig_level, ") pooled estimate")
    } else {
      threshold_result <- paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will not give a statistically significant (alpha = ", sig_level, ") pooled estimate")
    }
  # Where focus is on a specific upper or lower confidence band
  } else if (is.na(zero) & !is.na(lower) & is.na(upper)) {
    if (new_lower > lower) {
      threshold_result <- paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will give a ", round(100*(1 - sig_level),1), "% CI that is higher than ", ifelse(outcome %in% c('OR', 'RR'), exp(lower), lower))
    } else {
      threshold_result <- paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will not give a ", round(100*(1 - sig_level),1), "% CI that is higher than ", ifelse(outcome %in% c('OR', 'RR'), exp(lower), lower))
    }
  } else if (is.na(zero) & is.na(lower) & !is.na(upper)) {
    if (new_upper < upper) {
      threshold_result <- paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will give a ", round(100*(1 - sig_level),1), "% CI that is lower than ", ifelse(outcome %in% c('OR', 'RR'), exp(upper), upper))
    } else {
      threshold_result <- paste0("Addition of new ", ifelse(length(SSnew) > 1, "studies", "study"), " will not give a ", round(100*(1 - sig_level),1), "% CI that is lower than ", ifelse(outcome %in% c('OR', 'RR'), exp(upper), upper))
    }
  } else {
    print ("(Only) one of zero, lower, or upper need to be given a value")
  }
  
  #------------------#
  # Summary diamonds #
  #------------------#
  
  if (summ_current) {
    summary_diamond_current <- data.frame(
      xsumm = c(current_meta$ci.lb, current_meta$b, current_meta$ci.ub, current_meta$b),
      ysumm = c(ylim[2] - 0.10 * axisdiff + summ_pos, ylim[2] - 0.07 * axisdiff + summ_pos, 
                ylim[2] - 0.10 * axisdiff + summ_pos, ylim[2] - 0.13 * axisdiff + summ_pos)
    )
    
    if (pred_interval) {	
      if (method == 'random') {
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
       if (method == 'random') {
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
    legendmat.col$shape[3] <- 23
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
  }
  
  # drop rows that are not included (based on inputs)
  legendmat.col <- legendmat.col[!is.na(legendmat.col$labels), ]
  
  #-------------------#
  # Put together plot #
  #-------------------#
  
  if (draw_plot) {
  
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
  if (method == 'fixed') {
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
    }
  }
    if (method == 'random') {
      plot <- plot +
        geom_raster(data = contour_tiles, 
          aes(x = cSS, y = csize, fill = code))
    }
  
    # Pooled effect lines
    if (plot_summ_current) {
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
        geom_vline(aes(xintercept = ifelse(!is.na(zero), zero, ifelse(!is.na(lower), lower, upper)), color = "threshold_col"))
    }
  
    # Study points
    if (points) {
      plot <- plot +
        geom_point(aes(color = "point_col"), shape = 19)
    }
    
    if (new_points) {
      plot <- plot +
        geom_point(data = data.frame(x = SSnew, y = seSSnew), 
                   aes(x = x, y = y, color = "new_point_col"), 
                   shape = 21, fill = "lightblue")
    }
    
    if (new_points & length(SSnew) > 1) {
      plot <- plot +
        geom_point(data = data.frame(x = new_avg_est, y = new_avg_se),
          aes(x = x, y = y, color = "new_point_avg_col"), 
          shape = 23, fill = "lightblue", size = 3)
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
    threshold_plot <- "'draw_plot' was not selected to be TRUE"  # message for when users want the plot but didn't select for it
  }
  
  return(list(threshold_result = threshold_result,
              threshold_plot = threshold_plot))
  
}

#-----------------------------------------#
#             End of function             #
#-----------------------------------------#

# Tests #

# Study points, new points, both summary diamonds, and zero line
result <- InterpretationThreshold(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, 
  SSnew = MAdata_bin_new$yi, seSSnew = MAdata_bin_new$sei,
  method = 'fixed', outcome = 'RR', zero = 1)

# Above but with different significance level
result <- InterpretationThreshold(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, 
  SSnew = MAdata_bin_new$yi, seSSnew = MAdata_bin_new$sei,
  method = 'fixed', outcome = 'RR', zero = 1, sig_level = 0.0025)

# Threshold is now based on lower bound of confidence interval
result <- InterpretationThreshold(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, 
  SSnew = MAdata_bin_new$yi, seSSnew = MAdata_bin_new$sei,
  method = 'fixed', outcome = 'RR', lower = 1.1)

# Threshold is now based on upper bound of confidence interval
result <- InterpretationThreshold(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, 
  SSnew = MAdata_bin_new$yi, seSSnew = MAdata_bin_new$sei,
  method = 'fixed', outcome = 'RR', upper = 0.9)

# Multiple new studies
result <- InterpretationThreshold(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, 
  SSnew = MAdata_bin_new_multiple$yi, seSSnew = MAdata_bin_new_multiple$sei,
  method = 'fixed', outcome = 'RR', zero = 1, sig_level = 0.0005)

# Random effects with one study and statistical significance
result <- InterpretationThreshold(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, 
  SSnew = MAdata_bin_new$yi, seSSnew = MAdata_bin_new$sei,
  method = 'random', outcome = 'RR', zero = 1, sig_level = 0.1, contour_points = 100)

# Random effects with one study and clinical significance
result <- InterpretationThreshold(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, 
  SSnew = MAdata_bin_new$yi, seSSnew = MAdata_bin_new$sei,
  method = 'random', outcome = 'RR', upper = 0.95, sig_level = 0.15, contour_points = 100)
