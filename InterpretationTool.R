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


### Test Data antibiotics vs. control for the common cold to alleviate symptoms by 7 days ###
raw_data_bin <- data.frame(StudyID = c(1, 2, 3, 4, 5, 6), Study = c("Herne_1980", "Hoaglund_1950", "Kaiser_1996", "Lexomboon_1971", "McKerrow_1961", "Taylor_1977"),
                   R.1 = c(7, 39, 97, 8, 5, 12), N.1 = c(7+39, 39+115, 97+49, 8+166, 5+10, 12+117), T.1 = rep("Treatment", 6),
                   R.2 = c(10, 51, 94, 4, 8, 3), N.2 = c(10+12, 51+104, 94+48, 4+83, 8+10, 3+56), T.2 = rep("Control", 6))
raw_data_new <- data.frame(StudyID = c(7), Study = c("Fake_1"),
                           R.1 = c(9), N.1 = c(9+36), T.1 = rep("Treatment", 1),
                           R.2 = c(11), N.2 = c(11+15), T.2 = rep("Control", 1))
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
# MAdata_con <- escalc(measure = "MD", m1i = Mean.1, m2i = Mean.2, sd1i = SD.1, sd2i = SD.2, n1i = N.1, n2i = N.2, data = raw_data_con)
MAdata_bin$sei <- sqrt(MAdata_bin$vi)  # Calculate standard errors
MAdata_bin_new$sei <- sqrt(MAdata_bin_new$vi)
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
#' @param plot_zero TRUE/FALSE plot the null effect vertical line as defined by the argument 'zero'	
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
    xlab = "Effect",
    ylab = "Standard Error" ,
    plot_zero = TRUE,
    plot_summ_current = FALSE,
    plot_sum_updated = FALSE,
    legendpos = NULL,
    summ_current = TRUE,
    summ_updated = TRUE,
    summ_pos = 0,
    new_points = TRUE,
    points = TRUE,
    pred_interval = FALSE,
    rand_load = 100
) {
  
  
  #----------------------------------------#
  # Calculate initial arguments/parameters #
  #----------------------------------------#
  
  #Convert the significance level into 'z' from the normal distribution (i.e. 0.05 -> 1.96)
  ci <- qnorm(1 - sig_level / 2)
  
  #Calculate current meta-analysis (using rma from {metafor})
  current_meta <- metafor::rma(yi = SS, sei = seSS, method = ifelse(method == 'random', "PM", "FE"), level = (1 - sig_level), measure = outcome)
  current_tau2 <- current_meta$tau2
  
  #Calculate updated meta-analysis (using rma from {metafor})
  updated_meta <- metafor::rma(yi = c(SS, SSnew), sei = c(seSS, seSSnew), method = ifelse(method == 'random', "PM", "FE"), level = (1 - sig_level), measure = outcome)
  updated_tau2 <- updated_meta$tau2
  
  
  #-------------------------------#
  # CURRENT WEIGHTINGS of studies #
  #-------------------------------#
  
  # Weights are calculated the same whether 'current' or 'new'
  if (method == "random") {
    size_current <- 1 / ((seSS^2) + current_tau2)  # standard inverse-variance weighting
    size_new <- 1 / ((seSSnew^2) + updated_tau2)
  } else {
    size_current <- 1 / (seSS^2)
    size_new <- 1 / (seSSnew^2)
  }
  
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
  if (method == "random")  {     # code had been developed from Langan et al & Florian Teichert
    print("A less computationally expensive method is required.")
    # tibble of every point on plot (i.e. every possible combo)
    #contour_tiles <- expand_grid(cSS = cSS, csize = csize) %>%
    #        mutate(id = row_number(), code = NA) %>%
    #        select(id, everything())
    #for (i in 1:nrow(contour_tiles)) {        # need to increase speed....
    #  if (rand.load>0) {
    #    roundi <- i/rand.load
    #    flush.console()
    #    if (roundi == round(roundi, 0)) {
    #      perc_complete <- (i/(contour.points*contour.points))*100
    #      cat(perc_complete, "%")
    #    } else cat(".")
    #  }
    #  metacont <- rma(yi = c(SS, contour_tiles$cSS[i]), sei = c(seSS, contour_tiles$csize[i]), method = ifelse(method == 'random', "PM", "FE"), level = (1-sig.level))
    #  lc <- metacont$ci.lb
    #  uc <- metacont$ci.ub
    #  # code according to significance
    #  if (lc < zero & uc < zero) contour_tiles$code[i] <- "sigless_col"   # sig < 0
    #  if (lc < zero & uc > zero) contour_tiles$code[i] <- "nosig_col"   # not sig
    #  if (lc > zero & uc > zero) contour_tiles$code[i] <- "sigmore_col"   # sig > 0
    #}
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
  legendmat.col <- data.frame(labels = rep(NA, 6), linetype = rep(NA, 6), shape = rep(NA, 6), color = rep(NA, 6), fill = rep(NA, 6))
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
    legendmat.col$labels[2] <- "New studies"
    legendmat.col$linetype[2] <- "blank"
    legendmat.col$shape[2] <- 23
    legendmat.col$color[2] <- "black"
    legendmat.col$fill[2] <- "blue"
  }
  
  if (pred_interval) {
    legendmat.col.values <- c(legendmat.col.values, "pred_col" = "black")
    legendmat.col$labels[3] <- "95% Predictive Interval"
    legendmat.col$linetype[3] <- "solid"
    legendmat.col$color[3] <- "black"
  }
  
  if (plot_summ_current) {
    legendmat.col.values <- c(legendmat.col.values, "summ_current_col" = "slategrey")
    legendmat.col$labels[4] <- "Current Pooled Effect"
    legendmat.col$linetype[4] <- "solid"
    legendmat.col$color[4] <- "slategrey"
  }
  
  if (plot_summ_updated) {
    legendmat.col.values <- c(legendmat.col.values, "summ_updated_col" = "cadetblue4")
    legendmat.col$labels[5] <- "Updated Pooled Effect"
    legendmat.col$linetype[5] <- "solid"
    legendmat.col$color[5] <- "cadetblue4"
  }
  
  if (plot_zero) {
    legendmat.col.values <- c(legendmat.col.values, "zero_col" = "lightgrey")
    legendmat.col$labels[6] <- "Null Effect"
    legendmat.col$linetype[6] <- "solid"
    legendmat.col$color[6] <- "lightgray"
  }
  
  if (summ_current) {
    legendmat.fill.values <- c(legendmat.fill.values, "diamond_fill_current" = "lavenderblush4")
    legendmat.fill.labels <- c(legendmat.fill.labels, "Current Pooled Result (diamond)")
  }
  
  if (summ_updated) {
    legendmat.fill.values <- c(legendmat.fill.values, "diamond_fill_updated" = "cornflowerblue")
    legendmat.fill.labels <- c(legendmat.fill.labels, "Updated Pooled Result (diamond)")
  }
  
  legendmat.fill.values <- c(legendmat.fill.values, "nosig_col" = "white", "sigless_col" = "gray91", "sigmore_col" = "gray72")
  legendmat.fill.labels <- c(legendmat.fill.labels, paste("Non Sig Effect (", sig_level*100, "% level)", sep = ""), paste("Sig Effect < NULL (", sig_level*100, "% level)", sep = ""), paste("Sig Effect > NULL (", sig_level*100, "% level)", sep = ""))
  
  # drop rows that are not included (based on inputs)
  legendmat.col <- legendmat.col[!is.na(legendmat.col$labels), ]
  
  #-------------------#
  # Put together plot #
  #-------------------#
  
  threshold_plot <- "'draw_plot' was not selected to be TRUE"  # message for when users want the plot but didn't select for it 
  
  if (draw_plot) {
  
    # empty frame
    plot <- ggplot(data = data.frame(x = SS, y = seSS), aes(x = x, y = y)) +
      labs(x = xlab, y = ylab) +
      theme_classic() + theme(aspect.ratio = 1, panel.background = element_rect(colour = "black")) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_reverse(expand = c(0, 0)) +
      coord_cartesian(xlim = xlim, ylim = ylim)
  
    # Specify axis ticks if specified or exponential
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
    }
    if (method == 'random') {
      plot <- plot +
        geom_tile(data = contour_tiles, aes(x = cSS, y = csize, fill = code))  # haven't tested this yet or its affect on the legend
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
    if (plot_zero) {
      plot <- plot +
        geom_vline(aes(xintercept = zero, color = "zero_col"))
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
                   shape = 23, fill = "blue")
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
                          labels = legendmat.fill.labels)
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
  }
  
  return(list(threshold_result = threshold_result,
              threshold_plot = threshold_plot))
  
}

#-----------------------------------------#
#             End of function             #
#-----------------------------------------#

# Tests #

# Study points & summary diamond  PASS
#extfunnel(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, method = 'fixed', outcome = 'OR',
#          ylim = c(0, 1), expxticks = c(0.25, 0.5, 1, 2, 4), xlab = "Odds Ratio", legend = TRUE)

# Study points & summary diamond with predictive interval  PASS
#extfunnel(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, method = 'random', outcome = 'OR',
#          ylim = c(0, 1), expxticks = c(0.25, 0.5, 1, 2, 4), xlab = "Odds Ratio", pred.interval = TRUE, legend = TRUE)

# Study points & summary diamond & effect line PASS
#extfunnel(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, method = 'fixed', outcome = 'OR',
#          ylim = c(0, 1), expxticks = c(0.25, 0.5, 1, 2, 4), xlab = "Odds Ratio", plot.summ = TRUE, legend = TRUE)

# Study points & summary diamond & significance contours  PASS
#extfunnel(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, method = 'fixed', outcome = 'OR',
#          ylim = c(0, 1), xlim = (log(c(0.1, 4))), expxticks = c(0.25, 0.5, 1, 2, 4), xlab = "Odds Ratio", contour = TRUE, legend = TRUE, legendpos = 'left')

# Study points, summary diamond, significance contours & simulated trials  PASS
#extfunnel(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, method = 'fixed', outcome = 'OR',
#          ylim = c(0, 1), expxticks = c(0.25, 0.5, 1, 2, 4), xlab = "Odds Ratio", contour = TRUE, sim.points = sims_bin$sim_study, legend = TRUE)

# Above but zoomed in on simulated studies  PASS
#extfunnel(SS = MAdata_bin$yi, seSS = MAdata_bin$sei, method = 'fixed', outcome = 'OR',
#          ylim = c(0.05, 0.15), xlim = log(c(0.4, 2.1)), expxticks = c(0.25, 0.5, 1, 2, 4), xlab = "Odds Ratio",
#          contour = TRUE, sim.points = sims_bin$sim_study, legend = TRUE)   # would be ideal to remove 'current studies' from legend if they are forced off from the plot

# # Mirror figure 2C in Langan et al  FAIL -> yes but big old triangle of colour missing at the bottom (I think it's going to have to stay as a bug for now...its only when the 'curtain' moves across to the opposite side to standard...)
# extfunnel(SS = MAdata_con$yi, seSS = MAdata_con$sei, method = 'fixed', outcome = 'MD',
#           ylim = c(0, 5.5), xlab = "Difference in means", contour = TRUE, plot.summ = TRUE, legend = TRUE)
#
# # Finding error of continuous not showing predictive interval (FIXED)
# extfunnel(SS = MAdata_con$yi, seSS = MAdata_con$sei, method = 'random', outcome = 'MD',
#           ylim = c(0, 5.5), xlab = "Difference in means", contour = FALSE, plot.summ = TRUE, legend = TRUE, pred.interval = TRUE)

