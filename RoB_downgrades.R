#---------------------------------------------------------#
# Tools for estimating any downgrading due to risk of bias #
#---------------------------------------------------------#

library(dplyr)

#' @param data Dataset containing information
#' @param rob_tool What type of risk of bias tool was used (1 or 2)
#' @param overall_rob Overall rating of risk of bias
#' @param random_selection Random sequence generation (selection bias) rating (RoB 1 only)
#' @param allocation_selection Allocation concealment (selection bias) rating (RoB 1 only)
#' @param performance Performance bias rating (RoB 1 only)
#' @param detection Detection bias rating (RoB 1 only)
#' @param attrition Attrition bias rating (RoB 1 only)
#' @param reporting Reporting bias rating (RoB 1 only)
#' @param other Other bias rating (RoB 1 only)
#' @param weights Meta-analysis weighting for each study
#' @return RoB_avg = weighted RoB value (1 = all low, 3 = all high)

weighted_RoB <- function(
    data, 
    rob_tool,
    overall_rob,
    random_selection = NA,
    allocation_selection = NA,
    performance = NA,
    detection = NA,
    attrition = NA,
    reporting = NA,
    other = NA,
    weights
) {
  
  #-------------------------------------------------------------------------------------#
  # Where they've used RoB 1 tool and there is no overall rating, calculate an estimate #
  #-------------------------------------------------------------------------------------#
  
  if (rob_tool == 1 & max(is.na(data[[overall_rob]])) == 1) {
    
    # Read in transformation matrices
    selection_bias_matrix <- readxl::read_excel("RoB matrices.xlsx", sheet = "selection_bias")
    overall_bias_matrix <- readxl::read_excel("RoB matrices.xlsx", sheet = "overall bias")
    
    # clean risk of bias answers
    selection_bias_matrix <- mutate_all(selection_bias_matrix, .funs = tolower)
    data <- data %>%
      mutate(across(all_of(c(random_selection, allocation_selection, performance, detection, attrition, reporting, other, overall_rob)),
                    ~ gsub(" risk", "", stringr::str_to_lower(.x))))
    
    
    # calculate selection bias probabilities
    # probability of selection bias being 'high'
    data$selection_high <- mapply(function(random_data, allocation_data) {
      idx <- which(selection_bias_matrix$`Random sequence generation` == random_data & 
                   selection_bias_matrix$`Allocation concealment` == allocation_data)
      if (length(idx) > 0) as.numeric(selection_bias_matrix$High[idx[1]]) else NA}, 
    data[[random_selection]], data[[allocation_selection]])
    # probability of selection bias being 'low'
    data$selection_low <- mapply(function(random_data, allocation_data) {
      idx <- which(selection_bias_matrix$`Random sequence generation` == random_data & 
                     selection_bias_matrix$`Allocation concealment` == allocation_data)
      if (length(idx) > 0) as.numeric(selection_bias_matrix$Low[idx[1]]) else NA}, 
      data[[random_selection]], data[[allocation_selection]])
    #probability of selection bias being 'unclear'
    data$selection_unclear <- mapply(function(random_data, allocation_data) {
      idx <- which(selection_bias_matrix$`Random sequence generation` == random_data & 
                     selection_bias_matrix$`Allocation concealment` == allocation_data)
      if (length(idx) > 0) as.numeric(selection_bias_matrix$Unclear[idx[1]]) else NA}, 
      data[[random_selection]], data[[allocation_selection]])
    
    # calculate probabilities of overall rating being low, some concerns, or high
    # probability being low
    data$overall_low <- mapply(function(selection_high_prob, selection_low_prob, selection_unclear_prob,
                                        performance_data, detection_data, attrition_data, reporting_data, other_data) {
      # Helper function: get match from bias matrix
      get_val <- function(sel_bias) {
        idx <- which(overall_bias_matrix$`Selection bias` == sel_bias &
                       performance_data == overall_bias_matrix$`Performance bias` &
                       detection_data == overall_bias_matrix$`Detection bias` &
                       attrition_data == overall_bias_matrix$`Attrition bias` &
                       reporting_data == overall_bias_matrix$`Reporting bias` &
                       other_data == overall_bias_matrix$Other)
        if (length(idx) > 0) {
          overall_bias_matrix$`Low risk`[idx[1]]
        } else {
          NA_real_
        }
      }
      
      # Get results for each type
      result_low     <- get_val("low")
      result_high    <- get_val("high")
      result_unclear <- get_val("unclear")
      
      # Weighted sum
      selection_high_prob    * result_high +
        selection_low_prob   * result_low +
        selection_unclear_prob * result_unclear
    },
    data$selection_high, data$selection_low, data$selection_unclear,
    data[[performance]], data[[detection]], data[[attrition]], data[[reporting]], data[[other]])
    # probability being some
    data$overall_some <- mapply(function(selection_high_prob, selection_low_prob, selection_unclear_prob,
                                        performance_data, detection_data, attrition_data, reporting_data, other_data) {
      # Helper function: get match from bias matrix
      get_val <- function(sel_bias) {
        idx <- which(overall_bias_matrix$`Selection bias` == sel_bias &
                       performance_data == overall_bias_matrix$`Performance bias` &
                       detection_data == overall_bias_matrix$`Detection bias` &
                       attrition_data == overall_bias_matrix$`Attrition bias` &
                       reporting_data == overall_bias_matrix$`Reporting bias` &
                       other_data == overall_bias_matrix$Other)
        if (length(idx) > 0) {
          overall_bias_matrix$`Some concerns`[idx[1]]
        } else {
          NA_real_
        }
      }
      
      # Get results for each type
      result_low     <- get_val("low")
      result_high    <- get_val("high")
      result_unclear <- get_val("unclear")
      
      # Weighted sum
      selection_high_prob    * result_high +
        selection_low_prob   * result_low +
        selection_unclear_prob * result_unclear
    },
    data$selection_high, data$selection_low, data$selection_unclear,
    data[[performance]], data[[detection]], data[[attrition]], data[[reporting]], data[[other]])
    # probability being high
    data$overall_high <- mapply(function(selection_high_prob, selection_low_prob, selection_unclear_prob,
                                        performance_data, detection_data, attrition_data, reporting_data, other_data) {
      # Helper function: get match from bias matrix
      get_val <- function(sel_bias) {
        idx <- which(overall_bias_matrix$`Selection bias` == sel_bias &
                       performance_data == overall_bias_matrix$`Performance bias` &
                       detection_data == overall_bias_matrix$`Detection bias` &
                       attrition_data == overall_bias_matrix$`Attrition bias` &
                       reporting_data == overall_bias_matrix$`Reporting bias` &
                       other_data == overall_bias_matrix$Other)
        if (length(idx) > 0) {
          overall_bias_matrix$`High risk`[idx[1]]
        } else {
          NA_real_
        }
      }
      
      # Get results for each type
      result_low     <- get_val("low")
      result_high    <- get_val("high")
      result_unclear <- get_val("unclear")
      
      # Weighted sum
      selection_high_prob    * result_high +
        selection_low_prob   * result_low +
        selection_unclear_prob * result_unclear
    },
    data$selection_high, data$selection_low, data$selection_unclear,
    data[[performance]], data[[detection]], data[[attrition]], data[[reporting]], data[[other]])
    
    # Add to overall column most likely rating
    data <- data %>%
      mutate(!!overall_rob := case_when(
        overall_low == pmax(overall_low, overall_some, overall_high) ~ "low",
        overall_some == pmax(overall_low, overall_some, overall_high) ~ "some concerns",
        overall_high == pmax(overall_low, overall_some, overall_high) ~ "high"
      ))
    
  } else if (rob_tool == 2 & max(is.na(data[[overall_rob]])) == 1) {
    print("A value is needed for the overall_rob column")
  } else if (max(is.na(data[[overall_rob]])) == 0 & rob_tool %in% c(1,2)) {
    cols <- c(random_selection, allocation_selection, performance, detection, attrition, reporting, other, overall_rob)
    cols <- cols[!is.na(cols)]   # remove any NA column names
    data <- data %>%
      mutate(across(any_of(cols),
                    ~ gsub(" risk", "", stringr::str_to_lower(.x))))
  } 
  else {
    print("The rob_tool column needs to be the numeric value 1 or 2")
  }
  
  
  #-----------------------------------------#
  # Calculate weighted average risk of bias #
  #-----------------------------------------#
  
  # Assign values to risk of bias ratings
  lookup <- c("not applicable" = 0,
              "low" = 1,
              "some concerns" = 2,
              "high" = 3)
  data <- data %>%
    mutate(RoB_value = lookup[data[[overall_rob]]])
  
  # Calculate weighted average
  data <- data %>%
    mutate(weighted_RoB = RoB_value * as.numeric(data[[weights]]))
  weighted_average <- sum(data[["weighted_RoB"]]) / sum(as.numeric(data[[weights]][data[["RoB_value"]] != 0]))
  
  
  #----------------#
  # Return RoB_avg #
  #----------------#
  
  return(unname(weighted_average))
  
}

#' @param RoB_avg Weighted average RoB value calculated from weighted_RoB
#' @param one_level_threshold Threshold value of weighted RoB value (1 = all low, 3 = all high) that means downgrading by one level
#' @param two_level_threshold Threshold value of weighted RoB value (1 = all low, 3 = all high) that means downgrading by two levels
#' @return levels = Number of levels the evidence is likely to be downgraded due to risk of bias (0, 1, or 2)

RoB_downgrades <- function(
    RoB_avg,
    one_level_threshold = 1.5,
    two_level_threshold = NULL
) {
  
  #---------------------------------------------------------------#
  # Decide how many levels to downgrade based on weighted average #
  #---------------------------------------------------------------#
  
  downgrade_levels <- 
    if (!is.null(two_level_threshold)) {
      if (RoB_avg < one_level_threshold) {
        0
      } else if (RoB_avg >= one_level_threshold & RoB_avg < two_level_threshold) {
        1
      } else if (RoB_avg >= two_level_threshold) {
        2
      }
    } else {
      if (RoB_avg < one_level_threshold) {
        0
      } else if (RoB_avg >= one_level_threshold) {
        1
      }
    }
  
  
  return(levels = downgrade_levels)
  
}
