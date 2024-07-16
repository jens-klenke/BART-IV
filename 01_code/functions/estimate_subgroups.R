estimate_subgroups <- function(inference, ...){
  
  # store matrix
  
  subgroup_estimates <- as.data.frame(matrix(NA, nrow = 2, ncol =7))
  names(subgroup_estimates) <- c("node", "CCACE", "pvalue", "Weak_IV_test", 
                       "Pi_obs", "ITT", "Pi_compliers")
  
  # subset -----
  positive_group_data <- inference %>%
    dplyr::filter(x1 == 0 & x2 == 0)
  
  negative_group_data <- inference %>%
    dplyr::filter(x1 == 1 & x2 == 1)
  
  # Run an IV Regression on each Subgroup ----  
  positive_group <- ivreg(y ~ w | z, data = positive_group_data)
  subgroup_estimates[1, ] <- get_iv_summary(positive_group, inference)
  
  negative_group <- ivreg(y ~ w | z, data = negative_group_data)
  subgroup_estimates[2, ] <- get_iv_summary(negative_group, inference)

  return(subgroup_estimates)
} # close function
