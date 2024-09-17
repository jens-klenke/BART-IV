### EVALUATION FUNCTIONS
metric_fun_df <- function(leaves_df, ...) {
  leaves_df %>%
    dplyr::mutate(difference = tau_true - tau_pred,
                  sq_difference = abs(tau_true - tau_pred)) %>%
    dplyr::summarise(bias  = mean(difference),
                     bias_rm.na = mean(difference, na.rm = TRUE),
                     abs_bias  = mean(abs(difference)),
                     abs_bias_rm.na = mean(abs(difference), na.rm = TRUE),
                     PEHE  = mean(sq_difference),
                     PEHE_rm.na = mean(sq_difference, na.rm = TRUE),
                     coverage = mean(coverage)
    )
}

# collect all pred and true taus in leaves
leaf_df_fun <- function(results, ...){
  df <- results %>%
    dplyr::filter(leaves == 1) %>%
    dplyr::select(pred) %>%
    tidyr::unnest(pred)
  
  return(df)
}

# check if we FALSE -> estimation problems!
check_fun <- function(results, ...){
  est_problem <- results %>%
    dplyr::select(est_problems) %>%
    dplyr::summarize(any_yes = any(est_problems == "yes", na.rm = TRUE)) %>%
    pull()
  return(est_problem)
}

# function
analysis_fun <- function(results, ...){
  
  # df
  df <- results %>%
    tibble::as_tibble() %>%
    dplyr::select(bcf_results, s_bcf_results) %>%
    tidyr::pivot_longer(cols = c('bcf_results', 's_bcf_results'), names_to = 'detect_model', values_to = 'results') %>%
    dplyr::mutate(methods = names(results)) %>% 
    dplyr::mutate(est_problems = purrr::pmap_lgl(., check_fun),
                  leaves_df = purrr::pmap(., leaf_df_fun)
    ) %>%
    dplyr::arrange(detect_model)
  
  df %<>%
    dplyr::mutate(purrr::pmap_df(., metric_fun_df))
  
  return(df)
  
}


bias_fun <- function(tau_true, tau_pred, ...) mean(tau_true - tau_pred) # not used right now
abs_bias_fun <- function(tau_true, tau_pred, ...) mean(abs(tau_true - tau_pred))
PEHE_fun <- function(tau_true, tau_pred, ...) mean( (tau_true - tau_pred)^2)
MC_se_fun <- function(x, B) qt(0.975, B-1)*sd(x)/sqrt(B)
coverage_95_fun <- function(obj, tau_true, bayes = FALSE, ...){

  # frequentistic
  if(!bayes){
    # estimate 
    est <- obj$coefficients[2, 'Estimate']
    # std error 
    std.error <- obj$coefficients[2, 'Std. Error']
    # degree of fredoom
    df_iv <- obj$df[2]
    
    # confidence interval
    low_ci <- est - std.error * qt(0.975, df_iv) 
    upper_ci <- est + std.error * qt(0.975, df_iv)
    # print(paste('low:', low_ci, '; high:', upper_ci))
    # checking coverage
    cove <- low_ci < tau_true & tau_true < upper_ci

    
  }
  
  if(bayes){
    cove <- obj[2, 3] < tau_true & tau_true < obj[2, 4]
  }
  
  return(cove)
}

## https://github.com/albicaron/SparseBCF/blob/main/Simulated%20examples/Section%205.1/AllR_Models_P25.R