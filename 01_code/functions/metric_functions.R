### EVALUATION FUNCTIONS
bias_fun <- function(tau_true, tau_pred, ...) mean(tau_true - tau_pred) 
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