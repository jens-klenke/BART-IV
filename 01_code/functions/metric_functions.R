### EVALUATION FUNCTIONS
bias_fun <- function(tau_true, tau_pred, ...) mean(tau_true - tau_pred) 
PEHE_fun <- function(tau_true, tau_pred, ...) mean( (tau_true - tau_pred)^2)
MC_se_fun <- function(x, B) qt(0.975, B-1)*sd(x)/sqrt(B)
coverage_95_fun <- function(CCACE_EST, CCACE) {
  quan = apply(CCACE_EST, 2, function(x) quantile(x, c(0.025, 0.975)))
  cove = sum((quan[1,] < CCACE) & (ITE < CCACE[2,]))/length(CCACE)
  
  return(cove)
}


## https://github.com/albicaron/SparseBCF/blob/main/Simulated%20examples/Section%205.1/AllR_Models_P25.R