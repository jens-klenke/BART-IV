### EVALUATION FUNCTIONS
bias_fun <- function(tau_true, tau_pred, ...) mean(tau_true - tau_pred) 
PEHE_fun <- function(tau_true, tau_pred, ...) mean( (tau_true - tau_pred)^2)
coverage_fun <- function(...) 'Coveage needs to be implementet'

## https://github.com/albicaron/SparseBCF/blob/main/Simulated%20examples/Section%205.1/AllR_Models_P25.R