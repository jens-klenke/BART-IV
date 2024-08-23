### EVALUATION FUNCTIONS
bias_fun <- function(y_hat, y, ...) mean(y_hat -y) 
PEHE_fun <- function(y_hat, y, ...) mean( (y_hat - y)^2)
coverage_fun <- function(...) 'Coveage needs to be implementet'

## https://github.com/albicaron/SparseBCF/blob/main/Simulated%20examples/Section%205.1/AllR_Models_P25.R