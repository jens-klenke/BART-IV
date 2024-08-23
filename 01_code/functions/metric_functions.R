### EVALUATION FUNCTIONS
bias_fun <- function(y_hat, y, ...) mean(y_hat -y) 
PEHE_fun <- function(y_hat, y, ...) mean( (y_hat - y)^2)
coverage_fun <- function(...) 'Coveage needs to be implementet'

