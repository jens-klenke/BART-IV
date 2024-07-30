library(AER)
library(brms)

set.seed(123456)  # For reproducibility

# Sample size
n <- 1000

# Generate the instrument (Z)
Z <- rnorm(n)

# Generate confounder (U)
U <- rnorm(n)
# Generate the treatment (X), affected by Z and U
X <- 0.5 * Z + 0.3 * U + rnorm(n)
# Generate the outcome (Y), affected by X and U
Y <- 2 * X + 0.5 * U + rnorm(n)

# Combine into a data frame ## needs to be changed! 
data <- data.frame(Y, X, Z, U)

# first stage -----
# first stage: regress X on Z
stan_model_first_stage <- brm(X ~ Z, data = data, chains = 4, iter = 2000, warmup = 1000, 
                              save_model = "05_stan_code/brms_first_stage.stan",  silent = 2, refresh = 0)

saveRDS(stan_model_first_stage, file = '05_stan_code/brms_first_stage.rds')

# first stage -----
stan_model_first_stage <- brm(X ~ Z, data = data, chains = 4, iter = 2000, warmup = 1000, 
                              save_model = "05_stan_code/brms_first_stage.stan",  silent = 2, refresh = 0)

saveRDS(stan_model_first_stage, file = '05_stan_code/brms_first_stage.rds')

# second stage -----
# Second stage: regress Y on X_hat
stan_model_second_stage <- brm(Y ~ X_hat, data = data, chains = 4, iter = 2000, warmup = 1000,
                               save_model = "05_stan_code/brms_second_stage.stan", silent = 2, refresh = 0)

saveRDS(stan_model_second_stage, file = '05_stan_code/brms_second_stage.rds')


# run function, function in iv_brms_func.R
parameters <- brms_iv_function(data, stan_model_first_stage, stan_model_second_stage)
