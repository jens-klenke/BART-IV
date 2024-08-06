library(AER)
library(brms)

# old example data 
#set.seed(123456)  # For reproducibility
## Sample size
#n <- 1000
## Generate the instrument (Z)
#Z <- rnorm(n)
## Generate confounder (U)
#U <- rnorm(n)
## Generate the treatment (X), affected by Z and U
#X <- 0.5 * Z + 0.3 * U + rnorm(n)
## Generate the outcome (Y), affected by X and U
#Y <- 2 * X + 0.5 * U + rnorm(n)
## Combine into a data frame ## needs to be changed! 
#data <- data.frame(Y, X, Z, U)

# example dataset 
dataset <- readRDS(here::here('00_sim_data/effect_2/baselinemu_ncovs_10/dataset_ncovs10_1'))
y <- dataset$y
w <- dataset$w
z <- dataset$z
x <- dataset$X
iv.data <- as.data.frame(cbind(y, w, z, x))

#try <- ivreg(y ~ w | z,  
#     data = iv.data) 
#summary(try)

#first_stage <- lm(w ~ z, data = iv.data)
#w_hat <- fitted.values(first_stage)

#data$w_hat <- w_hat

#second_stage <- lm(y ~ w_hat, data = iv.data)

# first stage -----
# first stage: regress X on Z
stan_model_first_stage <- brm(w ~ z, data = iv.data, chains = 4, iter = 2000, warmup = 1000, 
                              save_model = "05_stan_code/brms_first_stage.stan",  silent = 2, refresh = 0)

saveRDS(stan_model_first_stage, file = '05_stan_code/brms_first_stage.rds')

# add fitted values to data.frame 
iv.data$w_hat <- fitted(stan_model_first_stage)[,1]

# second stage -----
# Second stage: regress Y on X_hat
stan_model_second_stage <- brm(y ~ w_hat, data = iv.data, chains = 4, iter = 2000, warmup = 1000,
                               save_model = "05_stan_code/brms_second_stage.stan", silent = 2, refresh = 0)

saveRDS(stan_model_second_stage, file = '05_stan_code/brms_second_stage.rds')


# try ----
stan_model_first_stage <- readRDS('05_stan_code/brms_first_stage.rds')
stan_model_second_stage <- readRDS('05_stan_code/brms_second_stage.rds')

# run function, function in iv_brms_func.R
tictoc::tic()
parameters <- brms_iv_function(iv.data, stan_model_first_stage, stan_model_second_stage)
tictoc::toc()