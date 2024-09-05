brms_iv_function <- function(data, ...){
  # first stage 
  first_stage <- brms:::update.brmsfit(
    readRDS(here::here("05_stan_code/brms_first_stage.rds")),
    #model_first_stage, 
    newdata = data)
  
#  print(paste('inside brms function:', nrow(data)))
  
  # Extract fitted values (predicted X)
  data$w_hat <- fitted(first_stage)[,1]
  
  # second stage 
  second_stage <- brms:::update.brmsfit(
    readRDS(here::here("05_stan_code/brms_second_stage.rds")), #model_second_stage, 
    newdata = data)
  
  # return parameters
  tau_hat <- summary(second_stage)$fixed
  #  print(paste('effect:', summary(second_stage)$fixed[2, 1]))
  
  rm(first_stage, second_stage)
  
  return(tau_hat)
}
