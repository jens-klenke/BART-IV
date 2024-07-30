brms_iv_function <- function(data, model_first_stage, model_second_stage, ...){
  # first stage 
  first_stage <- update(model_first_stage, newdata = data)
  
  # Extract fitted values (predicted X)
  data$w_hat <- fitted(first_stage)[,1]
  
  # second stage 
  second_stage <- update(model_second_stage, newdata = data)
  
  # return parameters
  tau_hat <- summary(second_stage)$fixed
  
  return(tau_hat)
}
