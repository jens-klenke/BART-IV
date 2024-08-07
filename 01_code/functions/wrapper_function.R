wrapper_function <- function(path_in, row_num, stan_model_first_stage, 
                             stan_model_second_stage, ...){
  data <- readRDS(path_in)
  
  results <- own_bcf_iv(data$y, data$w, data$z, data$X)
  
  print(paste('Dataset', row_num, 'completed.'))
  
  base::return(results)
}
