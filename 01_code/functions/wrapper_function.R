wrapper_function <- function(path_in, row_num, ...){
  data <- readRDS(path_in)
  # function call y, w, z, x, tau_true, w1, w0, 
  results <- own_bcf_iv(data$y, data$w, data$z, data$X, data$tau_true, data$w1, data$w0)
  
#  print(paste('Dataset', row_num, 'completed.'))
  
  base::return(results)
}
