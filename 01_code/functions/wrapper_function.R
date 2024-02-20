wrapper_function <- function(path_in, ...){
  data <- readRDS(path_in)
  
  results <- own_bcf_iv(data$y, data$w, data$z, data$X)
  
  return(results)
  
}
