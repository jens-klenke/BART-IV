#'
#' 
wrapper_data_generation <- function(H = 100, n = 1000, p_vec, covariates, uncorrelated, effect_size, baseline_effect, compliance = 0.75, confounded) {
  # function head 
  
  baseline <- ifelse(baseline_effect, 'baseline.ef', 'no.baseline.ef')
  uncorrelated_data <- ifelse(uncorrelated, 'uncorrelated', 'correlated')
  conf  <- ifelse(confounded, 'confounded', '')
  
  for (i in seq_along(p_vec)) {
    # get number of covvariates  
    p <- p_vec[i]  
    # folder path 
    folder_path <- here::here('00_sim_data', 
                              paste0('effect.', effect_size),
                              paste0('compliance.', compliance),
                              uncorrelated_data,
                              paste0(baseline, conf),
                              paste0('ncov.', p))
    
    # name of the dataset
    data_name <- paste0('ef.', effect_size, '_',
                        'co.', compliance, '_',
                        baseline, '_',
                        uncorrelated_data, '_',
                        conf, '_',
                        'ncov.', p)
    
    if(!dir.exists(folder_path)){
      dir.create(folder_path, recursive = TRUE)
    } 
    
    ### generate data
    for (j in 1:H){
      generate_dataset(n = n, p = p, covariates = covariates, base_line_effect = baseline_effect, uncorrelated = uncorrelated,
                       effect_size = effect_size, confounded = confounded) %>%
        saveRDS(file = paste0(folder_path, '/', data_name, '_', j))
      # printing
      if(j %% 10 == 0)
        cat( j, " of 100 Dataset with", p, "covariates finished. \n")
    }
    
  }
}
