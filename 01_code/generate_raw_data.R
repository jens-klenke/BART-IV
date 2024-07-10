# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
       source))

# loop over parameters
p <- c(10, 50, 100)
covariates <- 'cont-cov'

tictoc::tic()

for (i in seq_along(p)) {
  p_i <- p[i]
  for (j in 1:100) {
    generate_dataset(n = 1000, p = p_i, covariates = covariates, base_line_effect = T) %>%
    saveRDS(file =here::here(
      paste0('00_sim_data/baseline_try_', p_i,'/dataset_ncovs', p_i, '_', j)
                )
    )
    if(j %% 10 == 0)
      cat( j, "100 Dataset of n covariates", p_i, "finished. \n")
  }
}

tictoc::toc()

dataset <- generate_dataset(n = 10, p = 10, base_line_effect = T)
