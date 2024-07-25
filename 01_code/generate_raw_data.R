# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
       source))

# loop over parameters
p <- c(10, 50, 100)
covariates <- 'cont-cov'
uncorrelated <- TRUE
effect_size <- 1

tictoc::tic()

for (i in seq_along(p)) {
  p_i <- p[i]
  for (j in 1:100) {
    generate_dataset(n = 1000, p = p_i, covariates = covariates, base_line_effect = T, uncorrelated = uncorrelated,
                     effect_size = effect_size) %>%
    saveRDS(file =here::here(
      paste0('00_sim_data/effect_1/baseline_uncorrelated_', p_i,'/dataset_ncovs_', p_i, '_uncorrelated_', j)
                )
    )
    if(j %% 10 == 0)
      cat( j, " of 100 Dataset with", p_i, "covariates finished. \n")
  }
}

tictoc::toc()

#dataset <- generate_dataset(n = 10, p = 10, base_line_effect = T)
