# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))

# source all files in the functions folder
sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
       source)


# loop over parameters
p <- c(10, 50, 100)

tictoc::tic()

for (i in seq_along(p)) {
  p_i <- p[i]
  for (j in 1:100) {
    generate_dataset(n = 1000, p = p_i) %>%
    saveRDS(file =here::here(
      paste0('00_sim_data/ncovs_', p_i,
             '/dataset_ncovs', p_i, '_', j)
                )
    )
  }
  cat("100 Dataset of n covariates", p_i, "finished. \n")
}

tictoc::toc()

dataset <- generate_dataset(n = 100, p = 10, continuous_covariates = TRUE)
