# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
       source))

i <- 1 
j <- 10
p <- c(10, 50, 100)
p_i <- p[i]

# reading in the data
data <- readRDS(
  here::here(
    paste0('00_sim_data/ncovs_', p_i,
           '/dataset_ncovs', p_i, '_', j)
    )
)

# own 
tictoc::tic()
own_bcf <- own_bcf_iv(data$y, data$w, data$z, data$X, n_burn = 1000)
tictoc::toc()

tictoc::tic()
bcf_package <- BayesIV::bcf_iv(data$y, data$w, data$z, data$X, n_burn = 1000)
tictoc::toc()

tictoc::tic()
first_bcf <- bcf_iv(data$y, data$w, data$z, data$X, n_burn = 1000)
tictoc::toc()

pryr::object_size(own_bcf)



