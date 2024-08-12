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
baseline_effect <- TRUE
effect_size <- 1 # effect size as a function?
compliance <- 0.75

# confounded missing

# running the function 
tictoc::tic()
wrapper_data_generation(p_vec = p, covariates = covariates, uncorrelated = uncorrelated, effect_size = effect_size, 
                        baseline_effect = TRUE, compliance = compliance)
tictoc::toc()
