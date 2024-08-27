# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
       source))

# loop over parameters
effect_size <- 2 # effect size as a function?
compliance <- 0.75
uncorrelated <- FALSE
confounded <- TRUE
path <- 'C:\\Users\\Jens Klenke\\sciebo - Klenke, Jens (snjeklen@uni-duisburg-essen.de)@uni-duisburg-essen.sciebo.de\\BART_IV-data\\00_sim_data'

# confounded missing

# running the function 
tictoc::tic()
wrapper_data_generation(p_vec = c(10, 50, 100), covariates = 'cont-cov', uncorrelated = uncorrelated, 
                        effect_size = effect_size, baseline_effect = TRUE, compliance = compliance, 
                        confounded = confounded)
tictoc::toc()
