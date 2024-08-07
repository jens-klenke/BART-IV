# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
       source))

# load stan models ----
stan_model_first_stage <- readRDS(here::here("05_stan_code/brms_first_stage.rds"))
stan_model_second_stage <- readRDS(here::here("05_stan_code/brms_second_stage.rds"))

# parallel plan
future::plan(multisession, workers = 10)
options(future.globals.maxSize = 2147483648) # 2GB  

# reading data 
data <- tibble::tibble(
  # path for loading data
  path_in = list.files(
    here::here("00_sim_data/effect.2"), recursive = TRUE, full.names = TRUE)
  ) %>%
#  dplyr::filter(str_detect(path_in, 'try')) %>%
  dplyr::mutate(ncov = readr::parse_number(stringr::str_extract(path_in, pattern = 'ncov.[0-9]*'), 
                                           locale =  readr::locale(decimal_mark = ",")),
                row_num = paste(dplyr::row_number(), 'of', max(dplyr::row_number()))) %>%
  dplyr::mutate(stan_model_first_stage = list(stan_model_first_stage),
                stan_model_second_stage = list(stan_model_second_stage))

try <- data %>%
  sample_n(1)
  
# 

tictoc::tic()
sim_results <- try %>%
  dplyr::mutate(results = purrr::pmap(., wrapper_function, .progress = TRUE))
tictoc::toc()

tictoc::tic()
sim_results <- try %>%
  dplyr::mutate(results = furrr::future_pmap(., wrapper_function, .progress = TRUE))
tictoc::toc()

# tictoc::tic()
# sim_results <- data %>%
#   dplyr::mutate(results = purrr::pmap(., wrapper_function, .progress = TRUE))
# tictoc::toc()

saveRDS(sim_results, here::here('03_sim_results/discrete_covariates_cost_func_n.rds'))

################################################################################
#####                          continuous data                              #####
################################################################################
# reading data 
data <- tibble::tibble(
  # path for loading data
  path_in = list.files(
    here::here("00_sim_data/"), recursive = TRUE, full.names = TRUE)
  ) %>%
  dplyr::filter(str_detect(path_in, 'cont-cov')) %>%
  dplyr::mutate(ncov = readr::parse_number(stringr::str_extract(path_in, pattern = 'ncovs_[0-9]*')),
                row_num = paste(dplyr::row_number(), 'of', max(dplyr::row_number())))

# 
tictoc::tic()
sim_results <- data %>%
  dplyr::mutate(results = furrr::future_pmap(., wrapper_function, .progress = TRUE))
tictoc::toc()

################################################################################
#####                                 save                                #####
################################################################################
saveRDS(sim_results, here::here('03_sim_results/continuous_covariates.rds'))

# 
tictoc::tic()
sim_results <- data %>%
  dplyr::mutate(results = purrr::pmap(., wrapper_function, .progress = TRUE))
tictoc::toc()

################################################################################
#####                                 save                                #####
################################################################################
saveRDS(sim_results, here::here('03_sim_results/continuous_covariates.rds'))


################################################################################
#####                          baseline mu data                            #####
################################################################################
# reading data 
data <- tibble::tibble(
  # path for loading data
  path_in = list.files(
    here::here("00_sim_data/"), recursive = TRUE, full.names = TRUE)
) %>%
  dplyr::filter(str_detect(path_in, 'baselinemu')) %>%
  dplyr::mutate(ncov = readr::parse_number(stringr::str_extract(path_in, pattern = 'ncovs_[0-9]*')),
                row_num = paste(dplyr::row_number(), 'of', max(dplyr::row_number())))

# 
tictoc::tic()
sim_results <- data %>%
  dplyr::mutate(results = furrr::future_pmap(., wrapper_function, .progress = TRUE))
tictoc::toc()

################################################################################
#####                                 save                                #####
################################################################################
saveRDS(sim_results, here::here('03_sim_results/baseline_mu.rds'))
saveRDS(sim_results, 'C:/Users/jens.klenke/Dropbox/jens/sim_data/baseline_mu.rds')










################################################################################
#####                                 save                                #####
################################################################################
saveRDS(sim_results, here::here('03_sim_results/discrete_covariates.rds'))

################################################################################
#####                                 tries                                #####
################################################################################






