# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
       source))

# load stan models ----
#stan_model_first_stage <- readRDS(here::here("05_stan_code/brms_first_stage.rds"))
#stan_model_second_stage <- readRDS(here::here("05_stan_code/brms_second_stage.rds"))

# parallel plan
future::plan(multisession, workers = 10)
options(future.globals.maxSize = 2147483648) # 2GB  

# reading data 
data <- tibble::tibble(
  # path for loading data
  path_in = list.files('C:\\Users\\Jens Klenke\\Documents\\BART-IV\\00_sim_data', # sim_data_path(),  
    recursive = TRUE, full.names = TRUE)
  ) %>%
#  dplyr::filter(str_detect(path_in, 'effect.1')) %>%
  dplyr::mutate(ncov = readr::parse_number(stringr::str_extract(path_in, pattern = 'ncov.[0-9]*'), 
                                           locale =  readr::locale(decimal_mark = ",")),
                row_num = paste(dplyr::row_number(), 'of', max(dplyr::row_number())))

try <- data %>%
  sample_n(1)
  
# 
tictoc::tic()
sim_results_pmap <- try %>%
  dplyr::mutate(results = purrr::pmap(., wrapper_function, .progress = TRUE))
tictoc::toc()

tictoc::tic()
sim_results <- try %>%
  dplyr::mutate(results = furrr::future_pmap(., wrapper_function, .progress = TRUE, .options = furrr_options(seed = T)))
tictoc::toc()

# tictoc::tic()
# sim_results <- data %>%
#   dplyr::mutate(results = purrr::pmap(., wrapper_function, .progress = TRUE))
# tictoc::toc()

# saveRDS(sim_results, here::here('03_sim_results/discrete_covariates_cost_func_n.rds'))


try <- data %>%
  dplyr::mutate(effect = readr::parse_number(stringr::str_extract(path_in, "effect.\\d+"),
                                             locale =  readr::locale(decimal_mark = ",")),
                compliance = ifelse(str_detect(path_in, '0.75'), 0.75, 0.5),
                uncorrelated = ifelse(str_detect(path_in, 'uncorrelated'), 'uncorrelated','correlated'),
                confounded = ifelse(str_detect(path_in, 'baseline.efconfounded'), 'confounded', 'unconfounded')
  ) %>%
  dplyr::group_by(effect, compliance, uncorrelated, confounded, ncov) %>%
  count()

