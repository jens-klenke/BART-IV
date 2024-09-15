# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
       source))

# load stan models ----
stan_model_first_stage <- readRDS(here::here("05_stan_code/brms_first_stage_4000_mac.rds"))
stan_model_second_stage <- readRDS(here::here("05_stan_code/brms_second_stage_4000_mac.rds"))

# parallel plan
future::plan(multisession, workers = 5)
#options(future.globals.maxSize = 2147483648) # 2GB  

# reading data 
data <- tibble::tibble(
  # path for loading data
  path_in = list.files(sim_data_path(), # 'C:\\Users\\Jens Klenke\\Documents\\BART-IV\\00_sim_data', #   
    recursive = TRUE, full.names = TRUE)
  ) %>%
  dplyr::filter(str_detect(path_in, 'ef.1_co.0.75_baseline.ef_correlated_confounded')) %>%
  dplyr::mutate(ncov = readr::parse_number(stringr::str_extract(path_in, pattern = 'ncov.[0-9]*'), 
                                           locale =  readr::locale(decimal_mark = ",")))

data_10 <- data %>%
  dplyr::filter(ncov == 10) %>%
  dplyr::mutate(row_num = paste(dplyr::row_number(), 'of', max(dplyr::row_number())))

data_50 <- data %>%
  dplyr::filter(ncov == 50) %>%
  dplyr::mutate(row_num = paste(dplyr::row_number(), 'of', max(dplyr::row_number())))

data_100 <- data %>%
  dplyr::filter(ncov == 100) %>%
  dplyr::mutate(row_num = paste(dplyr::row_number(), 'of', max(dplyr::row_number())))

#### estimation ----
# ncov = 10
tictoc::tic()
sim_10_results <- data_10 %>%
  dplyr::mutate(results = furrr::future_pmap(., wrapper_function, 
                                             .progress = TRUE, .options = furrr_options(seed = T)))
tictoc::toc()

# ncov = 50
tictoc::tic()
sim_50_results <- data_50 %>%
  dplyr::mutate(results = furrr::future_pmap(., wrapper_function, 
                                             .progress = TRUE, .options = furrr_options(seed = T)))
tictoc::toc()

# ncov = 100
tictoc::tic()
sim_100_results <- data_100 %>%
  dplyr::mutate(results = furrr::future_pmap(., wrapper_function, 
                                             .progress = TRUE, .options = furrr_options(seed = T)))
tictoc::toc()

#### combine -----
sim_results <- dplyr::bind_rows(
  sim_10_results,
  sim_50_results,
  sim_100_results
)

# clean up
unlink(tempdir(), recursive = TRUE, force = TRUE)

save(sim_results, file = here::here('03_sim_results/effect.1/ef.1_co.0.75_baseline.ef_correlated_confounded.RData'))

# try 
try <- data %>%
  dplyr::sample_n(2)

tictoc::tic()
sim_results <- try %>%
  dplyr::mutate(results = purrr::pmap(., wrapper_function, .progress = TRUE))
tictoc::toc()
