# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
                 source))

# parallel plan
#future::plan(multisession, workers = 50) # base::floor(parallel::detectCores()*0.75))
#options(future.globals.maxSize = 2147483648) # 2GB  

# reading data -----
data <- tibble::tibble(
  # path for loading data
  path_in = list.files(sim_data_path(), recursive = TRUE, full.names = TRUE)
) %>%
  dplyr::filter(str_detect(path_in, 'ef.2_co.0.75_baseline.ef_correlated__')) %>%
  dplyr::mutate(ncov = readr::parse_number(stringr::str_extract(path_in, pattern = 'ncov.[0-9]*'), 
                                           locale =  readr::locale(decimal_mark = ",")),
                row_num = paste(dplyr::row_number(), 'of', max(dplyr::row_number())))

try <- data %>%
  sample_n(1)


tictoc::tic()
sim_results_pmap <- try %>%
  dplyr::mutate(results = purrr::pmap(., wrapper_function, .progress = TRUE))
tictoc::toc()


save(sim_results, file =' ')
# save(sim_results, file = here::here('03_sim_results/effect.1/ef.1_co.0.5_baseline.ef_correlated.RData'))

quit()


