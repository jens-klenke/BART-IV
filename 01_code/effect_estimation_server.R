# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
                 source))

# parallel plan
future::plan(multisession, workers = 40) # base::floor(parallel::detectCores()*0.75))
options(future.globals.maxSize = 2147483648) # 2GB  

# reading data -----
data <- tibble::tibble(
  # path for loading data
  path_in = list.files(
    here::here("00_sim_data/effect.2"), recursive = TRUE, full.names = TRUE)
) %>%
  #  dplyr::filter(str_detect(path_in, 'try')) %>%
  dplyr::mutate(ncov = readr::parse_number(stringr::str_extract(path_in, pattern = 'ncov.[0-9]*'), 
                                           locale =  readr::locale(decimal_mark = ",")),
                row_num = paste(dplyr::row_number(), 'of', max(dplyr::row_number())))



#tictoc::tic()
#sim_results <- data %>%
#  dplyr::mutate(results = purrr::pmap(., wrapper_function, .progress = TRUE))
#tictoc::toc()

tictoc::tic()
sim_results <- data %>%
  dplyr::mutate(results = furrr::future_pmap(., wrapper_function, .progress = TRUE, .options = furrr_options(seed = T)))
tictoc::toc()

save(sim_results, file ='Z:/Data/bayes_iv-try.RData')

quit()

