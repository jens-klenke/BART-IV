# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
       source))

# parallel plan
future::plan(multisession, workers = parallel::detectCores()*.9)
options(future.globals.maxSize = 2147483648) # 2GB  

# # reading data 
data <- tibble::tibble(
  # path for loading data
  path_in = list.files(
    here::here("00_sim_data/"), recursive = TRUE, full.names = TRUE)
  ) %>%
  dplyr::mutate(ncov = readr::parse_number(stringr::str_extract(path_in, pattern = 'ncovs_[0-9]*')),
                row_num = paste(dplyr::row_number(), 'of', max(dplyr::row_number())))

# example data set
try <- data %>%
  dplyr::slice_head(n = 2)

tictoc::tic()
try_1 <- try %>%
  dplyr::mutate(results = purrr::pmap(., wrapper_function, .progress = TRUE))
tictoc::toc()

# 
tictoc::tic()
parallel_try_1 <- try %>%
  dplyr::mutate(results = furrr::future_pmap(., wrapper_function, .progress = TRUE))
tictoc::toc()




# 1 and 1 works
i <- 1 
j <- 1
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
own_bcf <- own_bcf_iv(data$y, data$w, data$z, data$X)
tictoc::toc()

tictoc::tic()
bcf_package <- BayesIV::bcf_iv(data$y, data$w, data$z, data$X, n_burn = 1000)
tictoc::toc()

tictoc::tic()
first_bcf <- bcf_iv(data$y, data$w, data$z, data$X, n_burn = 1000)
tictoc::toc()

pryr::object_size(own_bcf)



