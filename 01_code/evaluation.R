# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
                 source))

# parallel plan
future::plan(multisession, workers = 10)
options(future.globals.maxSize = 2147483648) # 2GB  

# load raw estimations
estimation_results <- readRDS(here::here('03_sim_results/discrete_covariates.rds'))

names(estimation_results)
# unnest estimation results

result <- tibble::as_tibble(estimation_results$results[1], .name_repair = 'minimal') %>%
  dplyr::rename('results' = 1) %>%
  dplyr::mutate(names = names(.$results)) %>%
  tidyr::pivot_wider(names_from = names, values_from = results)

get_sub_group_results <- function(results, ...){
  results %>%
    .$bcf_results %>%
    dplyr::filter(node == 'x1>=0.5 & x2>=0.5')
}

get_bcf_results <- function(results, ...){
  results %>%
    .$bcf_results
}

try <- estimation_results %>%
  dplyr::slice_head(n = 1000)

try_1 <- try %>%
  dplyr::mutate(bcf_result = purrr::pmap(., get_bcf_results, .progress = TRUE)) 


try_1 %<>%
  dplyr::select(ncov, row_num, bcf_result)



result %>%
  dplyr::select(bcf_results) %>%
  tidyr::unnest(bcf_results) %>%
  dplyr::filter(node == 'x1>=0.5 & x2>=0.5') %>%
  dplyr::select(node, CCACE)
