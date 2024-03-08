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
estimation_results_raw <- readRDS(here::here('03_sim_results/discrete_covariates.rds'))
estimation_results_raw <- readRDS(here::here('03_sim_results/discrete_covariates_cost_func.rds'))


##############################################
#####              wrangling             #####
##############################################

estimation_results <- estimation_results_raw %>%
  # unnest results
  tidyr::unnest(results) %>%
  # # name results 
  dplyr::mutate(result_names = names(results)) %>%
  # parse ID 
  dplyr::mutate(ID = readr::parse_number(str_extract(path_in, "_(?:.(?!_))+$"))) %>%
  # delselect unimportant variables
  dplyr::select(-c(row_num, path_in)) %>%
  #  filter for main results
  dplyr::filter(result_names %in% c('bcf_results', 's_bcf_results')) %>%
  # force same datatype
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, node = as.character(node)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, CCACE = as.double(CCACE)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, pvalue = as.double(pvalue)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, Weak_IV_test = as.double(Weak_IV_test)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, Pi_obs = as.double(Pi_obs)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, ITT = as.double(ITT)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, Pi_compliers = as.double(Pi_compliers)))) %>%
  # unnest bcf and bcfs results
  tidyr::unnest(results) %>%
  dplyr::arrange(ncov, ID) %>%
  dplyr::mutate(subgroup = dplyr::case_when(
    node == 'x1>=0.5 & x2>=0.5' | node == 'x2>=0.5 & x1>=0.5' |
    node == 'x2> 0.5 & x1> 0.5' | node == 'x1> 0.5 & x2> 0.5' ~ 'negative effect',
    node == 'x2< 0.5 & x1< 0.5' | node == 'x1< 0.5 & x2< 0.5' |
    node == 'x2<=0.5 & x1<=0.5' | node == 'x1<=0.5 & x2<=0.5' ~ 'positive effect',
    .default = NA)
    )

##############################################
#####    investigating the subgroups     #####
##############################################


estimation_results_ncov_10 <- estimation_results %>%
  dplyr::filter(ncov == 10) %>%
  dplyr::group_by(result_names, subgroup)

tibble::tibble(ID = 1:100) %>%
  dplyr::left_join(
    estimation_results_ncov_10 %>%
      dplyr::select(subgroup, CCACE, result_names, ID), 
    by = 'ID') %>%
  dplyr::filter(subgroup %in% c('negative effect', 'positive effect')) %>%
  ggplot2::ggplot(aes(x = ID, y = CCACE, color = result_names)) +
  ggplot2::geom_point(size = 2, alpha = 0.95) +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap('subgroup')+
  ggplot2::geom_smooth(method = lm, formula = y ~ 1, se = FALSE)

estimation_results_ncov_50 <- estimation_results %>%
  dplyr::filter(ncov == 50) %>%
  dplyr::group_by(result_names, subgroup)

tibble::tibble(ID = 1:100) %>%
  dplyr::left_join(
    estimation_results_ncov_50 %>%
      dplyr::select(subgroup, CCACE, result_names, ID), 
    by = 'ID') %>%
  dplyr::filter(subgroup %in% c('negative effect', 'positive effect')) %>%
  ggplot2::ggplot(aes(x = ID, y = CCACE, color = result_names)) +
  ggplot2::geom_point(size = 2, alpha = 0.95) +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap('subgroup')+
  ggplot2::geom_smooth(method = lm, formula = y ~ 1, se = FALSE)

estimation_results_ncov_100 <- estimation_results %>%
  dplyr::filter(ncov == 100) %>%
  dplyr::group_by(result_names, subgroup)

tibble::tibble(ID = 1:100) %>%
  dplyr::left_join(
    estimation_results_ncov_100 %>%
      dplyr::select(subgroup, CCACE, result_names, ID), 
    by = 'ID') %>%
  dplyr::filter(subgroup %in% c('negative effect', 'positive effect')) %>%
  ggplot2::ggplot(aes(x = ID, y = CCACE, color = result_names)) +
  ggplot2::geom_point(size = 2, alpha = 0.95) +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap('subgroup')+
  ggplot2::geom_smooth(method = lm, formula = y ~ 1, se = FALSE)

# table
estimation_results_ncov_10 %>%
  dplyr::filter(subgroup %in% c('negative effect', 'positive effect')) %>%
  dplyr::summarise(n = dplyr::n(), share = dplyr::n()/100, mean = mean(CCACE), sd = sd(CCACE)) %>%
  tidyr::pivot_wider(names_from = result_names, values_from = c(n, share, mean, sd)) %>%
  dplyr::relocate(subgroup, n_bcf_results, share_bcf_results, mean_bcf_results, sd_bcf_results)

try <- estimation_results_ncov_10 %>%
  dplyr::filter(subgroup %in% c('negative effect', 'positive effect')) %>%
  dplyr::group_by(result_names, ID) %>%
  dplyr::count() %>%
  dplyr::filter(n == 2) %>%
  dplyr::group_by(result_names) %>%
  dplyr::summarise(n = dplyr::n())



try


#### --------------- trash ----------------







####

###
try_1 <- estimation_results %>%
  # unnest results
  tidyr::unnest(results) %>%
  # # name results 
  dplyr::mutate(result_names = names(results)) %>%
  # parse ID 
  dplyr::mutate(ID = readr::parse_number(str_extract(path_in, "_(?:.(?!_))+$"))) %>%
  # delselect unimportant variables
  dplyr::select(-c(row_num, path_in)) %>%
  #  filter for main results
  dplyr::filter(result_names %in% c('bcf_results', 's_bcf_results')) %>%
  # force same datatype
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, node = as.character(node)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, CCACE = as.double(CCACE)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, pvalue = as.double(pvalue)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, Weak_IV_test = as.double(Weak_IV_test)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, Pi_obs = as.double(Pi_obs)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, ITT = as.double(ITT)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, Pi_compliers = as.double(Pi_compliers)))) %>%
  # unnest bcf and bcfs results
  tidyr::unnest(results) %>%
  dplyr::arrange(ncov, ID) %>%
  dplyr::filter(CCACE > 0)




  dplyr::select(bcf_result) %>%
  tidyr::unnest(bcf_result) %>%
  dplyr::filter(node == 'x1>=0.5 & x2>=0.5') %>%
  dplyr::select(node, CCACE)



###
result %>%
  dplyr::select(bcf_results) %>%
  tidyr::unnest(bcf_results) %>%
  dplyr::filter(node == 'x1>=0.5 & x2>=0.5') %>%
  dplyr::select(node, CCACE)















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

get_sbcf_results <- function(results, ...){
  results %>%
    .$s_bcf_results
}

try <- estimation_results %>%
  dplyr::slice_head(n = 1)

try_1 <- try %>%
  dplyr::mutate(bcf_result = purrr::pmap(., get_bcf_results, .progress = TRUE),
                s_bcf_result = purrr::pmap(., get_sbcf_results, .progress = TRUE)) 


try_1 %<>%
  dplyr::select(ncov, row_num, bcf_result)

###
try_1 %>%
  dplyr::select(bcf_result) %>%
  tidyr::unnest(bcf_result) %>%
  dplyr::filter(node == 'x1>=0.5 & x2>=0.5') %>%
  dplyr::select(node, CCACE)



###
result %>%
  dplyr::select(bcf_results) %>%
  tidyr::unnest(bcf_results) %>%
  dplyr::filter(node == 'x1>=0.5 & x2>=0.5') %>%
  dplyr::select(node, CCACE)



  
