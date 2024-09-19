# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
                 source))

# load example data
load(here::here('03_sim_results/effect.1/ef.1_co.0.75_baseline.ef_correlated_confounded.RData'))

ef.1_co.0.75_cor_conf <- sim_results

ef.1_co.0.75_cor_conf %<>%
  dplyr::mutate(dataset_num = as.numeric(str_extract(row_num, "^[^ of]+"))) %>%
  dplyr::mutate(leave_results = purrr::pmap(., analysis_fun)) %>%
  dplyr::rename('raw_data' = results) %>%
  dplyr::select(!path_in) %>%
  tidyr::unnest(leave_results)

pryr::object_size(ef.1_co.0.75_cor_conf)

# old ----
# per row (with map)
results_per_n <- sim_results$results[[1]]


try <- sim_results %>%
  dplyr::select(results)




try_analysis <-try %>%
  dplyr::mutate(purrr::pmap_df(., metric_fun_df))

row_try <- try %>%
  dplyr::slice(4)

results <- row_try %>% dplyr::select(results) %>% tidyr::unnest(results)






# inference and bayes
leaves_bcfiv_pred <- bcf_leave_results$bcfivResults %>%
  dplyr::filter(leaves == 1) %>%
  dplyr::select(pred, node) %>%
  tidyr::unnest(pred)

bcf_leave_results$bcfivResults %>%
  dplyr::filter(est_problems == 'yes') %>%
  dplyr::select(pred) %>%
  tidyr::unnest(pred)
  

leave_result <- tibble::tibble(
  # Model
  model = 'bcf_inference', 
  # Pehe leaves
  pehe_leaves = PEHE_fun(leaves_bcfiv_pred$tau_true, leaves_bcfiv_pred$tau_pred),
  # bias leaves
  bias_leaves = bias_fun(leaves_bcfiv_pred$tau_true, leaves_bcfiv_pred$tau_pred),
  # abs bias leaves
  abs_bias_leaves = abs_bias_fun(leaves_bcfiv_pred$tau_true, leaves_bcfiv_pred$tau_pred),
  # coverage leaves
  coverage_leaves = mean(leaves_bcfiv_pred$coverage)
)

bcf_leave_results.bayes <- bcf_leave_results$bayes_ivResults


# infernece
leaves_bcfiv_pred <- s_leave_results$bcfivResults %>%
  dplyr::filter(leaves == 1) %>%
  dplyr::select(pred) %>%
  tidyr::unnest(pred)

levae_result <- tibble::tibble(
  # Pehe leaves
  pehe_leaves = PEHE_fun(leaves_bcfiv_pred$tau_true, leaves_bcfiv_pred$tau_pred),
  # bias leaves
  bias_leaves = bias_fun(leaves_bcfiv_pred$tau_true, leaves_bcfiv_pred$tau_pred),
  # abs bias leaves
  abs_bias_leaves = abs_bias_fun(leaves_bcfiv_pred$tau_true, leaves_bcfiv_pred$tau_pred),
  # coverage leaves
  coverage_leaves = mean(leaves_bcfiv_pred$coverage)
)
