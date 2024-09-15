# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
                 source))

# load example data 
load(here::here('est_problems.RData'))

# 
sim_results <- try_sim_results
##
# per row (with map)
results_per_n <- sim_results$results[[1]]

# initialize Dataframe
df <- tibble::tibble(
  'detect_model' = rep(c('bcf', 's_bcf'), times = c(2, 2)) ,
  'est_method' = rep(c('inferential', 'bayes'), times = 2)
)



# bcf vs sbcf
bcf_leave_results <- results_per_n$bcf_results
s_bcf_leave_results <- results_per_n$s_bcf_results

inference_results <- bcf_leave_results$bcfivResults
bayes_results <- bcf_leave_results$bayes_ivResults

# check if we FALSE -> estimation problems!
check_inference <- ifelse(sum(inference_results$est_problems == 'yes') > 0, 
                          FALSE, TRUE)
check_bayes <- ifelse(sum(bayes_results$est_problems == 'yes') > 0, 
                          FALSE, TRUE)

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
