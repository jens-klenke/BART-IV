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

# extract and summarize the data for each MCMC (4 Cases: (bcf, s_bcf) * (inference, bayes))
ef.1_co.0.75_cor_conf %<>%
  dplyr::mutate(dataset_num = as.numeric(str_extract(row_num, "^[^ of]+"))) %>%
  dplyr::mutate(leave_results = purrr::pmap(., analysis_fun)) %>%
  dplyr::rename('raw_data' = results) %>%
  dplyr::select(!path_in) %>%
  tidyr::unnest(leave_results) %>%
  dplyr::mutate(detect_model = dplyr::recode(detect_model, 
                                      'bcf_results' = 'BCF',
                                      's_bcf_results' = 'SBCF' ), 
                methods = dplyr::recode(methods, 
                                        'bcfivResults' = 'frequentist',
                                        'bayes_ivResults' = 'bayesian' ))


# statistics
first_stats <- ef.1_co.0.75_cor_conf %>%
  dplyr::group_by(ncov, detect_model, methods) %>%
  dplyr::summarise(est_problems_count = sum(est_problems),
                   across(c(bias, bias_rm.na, abs_bias, abs_bias_rm.na, PEHE, PEHE_rm.na,coverage),
                          list(mean = ~mean(., na.rm = TRUE), 
                               sd = ~sd(., na.rm = TRUE)),
                          .names = "{.col}_{.fn}")
                   )
# pryr::object_size(ef.1_co.0.75_cor_conf)

# detect subgroups and stats for subgroups -----
detect_subgroups <- ef.1_co.0.75_cor_conf %>%
  dplyr::select(detect_model, methods, ncov, dataset_num, results) %>%
  tidyr::unnest(results) %>%
  # detect subgroups 
  dplyr::mutate(subgroup = dplyr::case_when(
    node == 'x1>=0.5 & x2>=0.5' | node == 'x2>=0.5 & x1>=0.5' | node == 'x2> 0.5 & x1> 0.5' | node == 'x1> 0.5 & x2> 0.5' ~ 'negative effect',
    node == 'x2< 0.5 & x1< 0.5' | node == 'x1< 0.5 & x2< 0.5' | node == 'x2<=0.5 & x1<=0.5' | node == 'x1<=0.5 & x2<=0.5' ~ 'positive effect',
    .default = NA)
    ) %>%
  dplyr::filter(subgroup %in% c('negative effect', 'positive effect'))

detect_subgroups_n.e <- detect_subgroups %>%
  dplyr::filter(subgroup == 'negative effect') 

detect_subgroups_n.e %>%
  ggplot2::ggplot(aes(x = CCACE)) +
  ggplot2::geom_density() +
  ggplot2::facet_grid(ncov ~ detect_model + methods) +
 # ggplot2::theme_bw() +
  # Add number of observations as text
  geom_text(data = get_counts(detect_subgroups_n.e), aes(x = Inf, y = Inf, label = paste("n =", n)),
            hjust = 1.1, vjust = 1.1, size = 3, color = '#004c93') +
  ggplot2::geom_vline(data = get_vline(detect_subgroups_n.e, 1), 
                      aes(xintercept = x_line), color = '#004c93', 
                      linewidth = 0.75, linetype = 'longdash') +
  own_theme()


detect_subgroups_p.e <-detect_subgroups %>%
  dplyr::filter(subgroup == 'positive effect')

detect_subgroups_p.e %>%
  ggplot2::ggplot(aes(x = CCACE)) +
  ggplot2::geom_density() +
  ggplot2::facet_grid(ncov ~ detect_model + methods) +
  # Add number of observations as text
  geom_text(data = get_counts(detect_subgroups_p.e), aes(x = Inf, y = Inf, label = paste("n =", n)),
            hjust = 1.1, vjust = 1.1, size = 3, color = '#004c93') +
  ggplot2::geom_vline(data = get_vline(detect_subgroups_p.e, 1), 
                      aes(xintercept = x_line), color = '#004c93', 
                      linewidth = 0.75, linetype = 'longdash') +
  own_theme()


# BCF poorly! SBCF better and almost unaffected by high number of covariates!
detect_subgroups %>%
  dplyr::group_by(subgroup, ncov, detect_model, methods) %>%
  dplyr::summarise(n = dplyr::n(), 
                   share = dplyr::n()/500,
                   median = median(CCACE), 
                   mean = mean(CCACE), 
                   sd = sd(CCACE)) %>% 
  dplyr::arrange(subgroup, ncov, detect_model, methods)



# try analysis  ----