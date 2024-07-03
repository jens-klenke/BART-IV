# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
                 source))

## load data
### Data to big for Github -> download from 
# data <- readRDS(here::here('03_sim_results/baseline_mu.rds'))
data <- readRDS(url('https://uni-duisburg-essen.sciebo.de/s/Rt7sNPFwaVw1lr3/download'))

df <- data %>% 
  tidyr::unnest(results) %>%
  # # name results 
  dplyr::mutate(result_names = names(results)) %>%
  # parse ID 
  dplyr::mutate(ID = readr::parse_number(
    str_extract(path_in, "_(?:.(?!_))+$"))) %>%
  # delselect unimportant variables
  dplyr::select(-c(row_num, path_in)) %>%
  dplyr::filter(stringi::stri_detect(result_names, regex = 'bcf_exp'))


df_1 <- df %>%
  # unnest results to get the individual estimates 
  tidyr::unnest(results) %>%
  # code subgroups
  dplyr::mutate(subgroups = dplyr::case_when(
    x1 == 0 & x2 == 0 ~ 'positive_effect',
    x1 == 1 & x2 == 1 ~ 'negative_effect',
    .default = NA
  )) %>%
  # filter for subgroups
  dplyr::filter(!is.na(subgroups)) %>%
  dplyr::mutate(tauhat = dplyr::case_when(
    result_names == 'bcf_exp' ~ bcf_tauhat,
    result_names == 's_bcf_exp' ~ s_bcf_tauhat,
  )) %>%
  # group by method and dataset
  dplyr::group_by(result_names, subgroups, ID, ncov) %>%
  dplyr::summarise(mean = mean(tauhat), sd = sd(tauhat),
                   median = median(tauhat), min = min(tauhat),
                   max = max(tauhat))

df_1 %>%
  dplyr::filter(result_names == 'bcf_exp') %>%
  ggplot2::ggplot(aes(x = median)) +
  ggplot2::geom_density() +
  ggplot2::facet_grid(ncov ~ subgroups) +
  ggplot2::theme_bw()


df_1 %>%
  dplyr::filter(result_names == 's_bcf_exp') %>%
  ggplot2::ggplot(aes(x = median)) +
  ggplot2::geom_density() +
  ggplot2::facet_grid(ncov ~ subgroups) +
  ggplot2::theme_bw()


