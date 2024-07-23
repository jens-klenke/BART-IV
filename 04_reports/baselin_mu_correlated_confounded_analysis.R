# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
                 source))

## load data
### Data to big for Github -> download from -> Sciebo/BART-IV-Data 
load(url('https://uni-duisburg-essen.sciebo.de/s/1kkmrWzX38EL7ID/download'))

df <- sim_results %>% 
  tidyr::unnest(results) %>%
  # # name results 
  dplyr::mutate(result_names = names(results)) %>%
  # parse ID 
  dplyr::mutate(ID = readr::parse_number(
    str_extract(path_in, "_(?:.(?!_))+$"))) %>%
  # delselect unimportant variables
  dplyr::select(-c(row_num, path_in)) 

theoretical_df <- df %>%
  dplyr::filter(stringi::stri_detect(result_names, regex = 'theretical_results')) %>%
  # unnest results to get the individual estimates 
  tidyr::unnest(results) %>%
  dplyr::mutate(CCACE = as.numeric(CCACE)) %>%
  # group by method and dataset
  dplyr::group_by(node, ncov) 

# Theoretical results (heterogenous groups hard coded not detected via binary tree)
theoretical_df %>%
  dplyr::summarise(mean = mean(CCACE), sd = sd(CCACE),
                   median = median(CCACE), min = min(CCACE),
                   max = max(CCACE))

theoretical_df %>%
  ggplot2::ggplot(aes(x = CCACE)) +
  ggplot2::geom_density() +
  ggplot2::facet_grid(ncov ~ node) +
  ggplot2::theme_bw()

### results with binary detection ----
binary_results <- df %>%
  #  filter for main results
  dplyr::filter(stringi::stri_detect(result_names, regex = 'bcf_results')) %>%
  # force character datatype
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, node = as.character(node)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, CCACE = as.character(CCACE)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, pvalue = as.character(pvalue)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, Weak_IV_test = as.character(Weak_IV_test)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, Pi_obs = as.character(Pi_obs)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, ITT = as.character(ITT)))) %>%
  dplyr::mutate(results = purrr::map(results, ~ dplyr::mutate(.x, Pi_compliers = as.character(Pi_compliers)))) %>%
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

binary_results %<>%
  dplyr::filter(subgroup %in% c('negative effect', 'positive effect')) %>%
  # force double datatype
  dplyr::mutate(CCACE = as.double(CCACE),
                pvalue = as.double(pvalue),
                Weak_IV_test = as.double(Weak_IV_test),
                Pi_obs = as.double(Pi_obs),
                ITT = as.double(ITT),
                Pi_compliers = as.double(Pi_compliers)
  )

binary_results %>%
  ggplot2::ggplot(aes(x = CCACE)) +
  ggplot2::geom_density() +
  ggplot2::facet_grid(ncov ~ subgroup + result_names) +
  ggplot2::theme_bw() +
  # Add number of observations as text
  geom_text(data = get_counts(binary_results), aes(x = Inf, y = Inf, label = paste("n =", n)),
            hjust = 1.1, vjust = 1.1, size = 3, color = '#004c93') +
  geom_vline(data = get_vline(binary_results), aes(xintercept = x_line), color = '#004c93')



# BCF poorly! (100 ncovs no detection with correlation) 
# SBCF better and almost unaffected by high number of covariates!
binary_results %>%
  dplyr::filter(subgroup %in% c('negative effect', 'positive effect')) %>%
  dplyr::group_by(subgroup, ncov, result_names) %>%
  dplyr::summarise(n = dplyr::n(), share = dplyr::n()/100,
                   median = median(CCACE), mean = mean(CCACE), sd = sd(CCACE)) %>% 
  dplyr::arrange(result_names, subgroup, ncov)
