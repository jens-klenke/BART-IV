both_detected <- function(df, ...){
  df %>%
    dplyr::filter(subgroup %in% c('negative effect', 'positive effect')) %>%
    dplyr::group_by(result_names, ID) %>%
    dplyr::add_count() %>%
    dplyr::filter(n == 2) %>%
    dplyr::group_by(result_names, subgroup) %>%
    dplyr::summarise(both = dplyr::n())
}

analyse_df <- function(data, ...){
  data %>%
    dplyr::filter(subgroup %in% c('negative effect', 'positive effect')) %>%
    dplyr::summarise(n = dplyr::n(), share = dplyr::n()/100,
                     mean = mean(CCACE), sd = sd(CCACE)) %>%
    dplyr::left_join(
      data %>%
        both_detected(),
      by = c('result_names', 'subgroup')
    ) %>%
    tidyr::pivot_wider(names_from = result_names,
                       values_from = c(n, share, mean, sd, both)) %>%
    dplyr::relocate(subgroup, n_bcf_results, share_bcf_results,
                    mean_bcf_results, sd_bcf_results, both_bcf_results)
}

## get counts for the plot 
get_counts <- function(data, group.by = c('ncov', 'subgroup', 'result_names')){
  # Convert character vector to symbols
  vars_group <- syms(group.by)
  
  data %<>% 
    # make grid
    tidyr::expand(!!!vars_group) %>%
    dplyr::left_join(
      data %>%
        dplyr::group_by(!!!vars_group) %>%
        dplyr::summarise(n = dplyr::n(), .groups = 'drop'),
      by = group.by
    ) %>%
    tidyr::replace_na(list(n = 0)) %>%
    dplyr::arrange(!!!vars_group)
  
  # return data
  return(data)
  
}

# get vlines in effects
get_vline <- function(data , group.by = c('ncov', 'subgroup', 'result_names')){
  # Convert character vector to symbols
  vars_group <- syms(group.by)
  
  data %<>%
    tidyr::expand(!!!vars_group) %>%
    dplyr::mutate(x_line = ifelse(subgroup == 'positive effect', 2, -2))
  
  return(data)
}