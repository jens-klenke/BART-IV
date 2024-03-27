# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
                 source))

## load data
# reading data 
data <- readRDS(here::here('03_sim_results/baseline_mu.rds'))

df <- data %>% 
  dplyr::slice_head(n = 1) %>%
  tidyr::unnest(results) %>%
  # # name results 
  dplyr::mutate(result_names = names(results)) %>%
  # parse ID 
  dplyr::mutate(ID = readr::parse_number(
    str_extract(path_in, "_(?:.(?!_))+$"))) %>%
  # delselect unimportant variables
  dplyr::select(-c(row_num, path_in))
  

names(df$results)
  
  
  
  
  
df <- cbind.data.frame(dataset$y, dataset$X)

colnames(df) <- c('y', paste0('x', 1:ncol(dataset$X)))

# mean of y
df %>%
  dplyr::summarise(y_bar = mean(y), y_sd = sd(y))

df %>%
  dplyr::mutate(subgroups = dplyr::case_when(
    x1 == 0 & x2 == 0 ~ 'positive_effect',
    x1 == 1 & x2 == 1 ~ 'negative_effect',
    .default = NA
  )) %>%
  dplyr::filter(!is.na(subgroups)) %>%
  dplyr::group_by(subgroups) %>%
  dplyr::summarise(y_bar = mean(y), y_sd = sd(y))

