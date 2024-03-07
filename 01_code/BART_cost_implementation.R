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

# # reading data 
data <- tibble::tibble(
  # path for loading data
  path_in = list.files(
    here::here("00_sim_data/"), recursive = TRUE, full.names = TRUE)
) %>%
  dplyr::mutate(ncov = readr::parse_number(stringr::str_extract(path_in, pattern = 'ncovs_[0-9]*')),
                row_num = paste(dplyr::row_number(), 'of', max(dplyr::row_number())))

# bcf_iv <- function(
    # Data 
dataset <- readRDS(data$path_in[1])

y <- dataset$y
w <- dataset$w
z <- dataset$z
x <- dataset$X

# Parameters
binary = FALSE
n_burn = 1000
n_sim = 1000
inference_ratio = 0.5
max_depth = 2
cp = 0.01
minsplit = 30
adj_method = "holm"
seed = 42

x_names <- paste0('x', 1:ncol(x))

######################################################
####         Step 0: Initialize the Data          ####
######################################################

    # Split data into Discovery and Inference
    set.seed(seed)
    index <- sample(nrow(x), nrow(x)*inference_ratio, replace=FALSE)

    # Initialize total dataset
    iv.data <- as.data.frame(cbind(y, w, z, x))
    names(iv.data) <- c("y", "w", "z", x_names)

    # Discovery and Inference Samples
    discovery <- iv.data[-index,]
    inference <- iv.data[index,]

######################################################
####  Step 1: Compute the Bayesian Causal Forest  ####
######################################################

  # Compute the Propensity Score though a Logistic Regression
  p.score <- glm(z ~ x[-index,],
               family = binomial,
               data = discovery)
  pihat <- predict(p.score, as.data.frame(x[-index,]))

  # Perform the Bayesian Causal Forest  to calculate the Proportion of Compliers (pic)
  bcf_pic_tree <- quiet(bartCause::bartc(w[-index], z[-index], x[-index,],
                                       n.samples = n_sim, n.burn = n_burn, 
                                       n.chains = 2L))

  bcf_tau_pic <- bartCause::extract(bcf_pic_tree, type = "ite")
  bcf_pic <- apply(bcf_tau_pic, 2, mean)
  # workaround! 
  bcf_pic[bcf_pic == 0] <- 1e-2 
  # mean(pic_bcf) / median(pic) == compliance

######################################################
####     Continuous and Discrete Outcomes         ####
######################################################

  # Perform the Bayesian Causal Forest for the ITT
  # BCF 
  tictoc::tic()
  bcf_itt.tree <- quiet(
    SparseBCF::SparseBCF(y[-index], z[-index], x[-index,], pihat = pihat,
                         nsim = n_sim, nburn = n_burn, sparse = F))
  
  bcf_tau_itt <- bcf_itt.tree$tau
  bcf_itt <- colMeans(bcf_tau_itt)
  tictoc::toc()
  

  # Get posterior of treatment effects
  bcf_tauhat <- bcf_itt/bcf_pic
  # driven by outliers
  # bcf_exp %>% dplyr::filter(x1 == 0 & x2 == 0) %>% dplyr::summarise(median(bcf_tauhat))
  bcf_exp <- as.data.frame(cbind(bcf_tauhat, x[-index,]))

  # repair names? !! by me 
  names(bcf_exp)[2:length(bcf_exp)] <- names(inference)[-(1:3)]
  
  ##---- Costs 
  # posterior splitting probabilities
  bcf_post_split_probs <- colMeans(bcf_itt.tree$varprb_tau)
  # check if x21 and x22 get recognized 
  plot(bcf_post_split_probs)
  # transfer probabilities into variable-specific costs
  # These costs are scalings to be applied when considering splits,
  # so the improvement on splitting on a variable is divided by its cost
  # in deciding which split to choose.
  # why between 10 and 1?
  bcf_var_costs <- scales::rescale(bcf_post_split_probs, to=c(10, 1))
  
  #####------ SBCF

  'sparse BCF'
  ##### sparse BCF -> BCF https://github.com/albicaron/SparseBCF
  tictoc::tic()
  s_bcf_itt.tree <- quiet(
    SparseBCF::SparseBCF(y[-index], z[-index], x[-index,], pihat = pihat,
                         nsim = n_sim, nburn = n_burn))

  s_bcf_tau_itt <- s_bcf_itt.tree$tau
  s_bcf_itt <- colMeans(s_bcf_tau_itt)
  tictoc::toc()

  # Get posterior of treatment effects
  s_bcf_tauhat <- s_bcf_itt/bcf_pic
  # s_bcf_exp %>% dplyr::filter(V2 == 0 & V3 == 0) %>% dplyr::summarise(median(s_bcf_tauhat))
  s_bcf_exp <- as.data.frame(cbind(s_bcf_tauhat, x[-index,]))
  # repair names? !! by me 
  names(s_bcf_exp)[2:length(s_bcf_exp)] <- names(inference)[-(1:3)]
  
  ##---- Costs 
  # posterior splitting probabilities
  s_bcf_post_split_probs <- colMeans(s_bcf_itt.tree$varprb_tau)
  # check if x21 and x22 get recognized 
  plot(s_bcf_post_split_probs)
  # transfer probabilities into variable-specific costs
  # These costs are scalings to be applied when considering splits,
  # so the improvement on splitting on a variable is divided by its cost
  # in deciding which split to choose.
  # why between 10 and 1?
  s_bcf_var_costs <- scales::rescale(s_bcf_post_split_probs, to=c(10, 1))


