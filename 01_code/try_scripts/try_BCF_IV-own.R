## try script for function own_bcf_iv '01_code/functions/BCF_IV-own.R' ----

# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
                 source))

# load stan models ----
stan_model_first_stage <- readRDS(here::here("05_stan_code/brms_first_stage.rds"))
stan_model_second_stage <- readRDS(here::here("05_stan_code/brms_second_stage.rds"))

# parallel plan
future::plan(multisession, workers = 10)
options(future.globals.maxSize = 2147483648) # 2GB  


# Data only try data! -> SIM Data in Sciebo
path_in <- list.files(
  here::here('00_sim_data'), recursive = TRUE, full.names = TRUE)

dataset <- readRDS(path_in[7])

y <- dataset$y
w <- dataset$w
z <- dataset$z
x <- dataset$X
tau_true <- dataset$tau_true

# Parameters,
binary = FALSE
n_burn = 1000
n_sim = 1000
inference_ratio = 0.5
max_depth = 2
cp = 0.01
minsplit = 30
adj_method = "holm"
seed = 42
cost = TRUE

x_names <- paste0('x', 1:ncol(x))
    
  ######################################################
  ####         Step 0: Initialize the Data          ####
  ######################################################
  
  # Split data into Discovery and Inference
  set.seed(seed)
  index <- sample(nrow(x), nrow(x)*inference_ratio, replace = FALSE)
  
  # Initialize total dataset
  iv.data <- as.data.frame(cbind(y, w, z, x))
  names(iv.data) <- c("y", "w", "z", x_names)
  
  # Discovery and Inference Samples
  discovery <- iv.data[-index,]
  inference <- iv.data[index,]
  
  inference_tau <- tau_true[index] #NEW
  
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
  
#  if (binary == FALSE){
    
    # Perform the Bayesian Causal Forest for the ITT
  # is that \mu and ITT ? 
    'bcf - sparse bcf' 
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
    
    
    'sparse BCF'
    ##### sparse BCF -> BCF https://github.com/albicaron/SparseBCF
    tictoc::tic()
    s_bcf_itt.tree <- quiet(
      SparseBCF::SparseBCF(y[-index], z[-index], x[-index,], pihat = pihat,
                           nsim = n_sim, nburn = n_burn)
      )
    
    s_bcf_tau_itt <- s_bcf_itt.tree$tau
    s_bcf_itt <- colMeans(s_bcf_tau_itt)
    tictoc::toc()
    
    # Get posterior of treatment effects
    s_bcf_tauhat <- s_bcf_itt/bcf_pic
    # s_bcf_exp %>% dplyr::filter(V2 == 0 & V3 == 0) %>% dplyr::summarise(median(s_bcf_tauhat))
    s_bcf_exp <- as.data.frame(cbind(s_bcf_tauhat, x[-index,]))
    
    # repair names? !! by me 
    names(s_bcf_exp)[2:length(s_bcf_exp)] <- names(inference)[-(1:3)]
    
    ######################################################
    ####  Step 2: Build a CART on the Unit Level CITT ####
    ######################################################
    'binary decision tree to discover, in an interpretable manner, the drivers of the heterogeneity??'
    if(!cost){
      bcf_fit.tree <- rpart(bcf_tauhat ~ .,
                            data = bcf_exp,
                            maxdepth = max_depth,
                            cp = cp,
                            minsplit = minsplit
                            )
      # plot tree
      rpart.plot::rpart.plot(bcf_fit.tree)
    
      # binary tree for sparse trees 
      s_bcf_fit.tree <- rpart(s_bcf_tauhat ~ .,
                              data = s_bcf_exp,
                              maxdepth = max_depth,
                              cp = cp,
                              minsplit = minsplit
                              )
      # plot tree
      rpart.plot::rpart.plot(s_bcf_fit.tree)
    }
    
    # with cost function
    if(cost){
      bcf_fit.tree <- rpart(bcf_tauhat ~ .,
                            data = bcf_exp,
                            maxdepth = max_depth,
                            cp = cp,
                            minsplit = minsplit,
                            cost = scales::rescale(
                              colMeans(bcf_itt.tree$varprb_tau), 
                              to = c(10, 1))
      )
      # plot tree
      rpart.plot::rpart.plot(bcf_fit.tree)
      
      # binary tree for sparse trees
      s_bcf_fit.tree <- rpart(s_bcf_tauhat ~ .,
                              data = s_bcf_exp,
                              maxdepth = max_depth,
                              cp = cp,
                              minsplit = minsplit,
                              cost = scales::rescale(
                                colMeans(s_bcf_itt.tree$varprb_tau), 
                                to = c(10, 1))
      )
      # plot tree
      rpart.plot::rpart.plot(s_bcf_fit.tree)
    }
    
    ######################################################
    ####  Step 3: Extract Rules and IV Estimation    ####
    ######################################################
    
    bcf_ivResults <- heterogeneous_treatment_estimation(bcf_fit.tree, inference = inference,
                                                        adj_method = adj_method,
                                                        stan_model_first_stage,
                                                        stan_model_second_stage, 
                                                        tau_true = tau_true) # NEW
    
    s_bcf_ivResults <- heterogeneous_treatment_estimation(s_bcf_fit.tree, inference = inference,
                                                          adj_method = adj_method,
                                                          stan_model_first_stage,
                                                          stan_model_second_stage,
                                                          tau_true = tau_true) # NEW
