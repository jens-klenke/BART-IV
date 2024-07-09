#'
#'
#'

own_bcf_iv <- function(y, w, z, x, binary = FALSE, n_burn = 3000, n_sim = 7000, 
                   inference_ratio = 0.5, max_depth = 2, cp = 0.01, 
                   minsplit = 30, adj_method = "holm", seed = 42, cost = TRUE) {
  
  ######################################################
  ####         Step 0: Initialize the Data          ####
  ######################################################
  
  # Split data into Discovery and Inference
  set.seed(seed)
  index <- sample(nrow(x), nrow(x)*inference_ratio, replace=FALSE)
  
  # Initialize total dataset
  iv.data <- as.data.frame(cbind(y, w, z, x))
  names(iv.data) <- c("y", "w", "z", paste0('x', 1:ncol(x))) # names for covariates
  
  # Discovery and Inference Samples
  discovery <- iv.data[-index,]
  # 'names not fitting' # binary tree starts with V2 and not V4 in covariates
  inference <- iv.data[index,] 
  
  # print('Step 0 completed')
  ######################################################
  ####  Step 1: Compute the Bayesian Causal Forest  ####
  ######################################################
  
  # Compute the Propensity Score though a Logistic Regression
  p.score <- glm(z ~ x[-index,],
                 family = binomial,
                 data = discovery)
  pihat <- predict(p.score, as.data.frame(x[-index,]))
  
  # Perform the Bayesian Causal Forest  to calculate the Proportion of Compliers (pic)
  pic_bcf_tree <- quiet(bartCause::bartc(w[-index], z[-index], x[-index,],
                                         n.samples = n_sim, n.burn = n_burn, 
                                         n.chains = 2L))
  tau_bcf_pic <- bartCause::extract(pic_bcf_tree, type = "ite")
  pic_bcf <- apply(tau_bcf_pic, 2, mean)
  # workaround! 
  pic_bcf[pic_bcf == 0] <- 1e-2 
  
  ######################################################
  ####     Continuous and Discrete Outcomes         ####
  ######################################################
  
  # Perform the Bayesian Causal Forest for the ITT
  bcf_itt.tree <- quiet(
    SparseBCF::SparseBCF(y[-index], z[-index], x[-index,], pihat = pihat,
                         nsim = n_sim, nburn = n_burn, sparse = F))
  
  bcf_tau_itt <- bcf_itt.tree$tau
  bcf_itt <- colMeans(bcf_tau_itt)
  
  # posterior splitting probabilities
  bcf_post_split_probs <- colMeans(bcf_itt.tree$varprb_tau)
  
  # Get posterior of treatment effects
  bcf_tauhat <- bcf_itt/pic_bcf
  bcf_exp <- as.data.frame(cbind(bcf_tauhat, x[-index,]))
  
  # repair names for binary tree
  names(bcf_exp)[2:length(bcf_exp)] <- names(inference)[-(1:3)]
  
  ## SBART
  s_bcf_itt.tree <- quiet(
    SparseBCF::SparseBCF(y[-index], z[-index], x[-index,], pihat = pihat,
                         nsim = n_sim, nburn = n_burn)
  )
  
  s_bcf_tau_itt <- s_bcf_itt.tree$tau
  s_bcf_itt <- colMeans(s_bcf_tau_itt)
  
  # posterior splitting probabilities
  s_bcf_post_split_probs <- colMeans(s_bcf_itt.tree$varprb_tau)
  
  # Get posterior of treatment effects
  s_bcf_tauhat <- s_bcf_itt/pic_bcf
  # s_bcf_exp %>% dplyr::filter(V2 == 0 & V3 == 0) %>% dplyr::summarise(median(s_bcf_tauhat))
  s_bcf_exp <- as.data.frame(cbind(s_bcf_tauhat, x[-index,]))
  
  # repair names? !! by me 
  names(s_bcf_exp)[2:length(s_bcf_exp)] <- names(inference)[-(1:3)]
  
  # print('Step 1 completed')
  ######################################################
  ####  Step 2: Build a CART on the Unit Level CITT ####
  ######################################################
  
  # without cost function
  if(!cost){
    bcf_fit.tree <- rpart(bcf_tauhat ~ .,
                          data = bcf_exp,
                          maxdepth = max_depth,
                          cp = cp,
                          minsplit = minsplit)
    
    # binary tree for sparse trees
    s_bcf_fit.tree <- rpart(s_bcf_tauhat ~ .,
                            data = s_bcf_exp,
                            maxdepth = max_depth,
                            cp = cp,
                            minsplit = minsplit)
  }
  
  # with cost function
  if(cost){
    bcf_fit.tree <- rpart(bcf_tauhat ~ .,
                          data = bcf_exp,
                          maxdepth = max_depth,
                          cp = cp,
                          minsplit = minsplit,
                          cost = (max(bcf_post_split_probs)/bcf_post_split_probs)
                          )
    
    # binary tree for sparse trees
    s_bcf_fit.tree <- rpart(s_bcf_tauhat ~ .,
                            data = s_bcf_exp,
                            maxdepth = max_depth,
                            cp = cp,
                            minsplit = minsplit,
                            cost = (max(s_bcf_post_split_probs)/s_bcf_post_split_probs)
                            )
  }
  
  # print('Step 2 completed')
  ######################################################
  ####  Step 3: Extract the Causal Rules (Nodes)    ####
  ######################################################

  bcf_ivResults <- extract_causal_rules(bcf_fit.tree, inference = inference, adj_method = adj_method)
  
  s_bcf_ivResults <- extract_causal_rules(s_bcf_fit.tree, inference = inference, adj_method = adj_method)

  # print('Step 3 completed')
  ######################################################
  ####             Step 4: Return results           ####
  ######################################################
  #  print('before returning')
  
  return(
    list(
    'bcf_results' = bcf_ivResults,
    's_bcf_results' = s_bcf_ivResults,
    'pic' = pic_bcf,
    'bcf_itt' = bcf_itt,
    's_bcf_itt'= s_bcf_itt,
    'bcf_tauhat' = bcf_tauhat,
    's_bcf_tauhat' = s_bcf_tauhat,
    'bcf_exp' = bcf_exp,
    's_bcf_exp' = s_bcf_exp
    )
  )
}
