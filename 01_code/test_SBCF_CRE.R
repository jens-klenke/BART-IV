# Causal Rule Ensemble (CRE)
# https://github.com/NSAPH-Software/CRE

# playing around with discover_rules.R file of CRE-repo.
# they use randomForest (instead of single rpart tree) to discover (many more) suitable decision rules 
# afterwards, they filter out irrelevant rules 


# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)source(here::here('01_code/packages.R'))


############## Generate Data ################
# we need variables y, w, z, x

# num covariates
p <- 10
# num obs 
n <- 1000
# effects
effect_size <- 2
null <- 0



# generate covariates
get_features <- function(N, P, uncorrelated=T) {
  
  if(uncorrelated==F){
    # Generate correlated uniforms from a Gaussian Copula
    mysigma = matrix(1, P, P)
    
    for (i in 1:P) {
      for (j in 1:P) {
        mysigma[i, j] = 0.3^abs(i - j) + ifelse(i == j, 0, 0.1)
      }
    }
  }else if(uncorrelated==T){
    mysigma <- diag(1, nrow=P, ncol=P)
  }else{
    stop("check 'uncorrelated' arguemnt.")
  }
  
  
  
  mycop = MASS::mvrnorm(N, rep(0, P), Sigma = mysigma)
  unif = pnorm(mycop)
  
  
  # Transform in continuous and binary covariates
  X = matrix(NA, N, P)
  P_half <- round(P/2)
  X[, 1:P_half] = qnorm(unif[, (1:P_half)])
  X[, (P_half+1):P] = qbinom(unif[, ((P_half+1):P)], 1, 0.3)
  
  return(X)
  
}

X <- get_features(N=n, P=p, T)

# Simulate Pscore and Z
#und_lin = -0.5 + 0.2*X[, 1] + 0.1*X[, 2] + 0.4*X[, 21] + runif(n)/10
#pscore <- pnorm(und_lin)
# we stick to a perfectly random instrument
pscore <- 0.5
# Generate Random Instrument
z = rbinom(n, 1, pscore)

# Generate unit level observed exposure
compliance <- 0.75
w1 <- rbinom(n, 1, compliance)
w0 <- numeric(n)



# baseline effect
# 6 covariates are relevant for mu()
# X1, X2, X6, X7
mu = 3 + 1.5*sin(pi*X[, 1]) + 0.5*(X[, 2] - 0.5)^2 + 1.5*(2-abs(X[, 6])) + 0.5*(X[, 7] + 1) # without + 1.5*X[, 4]*(0.5*(X[, 21] + 1)

# Generate unit level potential outcome
y0 <- mu + rnorm(n)
y1 <- numeric(n)

# x6 and x7 needed for heterogeneous effects
x3 <- X[, 6]
x4 <- X[, 7]

# Generate Heterogeneity - only depends on x3 and x4  
y1[x3 == 0 & x4 == 0] <- y0[x3 == 0 & x4 == 0] + w1[x3 == 0 & x4 == 0] * effect_size
y1[x3 == 0 & x4 == 1] <- y0[x3 == 0 & x4 == 1] + w1[x3 == 0 & x4 == 1] * null
y1[x3 == 1 & x4 == 0] <- y0[x3 == 1 & x4 == 0] + w1[x3 == 1 & x4 == 0] * null
y1[x3 == 1 & x4 == 1] <- y0[x3 == 1 & x4 == 1] + w1[x3 == 1 & x4 == 1] * -effect_size

tau_true <- numeric(n)
tau_true[x3 == 0 & x4 == 0] <- w1[x3 == 0 & x4 == 0] * effect_size
tau_true[x3 == 0 & x4 == 1] <- w1[x3 == 0 & x4 == 1] * null
tau_true[x3 == 1 & x4 == 0] <- w1[x3 == 1 & x4 == 0] * null
tau_true[x3 == 1 & x4 == 1] <- w1[x3 == 1 & x4 == 1] * -effect_size


# num of obs. with effect 2: 
length(y1[x3 == 0 & x4 == 0])

# num of obs. with effect 0: 
length(y1[x3 == 0 & x4 == 1]) + length(y1[x3 == 1 & x4 == 0])

# num of obs. with effect -2: 
length(y1[x3 == 1 & x4 == 1])


# Unit level observed exposure and observed response
w <- z * w1 + (1 - z) * w0
y <- z * y1 + (1 - z) * y0
x <- X

# function params
binary = FALSE
n_burn = 3000
n_sim = 7000
inference_ratio = 0.5
max_depth = 2
cp = 0.01
minsplit = 30
adj_method = "holm"
seed = 42




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

#print('Step 0 completed')
######################################################
####  Step 1: Compute the Bayesian Causal Forest  ####
######################################################

# Compute the Propensity Score though a Logistic Regression
p.score <- glm(z ~ x[-index,],
               family = binomial,
               data = discovery)
pihat <- predict(p.score, as.data.frame(x[-index,]))

# Perform the Bayesian Causal Forest  to calculate the Proportion of Compliers (pic)
pic_bcf_tree <- SimDesign::quiet(bartCause::bartc(w[-index], z[-index], x[-index,],
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
itt_bcf <- SparseBCF::SparseBCF(y[-index], z[-index], x[-index,], 
                                pihat = pihat,
                                nsim = n_sim, nburn = n_burn,
                                sparse=F
)
bcf_tau_itt <- itt_bcf$tau
bcf_itt <- colMeans(bcf_tau_itt)

# Get posterior of treatment effects
bcf_tauhat <- bcf_itt/pic_bcf
bcf_exp <- as.data.frame(cbind(bcf_tauhat, x[-index,]))

# repair names for binary tree
names(bcf_exp)[2:length(bcf_exp)] <- names(inference)[-(1:3)]

## SBART
s_bcf_itt.tree <- SparseBCF::SparseBCF(y[-index], z[-index], x[-index,], 
                                       pihat = pihat,
                       nsim = n_sim, nburn = n_burn,
                       sparse=T,
                       OOB=T,
                       x_pred_mu = x[-index,],
                       pi_pred=pihat,
                       x_pred_tau = x[-index,],
                       save_trees_mu_dir = paste(gewtd(), "/01_code/trees/mu_trees.txt", sep="")
                       )



s_bcf_tau_itt <- s_bcf_itt.tree$tau
s_bcf_itt <- colMeans(s_bcf_tau_itt)

# Get posterior of treatment effects
s_bcf_tauhat <- s_bcf_itt/pic_bcf
# s_bcf_exp %>% dplyr::filter(V2 == 0 & V3 == 0) %>% dplyr::summarise(median(s_bcf_tauhat))
s_bcf_exp <- as.data.frame(cbind(s_bcf_tauhat, x[-index,]))

# repair names? !! by me 
names(s_bcf_exp)[2:length(s_bcf_exp)] <- names(inference)[-(1:3)]



######################################################
####  Step 2: Build a CART on the Unit Level CITT ####
######################################################
bcf_fit.tree <- rpart(bcf_tauhat ~ .,
                      data = bcf_exp,
                      maxdepth = max_depth,
                      cp=cp,
                      minsplit=minsplit)
rpart.plot(bcf_fit.tree)



# posterior splitting probabilities
post_split_probs <- colMeans(s_bcf_itt.tree$varprb_tau)
# check if x3 and x4 get recognized 
plot(post_split_probs)
# transfer probabilities into variable-specific costs
# These costs are scalings to be applied when considering splits,
# so the improvement on splitting on a variable is divided by its cost
# in deciding which split to choose.
var_costs <- scales::rescale(post_split_probs , to=c(10, 1))
# default cost value is 1 in rpart
# scaling between 10 and 1 is chosen arbitrary here, guidance?
var_costs2 <- max(post_split_probs)/post_split_probs

# binary tree for sparse trees
s_bcf_fit.tree <- rpart(s_bcf_tauhat ~ .,
                        data = s_bcf_exp,
                        #data = s_bcf_exp_relevant,
                        maxdepth = max_depth,
                        cp=cp,
                        minsplit=minsplit,
                        cost=var_costs2)
rpart.plot(s_bcf_fit.tree)




##### Discovery with randomForest as in CRE #### 

max_depth = 2
max_rules = 50
forest <- randomForest::randomForest(y=s_bcf_tauhat,
                           x=s_bcf_exp[,-1],
                           sampsize=0.5 * dim(s_bcf_exp[,-1])[1], 
                           ntree = 10,
                           maxnodes=2^max_depth,
                           nodesize=20,
                           mtry=ncol(s_bcf_exp[,-1])*2/3)

                           
treelist <- inTrees::RF2List(forest)

#### initial set of decision rules ####

extract_rules <- function(treelist, X, max_depth, digits = 2) {
  
  if (is.numeric(digits)) digits <- as.integer(abs(digits))
  levelX <- list()
  for (iX in 1:ncol(X)) levelX <- c(levelX, list(levels(X[, iX])))
  ntree <- min(treelist$ntree)
  allRulesList <- list()
  for (iTree in 1:ntree) {
    rule <- list()
    count <- 0
    rowIx <- 1
    tree <- treelist$list[[iTree]]
    if (nrow(tree) <= 1) next # skip if there is no split
    ruleSet <- vector("list", length(which(tree[, "status"] == -1)))
    for (max_length in 1:max_depth) {
      res <- inTrees::treeVisit(tree,
                                rowIx = rowIx,
                                count,
                                ruleSet,
                                rule,
                                levelX,
                                length = 0,
                                max_length = max_length,
                                digits = digits)
      allRulesList <- c(allRulesList, res$ruleSet)
    }
  }
  
  allRulesList <- allRulesList[!unlist(lapply(allRulesList, is.null))]
  rules <- inTrees::ruleList2Exec(X, allRulesList)
  return(rules)
}
rules <- extract_rules(treelist, s_bcf_exp[,-1], max_depth)
rule_counts <- table(unlist(rules))
M <- min(max_rules, length(rule_counts))
rules <- names(sort(rule_counts, decreasing = TRUE)[1:M])
M_initial <- length(rules)


#### filter out irrelevant rules #####

t_decay = 0.025
filter_irrelevant_rules <- function(rules, X, ite, t_decay) {
  
  #logger::log_debug("Filtering irrelevant rules...")
  ite_ <- ite - mean(ite)
  
  rules_matrix <- matrix(rules)
  colnames(rules_matrix) <- "condition"
  metric <- inTrees::getRuleMetric(rules_matrix,
                                   X,
                                   ite_)
  
  pruned <- inTrees::pruneRule(rules = metric,
                               X = X,
                               target = ite_,
                               maxDecay = t_decay)
  rules <- unique(pruned[, 4])
  
  #logger::log_debug("Done with filtering irrelevant rules.")
  
  return(unique(rules))
}
rules <- filter_irrelevant_rules(rules, s_bcf_exp[,-1], s_bcf_tauhat,t_decay)
M_filter1 <-length(rules)

generate_rules_matrix <- function(X, rules_list) {
  
  # Generate and Rules Matrix
  samplesize <- dim(X)[1]
  nrules <- length(rules_list)
  rules_matrix <- matrix(0, nrow = samplesize, ncol = nrules)
  for (i in 1:nrules){
    rules_matrix[eval(parse(text = rules_list[i]), list(X = X)), i] <- 1
  }
  return(rules_matrix)
}
standardize_rules_matrix <- function(rules_matrix) {
  
  samplesize <- dim(rules_matrix)[1]
  nrules <- dim(rules_matrix)[2]
  mu_rules_matrix <- apply(rules_matrix, 2, mean)
  sd_rules_matrix <- apply(rules_matrix, 2, stats::sd)
  rules_matrix_std <- matrix(0, samplesize, nrules)
  for (l in 1:ncol(rules_matrix_std)) {
    rules_matrix_std[, l] <- ((rules_matrix[, l] - mu_rules_matrix[l]) /
                                sd_rules_matrix[l])
  }
  
  return(rules_matrix_std)
}
rules_matrix <- generate_rules_matrix(s_bcf_exp[,-1], rules)

t_ext = 0.025
filter_extreme_rules <- function(rules_matrix, rules_list, t_ext) {
  
  #logger::log_debug("Filtering extreme rules...")
  
  # Identify rules with too few or too many observations
  ind <- 1:dim(rules_matrix)[2]
  sup <- apply(rules_matrix, 2, mean)
  elim <- which((sup < t_ext) | (sup > (1 - t_ext)))
  if (length(elim) > 0) ind <- ind[-elim]
  
  # Remove rules with too few/too many observations
  rules_matrix <- rules_matrix[, ind, drop = FALSE]
  rules_list <- rules_list[ind]
  colnames(rules_matrix) <- rules_list
  
  #logger::log_debug("Done with filtering extreme rules.")
  
  return(rules_matrix)
}
rules_matrix <- filter_extreme_rules(rules_matrix, rules,t_ext)
rules <- colnames(rules_matrix)
M_filter2 <- length(rules)


t_corr = 1
filter_correlated_rules <- function(rules_matrix, rules_list, t_corr) {
  
  #logger::log_debug("Filtering correlated rules...")
  
  # Identify correlated rules
  nrules <- length(rules_list)
  ind <- 1:nrules
  C <- stats::cor(rules_matrix)
  elim <- c()
  for (i in 1:(nrules - 1)) {
    elim <- c(elim,
              which(round(abs(C[i, (i + 1):nrules]), digits = 4) >= t_corr) + i)
  }
  if (length(elim) > 0) ind <- ind[-elim]
  
  # Remove correlated rules
  rules_matrix <- rules_matrix[, ind, drop = FALSE]
  rules_list <- rules_list[ind]
  colnames(rules_matrix) <- rules_list
  
  #logger::log_debug("Done with filtering correlated rules.")
  
  return(rules_matrix)
}
rules_matrix <- filter_correlated_rules(rules_matrix, rules, t_corr)
rules <- colnames(rules_matrix)
M_filter3 <- length(rules)


#### select relevant rules via post lasso ####

stability_selection = "vanilla"
cutoff = 0.6
pfer = 1
B = 20
select_rules <- function(rules_matrix, rules, ite,
                         stability_selection, cutoff, pfer, B) {
  
  #logger::log_debug("Selecting rules...")
  
  "%>%" <- magrittr::"%>%"
  
  rules_weight <- c()
  for (rule in rules) {
    rule_length <- lengths(regmatches(rule, gregexpr("&", rule))) + 1
    rules_weight <- append(rules_weight, rule_length)
  }
  R <- t(t(rules_matrix) / rules_weight)
  M <- ncol(R)
  
  if (length(rules) > 1) {
    
    if (stability_selection=="vanilla") {
      # Vanilla Stability Selection
      stability_scores <- rep(0, M)
      ite_mean <- mean(ite)
      for (i in 1:B) {
        subsample <- 0.5
        indices <- sample(1:nrow(R),
                          size = round(nrow(R) * subsample),
                          replace = FALSE)
        lasso <- glmnet::cv.glmnet(x = as.matrix(R[indices, ]),
                                   y = ite[indices] - ite_mean,
                                   alpha = 1,
                                   nfolds = 5,
                                   gamma = c(0.01, 0.05, 0.1, 0.5, 1, 5, 10),
                                   intercept = FALSE)
        non_zero_indices <- which(as.matrix(coef(lasso)) != 0) - 1
        stability_scores[non_zero_indices] <- stability_scores[non_zero_indices] + 1
      }
      stability_scores <- stability_scores / B
      rules <- colnames(R)[stability_scores >= cutoff]
      
    } else if (stability_selection=="error_control") {
      # Stability Selection with Error Control
      stab_mod <- tryCatch(
        {
          stabs::stabsel(x = as.data.frame(R),
                         y = ite - mean(ite),
                         intercept = FALSE,
                         fitfun = "glmnet.lasso",
                         cutoff = cutoff,
                         PFER = pfer)
        },
        error = function(e) {
          stop(paste(
            "Combination of `cutoff` and `pfer` not allowed. ",
            "Try to decrease the `cutoff` or increase the `pfer`. ",
            "See Stability Selection documentation for further details.",
            "\n\nOriginal Error message:", e))
        }
      )
      rules <- rules[stab_mod$selected]
      
    } else if (stability_selection=="no") {
      # LASSO
      cv_lasso <- glmnet::cv.glmnet(x = rules_matrix,
                                    y = ite - mean(ite),
                                    alpha = 1,
                                    intercept = FALSE)
      aa <- stats::coef(cv_lasso, s = cv_lasso$lambda.1se)
      index_aa <- which(aa[-1, 1] != 0)
      rule_LASSO <- data.frame(rules = rules[index_aa],
                               val = aa[index_aa + 1, 1])
      rule_LASSO <- rule_LASSO[order(-rule_LASSO[, 2]), ]
      rule_LASSO <- rule_LASSO[!is.na(rule_LASSO$rules), ]
      rules <- rule_LASSO$rules
    
      }
    
   else if (stability_selection=="softbart") {
    # SoftBart
    sbart_model <- SoftBart::softbart(X=rules_matrix,
                                      X_test=rules_matrix,
                                      Y=ite - mean(ite))
    p_probs <- SoftBart::posterior_probs(sbart_model)
    median_prob_model <- p_probs$median_probability_model
    
    ordered_rules <- data.frame(rules=rules, post_probs=p_probs$post_probs)
    rules <- list(selected_rules=rules[median_prob_model], 
                  post_probs_all=ordered_rules[order(ordered_rules$post_probs, decreasing = T),])
   }
  } else {
    rules <- NULL
  }
  
  #logger::log_debug("Done with selecting rules.")
  return(rules)
}


rules <- select_rules(rules_matrix, rules, s_bcf_tauhat,
                      stability_selection="softbart", cutoff, pfer, B)




##### Discovery with randomForest as in CRE and SoftBart but without filter-functions#### 

max_depth = 2
max_rules = 50
forest <- randomForest::randomForest(y=s_bcf_tauhat,
                                     x=s_bcf_exp[,-1],
                                     sampsize=0.5 * dim(s_bcf_exp[,-1])[1], 
                                     ntree = 10,
                                     maxnodes=2^max_depth,
                                     nodesize=20,
                                     mtry=ncol(s_bcf_exp[,-1])*2/3)


treelist <- inTrees::RF2List(forest)

#### initial set of decision rules ####

extract_rules <- function(treelist, X, max_depth, digits = 2) {
  
  if (is.numeric(digits)) digits <- as.integer(abs(digits))
  levelX <- list()
  for (iX in 1:ncol(X)) levelX <- c(levelX, list(levels(X[, iX])))
  ntree <- min(treelist$ntree)
  allRulesList <- list()
  for (iTree in 1:ntree) {
    rule <- list()
    count <- 0
    rowIx <- 1
    tree <- treelist$list[[iTree]]
    if (nrow(tree) <= 1) next # skip if there is no split
    ruleSet <- vector("list", length(which(tree[, "status"] == -1)))
    for (max_length in 1:max_depth) {
      res <- inTrees::treeVisit(tree,
                                rowIx = rowIx,
                                count,
                                ruleSet,
                                rule,
                                levelX,
                                length = 0,
                                max_length = max_length,
                                digits = digits)
      allRulesList <- c(allRulesList, res$ruleSet)
    }
  }
  
  allRulesList <- allRulesList[!unlist(lapply(allRulesList, is.null))]
  rules <- inTrees::ruleList2Exec(X, allRulesList)
  return(rules)
}
rules <- extract_rules(treelist, s_bcf_exp[,-1], max_depth)
rule_counts <- table(unlist(rules))
M <- min(max_rules, length(rule_counts))
rules <- names(sort(rule_counts, decreasing = TRUE)[1:M])
M_initial <- length(rules)

generate_rules_matrix <- function(X, rules_list) {
  
  # Generate and Rules Matrix
  samplesize <- dim(X)[1]
  nrules <- length(rules_list)
  rules_matrix <- matrix(0, nrow = samplesize, ncol = nrules)
  for (i in 1:nrules){
    rules_matrix[eval(parse(text = rules_list[i]), list(X = X)), i] <- 1
  }
  return(rules_matrix)
}
standardize_rules_matrix <- function(rules_matrix) {
  
  samplesize <- dim(rules_matrix)[1]
  nrules <- dim(rules_matrix)[2]
  mu_rules_matrix <- apply(rules_matrix, 2, mean)
  sd_rules_matrix <- apply(rules_matrix, 2, stats::sd)
  rules_matrix_std <- matrix(0, samplesize, nrules)
  for (l in 1:ncol(rules_matrix_std)) {
    rules_matrix_std[, l] <- ((rules_matrix[, l] - mu_rules_matrix[l]) /
                                sd_rules_matrix[l])
  }
  
  return(rules_matrix_std)
}
rules_matrix <- generate_rules_matrix(s_bcf_exp[,-1], rules)


#### select relevant rules via post lasso ####

#stability_selection = "vanilla"
cutoff = 0.6
pfer = 1
B = 20
select_rules <- function(rules_matrix, rules, ite,
                         stability_selection, cutoff, pfer, B) {
  
  #logger::log_debug("Selecting rules...")
  
  "%>%" <- magrittr::"%>%"
  
  rules_weight <- c()
  for (rule in rules) {
    rule_length <- lengths(regmatches(rule, gregexpr("&", rule))) + 1
    rules_weight <- append(rules_weight, rule_length)
  }
  R <- t(t(rules_matrix) / rules_weight)
  M <- ncol(R)
  
  if (length(rules) > 1) {
    
    if (stability_selection=="softbart_vanilla") {
      # Softbart Vanilla Stability Selection
      stability_scores <- rep(0, M)
      ite_mean <- mean(ite)
      for (i in 1:B) {
        subsample <- 0.5
        indices <- sample(1:nrow(R),
                          size = round(nrow(R) * subsample),
                          replace = FALSE)
        sbart_model <- SoftBart::softbart(X=as.matrix(R[indices, ]),
                           X_test=as.matrix(R[-indices, ]),
                           Y=ite[indices] - mean(ite))
        p_probs <- SoftBart::posterior_probs(sbart_model)
        median_prob_model <- p_probs$median_probability_model
        stability_scores[median_prob_model] <- stability_scores[median_prob_model] + 1
      }
      stability_scores <- stability_scores / B
      rules <- rules[stability_scores >= cutoff]
      
    } else if (stability_selection=="error_control") {
      # Stability Selection with Error Control
      stab_mod <- tryCatch(
        {
          stabs::stabsel(x = as.data.frame(R),
                         y = ite - mean(ite),
                         intercept = FALSE,
                         fitfun = "glmnet.lasso",
                         cutoff = cutoff,
                         PFER = pfer)
        },
        error = function(e) {
          stop(paste(
            "Combination of `cutoff` and `pfer` not allowed. ",
            "Try to decrease the `cutoff` or increase the `pfer`. ",
            "See Stability Selection documentation for further details.",
            "\n\nOriginal Error message:", e))
        }
      )
      rules <- rules[stab_mod$selected]
      
    } else if (stability_selection=="no") {
      # LASSO
      cv_lasso <- glmnet::cv.glmnet(x = rules_matrix,
                                    y = ite - mean(ite),
                                    alpha = 1,
                                    intercept = FALSE)
      aa <- stats::coef(cv_lasso, s = cv_lasso$lambda.1se)
      index_aa <- which(aa[-1, 1] != 0)
      rule_LASSO <- data.frame(rules = rules[index_aa],
                               val = aa[index_aa + 1, 1])
      rule_LASSO <- rule_LASSO[order(-rule_LASSO[, 2]), ]
      rule_LASSO <- rule_LASSO[!is.na(rule_LASSO$rules), ]
      rules <- rule_LASSO$rules
      
    }
    
    else if (stability_selection=="softbart") {
      # SoftBart
      sbart_model <- SoftBart::softbart(X=rules_matrix,
                                        X_test=rules_matrix,
                                        Y=ite - mean(ite))
      p_probs <- SoftBart::posterior_probs(sbart_model)
      median_prob_model <- p_probs$median_probability_model
      
      ordered_rules <- data.frame(rules=rules, post_probs=p_probs$post_probs)
      rules <- list(selected_rules=rules[median_prob_model], 
                    post_probs_all=ordered_rules[order(ordered_rules$post_probs, decreasing = T),])
    }
  } else {
    rules <- NULL
  }
  
  #logger::log_debug("Done with selecting rules.")
  return(rules)
}


rules <- select_rules(rules_matrix, rules, s_bcf_tauhat,
                      stability_selection="softbart_vanilla", cutoff, pfer, B)







