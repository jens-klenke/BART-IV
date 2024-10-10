# Changing 2(b) of the Algorithm: softbart_probit for prop. of complier estimation


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
x3 <- X[, p/2+1]
x4 <- X[, p/2+1]

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

#### new testing here: calculate proportion of compliers with softbart
# prepare dataset
dat_test <- data.frame(w=factor(w[-index]), z=z[-index], x[-index,])
# Procedure according to Hill (2011) but with SoftBART-Probit instead of BART
# see https://muse.jhu.edu/pub/56/article/793357/pdf, http://cran.nexr.com/web/packages/BART/vignettes/bincat.pdf, https://arxiv.org/pdf/2210.16375
# fit BART to the observed data (W given Z and X) 
soft_test <- SoftBart::softbart_probit(w ~. , data=dat_test, test_data =  dat_test)
# Next make predictions for two datasets (Hill 2011). X is kept intact for both, however, in one
# all treatment values are set to 0, and in the other they are all set to 1. This allows BART 
# to draw from the posterior distribution  
# for E[W (1) | X]
val_covs_1 <- data.frame(w=factor(w[-index]), z=ifelse(dat_test$z==1, 1, 1), x[-index,])
soft_test_w1 <- predict(soft_test, newdata=val_covs_1)
# and E[W (0) | X]
val_covs_0 <- data.frame(w=factor(w[-index]), z=ifelse(dat_test$z==0, 0, 0), x[-index,])
soft_test_w0 <- predict(soft_test, newdata=val_covs_0)
# implying we can also obtain draws from E[W (1)−W (0) | X] for each person.
pic_sbart <- soft_test_w1$p_mean - soft_test_w0$p_mean

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
s_bcf_tauhat_2 <- s_bcf_itt/pic_sbart
# s_bcf_exp %>% dplyr::filter(V2 == 0 & V3 == 0) %>% dplyr::summarise(median(s_bcf_tauhat))
s_bcf_exp <- as.data.frame(cbind(s_bcf_tauhat, x[-index,]))
s_bcf_exp_2 <- as.data.frame(cbind(s_bcf_tauhat_2, x[-index,]))

# repair names? !! by me 
names(s_bcf_exp_2)[2:length(s_bcf_exp_2)] <- names(inference)[-(1:3)]



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

# binary tree for sparse trees
s_bcf_fit.tree_2 <- rpart(s_bcf_tauhat_2 ~ .,
                          data = s_bcf_exp_2,
                          #data = s_bcf_exp_relevant,
                          maxdepth = max_depth,
                          cp=cp,
                          minsplit=minsplit,
                          cost=var_costs2)


rpart.plot(bcf_fit.tree)
rpart.plot(s_bcf_fit.tree)
rpart.plot(s_bcf_fit.tree_2)



