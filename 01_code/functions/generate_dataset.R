############## Packages ################
# @title
# Generate Dataset for BCF-IV Example
#
# @description
# Function for generating a dataset for the discovery and estimation of heterogeneity 
# in the Complier Average Causal Effect in Instrumental Variable settings.
#
# @param n number of data points (default: 1000)
# @param p number of covariates (default: 10) # not correctly implemented? 
# @param rho correlation within the covariates (default: 0)
# @param null effect size for null condition (default: 0)
# @param seq effect size (default: 2) # why ? 
# @param compliance compliance rate (default: 0.75) # why ? 
# @param base_line_effect; implementation of the DGP by Caron et al. (2021)
# @param uncorrelated; implementation of the DGP by Caron et al. (2021) 
#
# @return
# A list containing the different variables in the generated dataset (y,z,w,X).
# 
# @examples
# dataset <- generate_dataset()
#
# @export
generate_dataset <- function(n = 1000, p = 100, rho = 0, null = 0, 
                             effect_size = 2, compliance = 0.75, 
                             covariates = 'cont-cov', base_line_effect = TRUE, 
                             share_d = 0.4, uncorrelated = T, confounded) {
  
  # number of discrete variables 
  n_d <- floor(p*share_d)
  
  # 
  if(!base_line_effect){
    # descrete case
    if(!covariates == 'cont-cov') {
      # Generate Variables
      mu <- rep(0, p)
      Sigma <- matrix(rho, nrow = p, ncol = p) + diag(p) * (1 - rho)
      rawvars <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
      pvars <- pnorm(rawvars)
      binomvars <- qbinom(pvars, 1, 0.5) 
      X <- binomvars
    }
    
    # continuous case
    if(covariates == 'cont-cov') {
      # Generate Variables
      mu <- rep(0, p)
      Sigma <- matrix(rho, nrow = p, ncol = p) + diag(p) * (1 - rho)
      rawvars <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
      pvars <- pnorm(rawvars)
      binomvars <- qbinom(pvars, 1, 0.5) 
      X <- cbind(binomvars[, 1:2], rawvars[, -c(1, 2)])
      }
  }
  
  if(base_line_effect){
    X <- get_features(N = n, P = p, n_d = n_d, uncorrelated = uncorrelated)
  }
  
  # x1 and x2 needed for heterogeneous effects
  x1 <- X[, 1]
  x2 <- X[, 2]
  
  # Generate unit level observed exposure
  w1 <- rbinom(n, 1, compliance)
  w0 <- numeric(n)
  
  if(!base_line_effect){
    # Generate unit level potential outcome
    y0 <- rnorm(n)
    y1 <- numeric(n)
  }
  
  if(base_line_effect){
    # Shrinkage Bayesian Causal Forests for Heterogeneous Treatment Effects Estimation
    # Alberto Caron, Gianluca Baio & Ioanna Manolopoulou 2022
    # p. 1208 sec. 5.1. Equation 17. # effect not correctly specified !
    
    if(!confounded){
      mu <-
        # n_d number of discrete variables 
        3 + 1.5*sin(pi*X[, n_d + 1]) + 0.5*(X[, n_d + 2] - 0.5)^2 + 
        1.5*(2-abs(X[, n_d + 3]))  # +
      # first and second discrete variable
      #  1.5*X[, n_d + 4]*(X[, 1] + 1) # can not find this part
      
    }
    
    if(confounded){
      mu <-
        # n_d number of discrete variables 
        3 + 1.5*sin(pi*X[, n_d + 1]) + 0.5*(X[, n_d + 2] - 0.5)^2 + 
        1.5*(2-abs(X[, n_d + 3]))  +
        # first and second discrete variable
        1.5*X[, n_d + 4]*(X[, 1] + 1) # can not find this part
    }
    
    
    # Generate unit level potential outcome
    y0 <- mu + rnorm(n)
    y1 <- numeric(n)
  }

  # Generate Heterogeneity - only depends on x1 and x2  
  y1[x1 == 0 & x2 == 0] <- y0[x1 == 0 & x2 == 0] + w1[x1 == 0 & x2 == 0] * effect_size
  y1[x1 == 0 & x2 == 1] <- y0[x1 == 0 & x2 == 1] + w1[x1 == 0 & x2 == 1] * null
  y1[x1 == 1 & x2 == 0] <- y0[x1 == 1 & x2 == 0] + w1[x1 == 1 & x2 == 0] * null
  y1[x1 == 1 & x2 == 1] <- y0[x1 == 1 & x2 == 1] + w1[x1 == 1 & x2 == 1] * -effect_size
  
  # True tau
  tau_true <- numeric(n)
  tau_true[x1 == 0 & x2 == 0] <- w1[x1 == 0 & x2 == 0] * effect_size
  tau_true[x1 == 0 & x2 == 1] <- w1[x1 == 0 & x2 == 1] * null
  tau_true[x1 == 1 & x2 == 0] <- w1[x1 == 1 & x2 == 0] * null
  tau_true[x1 == 1 & x2 == 1] <- w1[x1 == 1 & x2 == 1] * -effect_size
  
  # Generate Random Instrument
  z <- rbinom(n, 1, 0.5)
  
  # Unit level observed exposure and observed response
  w <- z * w1 + (1 - z) * w0
  y <- z * y1 + (1 - z) * y0
  
  # Observed data
  dataset <- list(y = y, z = z, w = w, w1=w1, w0=w0, X = X, tau_true = tau_true)
  
  return(dataset)
}
