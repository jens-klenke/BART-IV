############## Packages ################
source(here::here('01_code/packages.R'))
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
#
# @return
# A list containing the different variables in the generated dataset (y,z,w,X).
# 
# @examples
# dataset <- generate_dataset()
#
# @export
generate_dataset <- function(n = 1000, p = 100, rho = 0, null = 0, effect_size = 2, compliance = 0.75) {
  
  # Generate Variables
  mu <- rep(0, p)
  Sigma <- matrix(rho, nrow = p, ncol = p) + diag(p) * (1 - rho)
  rawvars <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  pvars <- pnorm(rawvars)
  binomvars <- qbinom(pvars, 1, 0.5) 
  X <- binomvars
  
  ########## just vectoriced? ########## 
  'why only the first 10 variables? Mistake?' 
  # x1 and x2 needed for heterogenes effects
  x1 <- X[, 1]
  x2 <- X[, 2]

  # Generate unit level observed exposure
  w1 <- rbinom(n, 1, compliance)
  w0 <- numeric(n)
  
  # Generate unit level potential outcome
  y0 <- rnorm(n)
  y1 <- numeric(n)
  
  # Generate Heterogeneity - only depends on x1 and x2  
  y1[x1 == 0 & x2 == 0] <- y0[x1 == 0 & x2 == 0] + w1[x1 == 0 & x2 == 0] * effect_size
  y1[x1 == 0 & x2 == 1] <- y0[x1 == 0 & x2 == 1] + w1[x1 == 0 & x2 == 1] * null
  y1[x1 == 1 & x2 == 0] <- y0[x1 == 1 & x2 == 0] + w1[x1 == 1 & x2 == 0] * null
  y1[x1 == 1 & x2 == 1] <- y0[x1 == 1 & x2 == 1] + w1[x1 == 1 & x2 == 1] * -effect_size
  
  # Generate Random Instrument
  z <- rbinom(n, 1, 0.5)
  
  # Unit level observed exposure and observed response
  w <- z * w1 + (1 - z) * w0
  y <- z * y1 + (1 - z) * y0
  
  # Observed data
  dataset <- list(y = y, z = z, w = w, X = X)
  
  return(dataset)
}

# generating
dataset <- generate_dataset()


