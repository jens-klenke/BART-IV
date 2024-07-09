# generate covariates
#' N number of observations
#' P number of covariates
#' share_p share of discrete variables
#' uncorellates bollean if uncorrealted
get_features <- function(N, P, n_d = n_d, uncorrelated = T) {
  
  if(!uncorrelated){
    # Generate correlated uniforms from a Gaussian Copula
    mysigma = matrix(1, P, P)
    
    for (i in 1:P) {
      for (j in 1:P) {
        mysigma[i, j] = 0.3^abs(i - j) + ifelse(i == j, 0, 0.1)
      }
    }
  }else if(uncorrelated){
    mysigma <- diag(1, nrow=P, ncol=P)
  }else{
    stop("check 'uncorrelated' arguemnt.")
  }
  
  mycop = MASS::mvrnorm(N, rep(0, P), Sigma = mysigma)
  unif = pnorm(mycop)
  
  
  # Transform in continuous and binary covariates
  #X = matrix(NA, N, 10)
  D = qbinom(unif[, 1:n_d], 1, 0.3) # Very low
  C = qnorm(unif[, (n_d+1):P])
  
  
  return(cbind(D, C))
  
}
