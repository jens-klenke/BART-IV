# packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tibble,
  ggplot2,
  furrr,
  Hmisc,
  magrittr,
  MASS,
  stats,
  invgamma,
  splines,
  MCMCpack,
  BayesTree,
  dbarts,
  bcf,
  crayon,
  future,
  rpart,
  AER,
  bartCause,
  pryr,
  stringr,
  kableExtra,
  rpart.plot,
  BayesIV, # https://github.com/fbargaglistoffi/BCF-IV/tree/master
  SparseBCF # https://github.com/albicaron/SparseBCF
)


# devtools::install_github("jaredsmurray/bcf")