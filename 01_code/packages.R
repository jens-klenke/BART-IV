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
  dpylr,
  tidyr,
  stringr,
  kableExtra,
  rpart.plot,
  devtools,
  BayesIV, # https://github.com/fbargaglistoffi/BCF-IV/tree/master install_github("fbargaglistoffi/BCF-IV", ref="master")
  SparseBCF # https://github.com/albicaron/SparseBCF   devtools::install_github("albicaron/SparseBCF")
)


# devtools::install_github("jaredsmurray/bcf")