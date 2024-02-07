# packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
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
  rpart.plot
)


# devtools::install_github("jaredsmurray/bcf")