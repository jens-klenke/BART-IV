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
  crayon,
  future
)