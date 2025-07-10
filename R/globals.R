# This file is used to declare global variables to satisfy R CMD check.
# These variables are used in dplyr and ggplot2 pipelines and are not
# actual global variables in the user's session.

utils::globalVariables(c(
   ".", ".data", "N_per_Arm", "N_per_Group", "N_per_Stratum", "Power",
   "mu0_hat", "Y_rmst_mean"
))

