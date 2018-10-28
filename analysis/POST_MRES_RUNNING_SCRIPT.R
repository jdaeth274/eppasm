###############################################################################
## So this will be my R script for adapting the CSVAR model for Brazil data ###
## Lets try and keep this fairly tidy eh ######################################
###############################################################################

require(devtools)
require(ggplot2)
devtools::load_all("~/Dropbox/jeff_hiv_work/eppasm")
load("~/Dropbox/jeff_hiv_work/r_data/brazil_spectrum_file", verbose = T)

brazil$fp$linear_diagnosis
brazil$fp$neg_binom <- FALSE

local_test_optim <- fitmod_csavr(brazil, incid_func ="ilogistic", B0 = 1e4, optfit = TRUE)
