###############################################################################
## So this will be my R script for adapting the CSVAR model for Brazil data ###
## Lets try and keep this fairly tidy eh ######################################
###############################################################################
require(didehpc)
require(devtools)
require(ggplot2)
setwd("~/HOMES_drive/")
options(didehpc.username = "jd2117",didehpc.home = "~/HOMES_drive/",
        didehpc.cluster = "fi--didemrchnb", didehpc.credentials = "~/.smbcredentials")
didehpc::didehpc_config(cores = 3,parallel = FALSE)
didehpc::didehpc_config()
devtools::load_all("~/Dropbox/jeff_hiv_work/eppasm")

context::context_log_start()

root <- "contexts_3"

ctx<- context::context_save(root, packages = c("epp","eppasm","anclik","ggplot2","splines", "buildr"), 
                            package_sources = provisionr::package_sources(local = c("~/Dropbox/jeff_hiv_work/eppasm",
                                                                                    "~/Dropbox/MRes/hiv_project/epp",
                                                                                    "~/Dropbox/MRes/hiv_project/anclik/anclik")))
config <- didehpc::didehpc_config(cores = 1, parallel = FALSE)
obj <- didehpc::queue_didehpc(ctx,config)
a <- obj$enqueue(sessionInfo())
a$wait(10)

################################################################################
## Only for when making changes to the eppasm folder ###########################
################################################################################

obj$provision(refresh_drat = TRUE)

a <- obj$enqueue(sessionInfo())
a$status()
a$result()

###############################################################################

devtools::load_all("~/Dropbox/jeff_hiv_work/eppasm")
load("~/Dropbox/jeff_hiv_work/eppasm/data/brazil_spectrum_file", verbose = T)

###############################################################################
## Some of the options with the Brazil fp object:                            ##                          
## - neg_binom = True or F for poisson                                       ##
## - diagnoses_uses = whether to use diags for loglk                         ##
## - linear_diagnoses = choose the diag model, can have "knot_linear" model, ##
##                      "spline" for 7 knot spline model or "gamma" for gamma #
###############################################################################

brazil$fp$linear_diagnosis

brazil$fp$neg_binom <- FALSE

local_test_optim <- fitmod_csavr(brazil, incid_func ="ilogistic", B0 = 1e4, optfit = TRUE)

test_optim_1_knot_neg <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0 = 1e4, optfit = TRUE),
                                     name = "double_log_knot_neg_diagn")
  