#####################################################################################
## Strem line with rw ###############################################################
#####################################################################################
setwd("X:")
options(didehpc.username = "jd2117",didehpc.home = "X:/simpleepp",didehpc.cluster = "fi--didemrchnb")

didehpc::didehpc_config(cores = 3,parallel = FALSE)
?didehpc::didehpc_config

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## !!!!!!!!!!!!!!!!!!!!!!!! Remember to turn on pulse secure at this point !!!!!!!!!!
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
didehpc::didehpc_config()

context::context_log_start()

root <- "contexts_2"

ctx<- context::context_save(root, packages = c("epp","eppasm","anclik"), 
                            package_sources = provisionr::package_sources(local = c("C:/Users/josh/Dropbox/hiv_project/eppasm",
                                                                                    "C:/Users/josh/Dropbox/hiv_project/epp",
                                                                                    "C:/Users/josh/Dropbox/hiv_project/anclik/anclik")))
config <- didehpc::didehpc_config(cores = 1, parallel = FALSE)
obj <- didehpc::queue_didehpc(ctx,config)


######################################## Refresh only !!!!! #######################################
obj$provision(refresh_drat = TRUE)

a <- obj$enqueue(sessionInfo())
a$status()
a$result()

#######################################################################################################
## Lets load up and prepare the fp file ###############################################################
#######################################################################################################

library(epp)

library(magrittr)
library(broom)
library(ggplot2)
devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")
devtools::build("C:/Users/josh/Dropbox/hiv_project/eppasm")

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/my_created_csvard_for_brazil_WITH_ART",verbose = T)

############################################################################################################
## NOw lets run the thing !!! ##############################################################################
############################################################################################################

brazil_pjnz <- "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/Brazil_2017_final.PJNZ"

brazil_fp <- prepare_directincid(brazil_pjnz)
brazil_fp$artmx_timerr <- rep(1.0, brazil_fp$ss$PROJ_YEARS)
brazil_fp$t_diagn_start <- 11L   # assume diagnoses starts in 1980


brazil_fp$relinfectART <- 0.3
brazil_fp$tsEpidemicStart <- 1970.5
brazil_fp$likelihood_cd4 <- TRUE
brazil_fp$artinit_use <- FALSE
brazil_fp$aidsdeaths <- TRUE
brazil_fp$time_at_which_get_cd4_counts <- 2001L
brazil_fp$diagn_rate <- array(0.8, c(dim(brazil_fp$cd4_mort), brazil_fp$ss$PROJ_YEARS))
brazil_fp$stages <- 4L

brazil <- list(fp = brazil_fp, csavrd = brazil_csvard_test)

#########################################################################################################
## Now lets get the right Spectrum ART data in there ####################################################
#########################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/spectrum_ART_numbers",verbose = T)

brazil$fp$art15plus_num <- spectrum_art_numbers

#########################################################################################################
## NOw lets get the corrected ART mortality numbers on there ############################################
#########################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/2015_vals_corrected_ART_mort_brazil", verbose = T)

brazil$fp$art_mort <- test_art_mort

#########################################################################################################
## Now lets get the diagnosis model in there ############################################################
#########################################################################################################

brazil$fp$linear_diagnosis <- "knot_linear"

##########################################################################################################
## Now we can get the art cd4 values in if we want to fit with these, (R MODEL ONLY)!! ###################
##########################################################################################################
load("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/cd4_at_ART_MATRIX", verbose = T)

brazil$fp$tara_art_cd4 <- art_cd4_counts
brazil$fp$use_cd4_art <- TRUE
brazil$fp$time_cd4_stop <- 2015L
brazil$fp$artcd4elig_idx <- rep(1L, 52)

##########################################################################################################
## Lets get our function for plotting the resulting trends ###############################################
##########################################################################################################

plot_undiagnosed <- function(optim_output,diag_start = 1980, art_start = 1996, model_labs, xlimits = c(1970,2015),
                             knot_linear = TRUE){
  tot_undiag_data <- NULL
  input_label <- model_labs
  tot_deaths <- NULL
  tot_diagnoses <- NULL
  tot_diag_rate <- NULL
  art_init_tot <- NULL
  tot_divided <- NULL
  tot_incid <- NULL
  
  for(i in 1:length(optim_output)){
    if(length(optim_output) == 1){
      list_version <- attributes(optim_output[[1]]$mod)
    }else{
      list_version <- attributes(optim_output[[i]]$mod)
    }
    undiag_percent <- (apply(list_version$hivpop,4,sum) - apply(list_version$diagnpop,4,sum)) / 
      (apply(list_version$hivpop,4,sum) + apply(list_version$artpop,5,sum)) * 100
    
    percent_dat <- cbind.data.frame(undiag_percent,rep(input_label[i],length(undiag_percent)),c(1970:2021))
    ## deaths
    
    deaths <- colSums(list_version$hivdeaths,,2)
    deaths_df <- cbind.data.frame(deaths, rep(input_label[i], length(deaths)), c(1970:2021))
    
    tot_deaths <- rbind.data.frame(tot_deaths, deaths_df)
    ## diagnoses
    
    diagnoses <- colSums(list_version$diagnoses,, 3)
    diagnoses_df <- cbind.data.frame(diagnoses, rep(input_label[i], length(diagnoses)), c(1970:2021))
    
    tot_diagnoses <- rbind.data.frame(tot_diagnoses, diagnoses_df)
    
    tot_undiag_data <- rbind.data.frame(tot_undiag_data,percent_dat)
    
    ## Diag rate 
    if(knot_linear == TRUE){
    diag_rate <- optim_output[[i]]$par[(length(optim_output[[i]]$par) - 4) : length(optim_output[[i]]$par)]
    knots <- c(1986, 1996, 2000, 2009, 2015)
    
    diagn_trend <- approx(knots, diag_rate, 1970:2021, rule = 2)$y
    }else{
      diag_rate <- optim_output[[i]]$par[(length(optim_output[[i]]$par) - 6) : length(optim_output[[i]]$par)]
      diag_rate <- c(0,0,diag_rate)
      nk <- 9 # number of splines
      dk <- diff(range(seq(1970,2021,0.1)))/(nk-3)
      knots <- c(1931.75,1944.5,1957.25,1970, 1980, 1985, 2001, 2009, 2015, 2021, 2033.75, 2046.5, 2059.25)#1970 + -3:nk*dk
      step_vector <- seq(1970,2021,0.1)
      
      Xsp <- splines::splineDesign(knots, step_vector , ord=4)
      diagn_trend <- Xsp %*% diag_rate
      diagn_trend <- diagn_trend[0:51*10 +1]
    }
    diagn_df <- cbind.data.frame(diagn_trend, rep(input_label[i], 52), c(1970:2021))
  
    tot_diag_rate <- rbind.data.frame(tot_diag_rate, diagn_df)
    
    ## ARt init cats
    get_dist_per_year <- array(0, dim = c(7,52))
    get_dist_4_cat <- array(0, dim = c(4, 52))
    
    for(j in 1:52){
      for(ii in 1:7){
        
        get_dist_per_year[ii, j] <- sum(list_version$artinits[ii,,,j])
        
      }
      get_dist_4_cat[1,j] <- get_dist_per_year[1,j]
      get_dist_4_cat[2, j] <- get_dist_per_year[2, j]
      get_dist_4_cat[3, j] <- get_dist_per_year[3, j] + get_dist_per_year[4, j]
      get_dist_4_cat[4, j] <- get_dist_per_year[5, j] + get_dist_per_year[6, j] + get_dist_per_year[7, j]
      
      
      
    }
    
    get_dist_4_cat <- t(get_dist_4_cat)
    get_dist_4_cat <- data.frame(get_dist_4_cat)
    get_dist_4_cat$time <- c(1970:2021)
    
    funky_dist_4 <- reshape2::melt(get_dist_4_cat, id = "time")
    
    funky_dist_4$input <- rep(input_label[i], nrow(funky_dist_4))
    
    
    art_init_tot <- rbind.data.frame(art_init_tot, funky_dist_4)
    
    ### Deaths by wether on ART on not
    
    deaths_by_hiv <- colSums(list_version$hivpopdeaths,,3)
    deaths_on_art <- colSums(list_version$artpopdeaths,,4)
    
    deaths_divided_df <- cbind.data.frame(deaths_by_hiv, deaths_on_art, c(1970:2021))
    names(deaths_divided_df) <- c("hiv","art","time")
    
    deaths_funky <- reshape2::melt(deaths_divided_df, id = "time")
    deaths_funky$input <- rep(input_label[i], nrow(deaths_funky))
    
    tot_divided <- rbind.data.frame(tot_divided, deaths_funky)
    
    ## INcidence 
    
    incid_df <- cbind.data.frame(list_version$incid15to49, rep(input_label[i], 52), c(1970:2021))
    tot_incid <- rbind.data.frame(tot_incid, incid_df)
    
    
    
  }
  
  names(tot_undiag_data) <- c("undiagnosed","input","year")
  names(tot_deaths) <- c("deaths", "input", "year")
  names(tot_diagnoses) <-c("diagnoses", "input", "year")
  names(tot_diag_rate) <-c("rate","input","year")
  names(art_init_tot) <- c("time","cat","val","input")
  names(tot_divided) <- c("time","class","deaths","input")
  names(tot_incid) <- c("incidence","input","time")
  
  art_init_tot$cat_input <- paste(art_init_tot$cat, art_init_tot$input, sep = "")
  tot_divided$class_input <- paste(tot_divided$class, tot_divided$input, sep ="")
  
  
  
  
  tot_deaths <- merge(tot_deaths, optim_output[[1]]$likdat$aidsdeaths, all = TRUE)
  tot_diagnoses <- merge(tot_diagnoses, optim_output[[1]]$likdat$diagnoses, all = TRUE)
  
  undiag_plot <- ggplot(data = tot_undiag_data,aes(x=year,y=undiagnosed,group=input))+geom_line(aes(colour=input),size=1.05)+
    labs(x="Time",y="Percent of HIV +ve population undiagnosed") + 
    geom_vline(xintercept = diag_start,col="midnightblue",size = 0.75) +
    geom_vline(xintercept = art_start,col="midnightblue",size=0.75) + coord_cartesian(ylim = c(0,100), xlim = xlimits)
  
  deaths_plot <- ggplot(data = tot_deaths, aes(x = year, y =deaths, group = input)) + geom_line(aes(colour = input), size = 1.05)+
    geom_point(aes(x = year, y = aidsdeaths), colour = "midnightblue", size = 1.5) +
    labs(x = "Time", y= "Number of AIDS deaths") + coord_cartesian(xlim = xlimits, ylim = c(0, 30000))
  
  diagnoses_plots <- ggplot(data = tot_diagnoses, aes(x = year, y = diagnoses, group = input)) +
    geom_line(aes(colour = input), size = 1.05)+
    geom_point(aes(x = year, y = total_cases), size = 1.5, colour = "midnightblue")+
    labs(x = "Time", y = "Diagnoses")  + coord_cartesian(xlim = xlimits, ylim = c(0, 80000))
  
  diag_rate_plot <- ggplot(data = tot_diag_rate, aes(x = year, y = rate, group = input)) + geom_line(aes(colour = input), size = 1.05)+
    labs(x = "time", y = "rate")  + coord_cartesian(xlim = xlimits, ylim = c(-20, 20))
  
  art_init_rate_plot <- ggplot(data = art_init_tot, aes(x = time, y = val, group = cat_input)) +
    geom_line(aes(colour = cat, linetype = input), size = 1.01) +
    labs(x = "time", y = "init numbers per class")  + coord_cartesian(xlim = xlimits, ylim = c(0, 20000))
  
  tot_divided_plot <- ggplot(data = tot_divided, aes(x = time, y = deaths, group = class_input)) +
    geom_line(aes(colour = class, linetype = input), size = 1.01) +
    labs(x = "time", y = "deaths from AIDS")  + coord_cartesian(xlim = xlimits, ylim = c(0, 20000))
  
  tot_incid_plot <- ggplot(data = tot_incid, aes(x = time, y = incidence, group = input)) +
    geom_line(aes(colour = input), size = 1.02) +
    labs(x = "time", y = "incidence") + coord_cartesian(xlim = xlimits, ylim = c(0, 1e-03))
  
  combined_plot <- ggpubr::ggarrange(undiag_plot, deaths_plot, diagnoses_plots,
                                     tot_incid_plot, art_init_rate_plot, tot_divided_plot,
                                     ncol = 2, nrow = 3, align = c("v"), labels = c("Undiagnosed (%)",
                                                                                    "AIDS deaths", "New diagnoses",
                                                                                    "Incidence", "ART initation by CD4 class",
                                                                                    "Deaths by treatment")) 
  
  return(list(undiag_df = tot_undiag_data,undiag_plot = undiag_plot, diagnoses_df = tot_diagnoses, diagnoses_plot = diagnoses_plots,
              deaths_df = tot_deaths, deaths_plot = deaths_plot, rate_plot = diag_rate_plot, combined_plot = combined_plot,
              art_inits = art_init_rate_plot, deaths_div = tot_divided_plot, incid_plot = tot_incid_plot))
  
}

##########################################################################################################
## Now let's test this out for the RW model ##############################################################
##########################################################################################################
brazil$fp$diagnoses_uses <- TRUE
brazil$fp$tARTstart <- 28L

save(brazil, file = "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/updated_brazil_fp_csavrd")

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/updated_brazil_fp_csavrd", verbose = T)

devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")
devtools::build("C:/Users/josh/Dropbox/hiv_project/eppasm")

brazil_opt_RW <- obj$enqueue(fitmod_csavr(brazil, eppmod="logrw", B0=1e4, optfit=TRUE),
                                         name = "logRW_model")
brazil_opt_RW$status()                                         
brazil_opt_RW$log()
brazilo_opt_log_r <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE),
                                 name = "optim_knot_linear_mod_3")
brazilo_opt_log_r$status()
brazilo_opt_log_r$log()
brazilo_opt_dooblay_incid <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0 = 1e4, optfit = TRUE),
                                         name = "incid double logistic")
brazilo_opt_dooblay_incid$status()
brazilo_opt_dooblay_incid$log()

brazil_r_log_RW <- obj$enqueue(fitmod_csavr(brazil, eppmod = "rlogrw", B0 = 1e4, optfit = TRUE),
                               name = "r log RW")

brazil_r_log_RW$status()
brazil_r_log_RW$log()
opt_rw_res <- brazil_opt_RW$result()
opt_r_log_res <- brazilo_opt_log_r$result()
opt_dooblay_incid <- brazilo_opt_dooblay_incid$result()

opt_res_list <- list(opt_rw_res, opt_r_log_res, opt_dooblay_incid)

output_test <- plot_undiagnosed(opt_res_list, model_labs = c("RW", "rlog","doubleincid"))

output_test$combined_plot

output_test$rate_plot


############################################################################################
## Running optim locally ###################################################################
############################################################################################

brazil$fp$incid_func <- "Null"
lof_r_local <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrw", B0 = 1e4, optfit = TRUE),
                           name = "RW_on_R_optim")
lof_r_local$status()

lof_rw_incid_loc <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logrw", B0 = 1e4, optfit = TRUE),
                                name = "RW_on_incid_optim")
lof_rw_incid_loc$status()
lof_rw_incid_loc$log()

lof_spline_r <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrspline", B0 = 1e4, optfit = TRUE),
                            name = "Spline_on_R_optim")
lof_spline_r$status()

lof_spline_incid <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logspline", B0 = 1e4, optfit = TRUE),
                                name = "Spline_on_incid_optim")
lof_spline_incid$status()

rw_on_inc_ouput <- list(and_here_i_lay, lof_r_local$result(),
                        lof_spline_incid$result(),lof_spline_r$result())

this_is_my_rock <- plot_undiagnosed(rw_on_inc_ouput, model_labs = c("RW_on_incid","RW_on_R","Spline_on_incid", "spline_on_r"),
                                    xlimits = c(1970,2021))
this_is_my_rock$combined_plot

devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")

and_here_i_lay <- fitmod_csavr(brazil, incid_func = "incid_logrw", B0 = 1e4, optfit = TRUE)

punctured_bicycle <- fitmod_csavr(brazil, eppmod = "logrw", B0 = 1e4, optfit = TRUE)

############################################################################################
## Running the updated RW's with IMIS ######################################################
############################################################################################

brazil_rw_incid <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logrw", B0=1e4,
                                            B=1e3, B.re=3e3, opt_iter=1:3*5),
                               name = "log_RW__INCID_IMIS_fit")
brazil_rw_incid$status()
brazil_rw_incid$log()
brazil_rw_incid_id <- brazil_rw_incid$id
save(brazil_rw_incid_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/IMIS_28_7_2018/RW_on_INCID_ID")

## RW on r(t)
brazil_rw_kappa <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrw", B0=1e4,
                                            B=1e3, B.re=3e3, opt_iter=1:3*5),
                               name = "log_RW_KAPPA_IMIS_fit")
brazil_rw_kappa$status()
brazil_rw_kappa$log()
brazil_rw_kappa_id <- brazil_rw_kappa$id
save(brazil_rw_kappa_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/IMIS_28_7_2018/RW_on_KAPPA_ID")

## Spline on incid 

brazil_sp_incid <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logspline", B0=1e4,
                                            B=1e3, B.re=3e3, opt_iter=1:3*5),
                               name = "log_SPLINE_INCID_IMIS_fit")
brazil_sp_incid$status()
brazil_sp_incid$log()
brazil_sp_incid_id <- brazil_sp_incid$id
save(brazil_sp_incid_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/IMIS_28_7_2018/SP_on_INCID_ID")

## Spline on kappa
brazil_sp_kappa <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrspline", B0=1e4,
                                            B=1e3, B.re=3e3, opt_iter=1:3*5),
                               name = "log_SPLINE_KAPPA_IMIS_fit")
brazil_sp_kappa$status()
brazil_sp_kappa$log()
brazil_sp_kappa_id <- brazil_sp_kappa$id
save(brazil_sp_kappa_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/IMIS_28_7_2018/SP_on_KAPPA_ID")

#### load up the results #####

id_path <- "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/IMIS_28_7_2018/"
ids <- list.files(id_path, full.names = TRUE)
for(i in 1:length(ids))
  load(ids[[i]], verbose = T)

rw_incid <- obj$task_get(brazil_rw_incid_id)
rw_kappa <- obj$task_get(brazil_rw_kappa_id)
sp_incid <- obj$task_get(brazil_sp_incid_id)
sp_kappa <- obj$task_get(brazil_sp_kappa_id)

rw_incid_res <- rw_incid$result()
rw_kappa_res <- rw_kappa$result()
sp_incid_res <- sp_incid$result()
sp_kappa_res <- sp_kappa$result()

brazil_out1_cd4 <- tidy(rw_incid_res) %>% data.frame(model = "RW_incid", .)
brazil_out2_cd4 <- tidy(rw_kappa_res) %>% data.frame(model = "RW_kappa", .)
brazil_out3_cd4 <- tidy(sp_incid_res) %>% data.frame(model = "Sp_incid", .)
brazil_out4_cd4 <- tidy(sp_kappa_res) %>% data.frame(model = "Sp_kappa", .)
brazil_out_cd4 <- rbind(brazil_out2_cd4, brazil_out3_cd4, brazil_out4_cd4)

save(brazil_out_cd4,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/IMIS_28_7_2018/tidied_results_for_sp_RW_NO_RW_incid")


knot_linear_diag_plot <- ggplot(subset(brazil_out_cd4, year %in% 1975:2017), 
                                aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line(size=1.05) + geom_ribbon(alpha=0.2) + 
  facet_wrap(~outcome, scales="free") +
  geom_point(aes(y=lik_data), col="darkred", size=0.5) +
  labs(title = "Spectrum ART numbers, RW and Spline Compo")

knot_linear_diag_plot

############################################################################################
## Now lets run some optim fits using different knotted splines and penalty orders #########
############################################################################################

## first lets run the logrspline at 7 and 15 knots, with penords 1 and 2 


sp_kappa_15_1 <- obj$enqueue(fitmod_csavr(brazil, eppmod="logrspline", numKnots = 15,
                                          rtpenord = 1, B0=2e4, optfit=TRUE),
                             name = "sp_kappa_15_knot_pen_1")
sp_kappa_15_1$status()
sp_kappa_15_1$log()

sp_kappa_15_2 <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrspline", numKnots = 15,
                                          B0 = 1e4, optfit = TRUE),
                             name = "sp_kappa_15_knot_pen_2")

sp_kappa_15_2$status()
sp_kappa_15_2$log()

sp_kappa_7_1 <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrspline", numKnots = 7,
                                         rtpenord = 1, B0 = 1e4, optfit = TRUE),
                            name = "sp_kappa_7_knot_pen_1")
sp_kappa_7_1$status()
sp_kappa_7_1$log()

sp_kappa_7_2 <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrspline", numKnots = 7,
                                         rtpenord = 2, B0 = 1e4, optfit = TRUE),
                            name = "sp_kappa_7_knot_pen_2")
sp_kappa_7_2$status()
sp_kappa_7_2$log()

splines_on_kappa_out <- list(sp_kappa_15_1$result(), sp_kappa_15_2$result(),
                             sp_kappa_7_1$result(), sp_kappa_7_2$result())
splines_on_kappa <- plot_undiagnosed(splines_on_kappa_out, model_labs = c("15_knot_pen_1","15_knot_pen_2","7_knot_pen_1",
                                                                          "7_knot_pen_2"))
splines_on_kappa$combined_plot

save(splines_on_kappa_out,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_29_7_2018/sp_15_and_7_on_kappa_results")



############################################################################################
## Now lets set up the same as above but for the incid runs ################################
############################################################################################

sp_incid_15_1 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logspline", numKnots = 15,
                                          rtpenord = 1, B0 = 2e4, optfit = TRUE),
                             name = "sp_incid_15_knot_pen_1")
sp_incid_15_1$status()
sp_incid_15_1$log()

sp_incid_15_2 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logspline", numKnots = 15,
                                          rtpenord = 2, B0 = 5e5, optfit = TRUE),
                             name = "sp_incid_15_knot_pen_2")
sp_incid_15_2$status()
sp_incid_15_2$log()
sp_incid_15_id <- sp_incid_15_2$id
save(sp_incid_15_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/opitm_sp_15_knot_2nd_order_id")
load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/opitm_sp_15_knot_2nd_order_id")
sp_incid_15_2 <- obj$task_get(sp_incid_15_id)
sp_incid_15_2$status()
sp_incid_15_2$log()

sp_incid_7_1 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logspline", numKnots = 7,
                                         rtpenord = 1, B0 = 1e4, optfit = TRUE),
                            name = "sp_incid_7_knot_pen_1")
sp_incid_7_1$status()
sp_incid_7_1$log()

sp_incid_7_2 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logspline", numKnots = 7,
                                         rtpenord = 2, B0 = 5e4, optfit = TRUE),
                            name = "sp_incid_7_knot_pen_2")
sp_incid_7_2$status()
sp_incid_7_2$log()

splines_on_incid_out <- list(sp_incid_15_1$result(), sp_incid_15_2$result(),
                             sp_incid_7_1$result(), sp_incid_7_2$result())
splines_on_incid <- plot_undiagnosed(splines_on_incid_out, model_labs = c("15_knot_pen_1","15_knot_pen_2",
                                                                          "7_knot_pen_1","7_knot_pen_2"))
splines_on_incid$combined_plot

save(splines_on_incid_out,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_29_7_2018/sp_15_and_7_on_INCID_results")

############################################################################################
## RUN IMIS fits on cluster, updated number of knots #######################################
############################################################################################

brazil_fit1_updated_art_knot_linear <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrw", B0=1e4,
                                                                B=1e3, B.re=3e3, opt_iter=1:3*5),
                                                   name = "log_RW_IMIS_fit")

brazil_fit1_updated_art_knot_linear$log()
brazil_fit1_updated_art_knot_linear$status()
brazil_fit1_updated_art_knot_linear_id <- brazil_fit1_updated_art_knot_linear$id
save(brazil_fit1_updated_art_knot_linear_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/SPECTRUM_ART_logrw")


brazil_fit2_updated_art_knot_linear <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic",
                                                                B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5),
                                                   name = "double_log_CD4_diag_deaths_knot_linear")
brazil_fit2_updated_art_knot_linear$status()
brazil_fit2_updated_art_knot_linear_id <- brazil_fit2_updated_art_knot_linear$id
save(brazil_fit2_updated_art_knot_linear_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/SPECTRUM_ART_double_log_incid")


## fit logistic model for transimssion rate (r(t))
brazil_fit3_updated_art_knot_linear <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5),
                                                   name = "rlog_CD4_diag_deaths_knot_linear")
brazil_fit3_updated_art_knot_linear$status()
brazil_fit3_updated_art_knot_linear_id <- brazil_fit3_updated_art_knot_linear$id
save(brazil_fit3_updated_art_knot_linear_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/SPECTRUM_ART_r_log")

path_to_victory <- "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/"
minor_slip_up <- list.files(path_to_victory, full.names = TRUE)
for(i in 1:length(minor_slip_up))
  load(minor_slip_up[[i]], verbose = TRUE)

brazil_fit1_updated_art_knot_linear <- obj$task_get(brazil_fit1_updated_art_knot_linear_id)
brazil_fit2_updated_art_knot_linear <- obj$task_get(brazil_fit2_updated_art_knot_linear_id)
brazil_fit3_updated_art_knot_linear <- obj$task_get(brazil_fit3_updated_art_knot_linear_id)

imis_log_rw <- brazil_fit1_updated_art_knot_linear$result()
immis_double_log <- brazil_fit2_updated_art_knot_linear$result()
imis_r_logistic <- brazil_fit3_updated_art_knot_linear$result()

brazil_out1_cd4 <- tidy(imis_log_rw) %>% data.frame(model = "R-logRW", .)
brazil_out2_cd4 <- tidy(immis_double_log) %>% data.frame(model = "double logistic", .)
brazil_out3_cd4 <- tidy(imis_r_logistic) %>% data.frame(model = "rlogistic", .)
brazil_out_cd4 <- rbind(brazil_out1_cd4, brazil_out2_cd4, brazil_out3_cd4)



knot_linear_diag_plot <- ggplot(subset(brazil_out_cd4, year %in% 1975:2017), 
                                aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line(size=1.05) + geom_ribbon(alpha=0.2) + 
  facet_wrap(~outcome, scales="free") +
  geom_point(aes(y=lik_data), col="darkred", size=0.5) +
  labs(title = "Spectrum ART numbers, RW compo")


##########################################################################################
## Testing out the spline diagnosis rate model ###########################################
##########################################################################################
devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")

local_test_optim <- fitmod_csavr(brazil, incid_func ="idbllogistic", B0 = 1e4, optfit = TRUE)

brazil$fp$linear_diagnosis <- "spline"

test_optim_1 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0 = 1e4, optfit = TRUE),
                            name = "double_log_spline_diagn")

test_optim_2 <- obj$enqueue(fitmod_csavr(brazil, eppmod = "rlogistic", B0 = 2e4, optfit = TRUE),
                            name = "rlog_spline_diagn")

test_optim_3 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logrw", B0 = 2e4, optfit = TRUE),
                            name = "rw_incid_spline_diagn")

test_optim_4 <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrw", B0 = 2e4, optfit = TRUE),
                            name = "rw_kappa_spline_diagn")

test_optim_5 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logspline", B0 = 1e4, optfit = TRUE),
                            name = "spline_incid_spline_diagn")

test_optim_6 <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrspline", B0 = 1e4, optfit = TRUE),
                            name = "spline_kappa_spline_diagn")


test_optim_1$status()
test_optim_2$status()
test_optim_3$status()
test_optim_4$status()
test_optim_5$status()
test_optim_6$status()

optims <- list(test_optim_1$result(), test_optim_2$result(), test_optim_3$result(), 
               test_optim_4$result(), test_optim_5$result(), test_optim_6$result())

save(optims,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_31_7_2018/6_inicd_meths_spline_diag")

### RW and splines

RW_and_splines <- optims[3:6]
rw_and_splines_rumble <- plot_undiagnosed(RW_and_splines, model_labs = c("incid_RW", "kappa_RW", "incid_spline","kappa_spline"),
                                          knot_linear = FALSE)
rw_and_splines_rumble$combined_plot
rw_and_splines_rumble$rate_plot
ggpubr::annotate_figure(rw_and_splines_rumble$combined_plot,
                        top = ggpubr::text_grob("Spline diagnosis rates, copmarison of RW and spline fitting methods ",
                                                color = "red", size = 14))
### RW only
rw_splines_diag <- optims[3:4]
rw_spline_rumble <- plot_undiagnosed(rw_splines_diag, model_labs = c("incid_RW", "kappa_RW"), knot_linear = FALSE)
ggpubr::annotate_figure(rw_spline_rumble$combined_plot,
                        top = ggpubr::text_grob("Spline diagnosis rates, copmarison of RW fitting methods ",
                                                color = "red", size = 14))

#### Spline only 
splines_splines_diag <- optims[5:6]
spline_spline_rumble <- plot_undiagnosed(splines_splines_diag, model_labs = c("incid_spline", "kappa_spline"), knot_linear = FALSE)
ggpubr::annotate_figure(spline_spline_rumble$combined_plot,
                        top = ggpubr::text_grob("Spline diagnosis rates, copmarison of Spline fitting methods ",
                                                color = "red", size = 14))


optos_rumble <- plot_undiagnosed(optims, model_labs = c("dblog","r_log","incid_RW", "kappa_RW","incid_spline",
                                                        "spline_kappa"),
                                 knot_linear = FALSE)
optos_rumble$combined_plot
optos_rumble$rate_plot
local_test_optim$par[6:14]

################################################################################################
## Re do knot linear ###########################################################################
################################################################################################

brazil$fp$linear_diagnosis <- "knot_linear"

test_optim_1_knot <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0 = 1e4, optfit = TRUE),
                            name = "double_log_knot_diagn")

test_optim_2_knot <- obj$enqueue(fitmod_csavr(brazil, eppmod = "rlogistic", B0 = 1e4, optfit = TRUE),
                            name = "rlog_knot_diagn")

test_optim_3_knot <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logrw", B0 = 2e4, optfit = TRUE),
                            name = "rw_incid_knot_diagn")

test_optim_4_knot <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrw", B0 = 2e4, optfit = TRUE),
                            name = "rw_kappa_knot_diagn")

test_optim_5_knot <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logspline", B0 = 1e4, optfit = TRUE),
                            name = "spline_incid_knot_diagn")

test_optim_6_knot <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrspline", B0 = 1e4, optfit = TRUE),
                            name = "spline_kappa_knot_diagn")


test_optim_1_knot$status()
test_optim_2_knot$status()
test_optim_3_knot$status()
test_optim_4_knot$status()
test_optim_5_knot$status()
test_optim_6_knot$status()

optims_knot <- list(test_optim_1_knot$result(), test_optim_2_knot$result(), test_optim_3_knot$result(),
                    test_optim_4_knot$result(), test_optim_5_knot$result(), test_optim_6_knot$result())

save(optims_knot,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_31_7_2018/6_inicd_meths_KNOT_diag")

### OPTIMS SPLINE AND RW 

sp_and_rw_knot_optims <- optims_knot[3:6]
optims_sp_and_rw <- plot_undiagnosed(sp_and_rw_knot_optims, model_labs = c("incid_RW", "kappa_RW", "incid_spline","kappa_spline"),
                                     knot_linear = T)
optims_sp_and_rw$combined_plot
ggpubr::annotate_figure(optims_sp_and_rw$combined_plot,
                        top = ggpubr::text_grob("Knot linear diagnosis rates, copmarison of RW and spline fitting methods ",
                                                color = "red", size = 14))

## Just RW methods 

rw_knot_optims <- optims_knot[3:4]
rw_knot_rumble <- plot_undiagnosed(rw_knot_optims, model_labs = c("incid_RW", "kappa_RW"), knot_linear = T)
ggpubr::annotate_figure(rw_knot_rumble$combined_plot,
                        top = ggpubr::text_grob("Knot linear diagnosis rates, copmarison of RW fitting methods ",
                                                color = "red", size = 14))

#### Spline only
spline_knot_optims <- optims_knot[5:6]
spline_knot_rumble <- plot_undiagnosed(spline_knot_optims, model_labs = c("incid_spline", "kappa_spline"), knot_linear = T)
ggpubr::annotate_figure(spline_knot_rumble$combined_plot,
                        top = ggpubr::text_grob("Knot linear diagnosis rates, copmarison of SPLINE fitting methods ",
                                                color = "red", size = 14))



optims_sp_and_rw$rate_plot

### 


optos_rumble_knot <- plot_undiagnosed(optims_knot, model_labs = c("dblog","r_log","incid_RW", "kappa_RW","incid_spline",
                                                                  "spline_kappa"), knot_linear = TRUE)
optos_rumble_knot$combined_plot
optos_rumble_knot$rate_plot
local_test_optim$par[6:14]


#####################################################################################
## testing out the neg binom ll model ###############################################
#####################################################################################

brazil$fp$neg_binom <- FALSE

devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")

local_test_optim <- fitmod_csavr(brazil, incid_func ="ilogistic", B0 = 1e4, optfit = TRUE)

brazil$fp$neg_binom <- TRUE

local_test_optim_binom <- fitmod_csavr(brazil, incid_func = "ilogistic", B0 = 1e4, optfit = TRUE)

local_test_optim$par
local_test_optim_binom$par

list_of_res <- list(local_test_optim, local_test_optim_binom)

outty <- plot_undiagnosed(list_of_res, model_labs = c("poisson","neg_binom"))


loccy_test_binom <- fitmod_csavr(brazil, incid_func ="ilogistic", B0 = 1e4, optfit = TRUE)

#####################################################################################
## Re do the spline and knot linear runs from previously as csavr updated ###########
#####################################################################################
brazil$fp$neg_binom <- FALSE
brazil$fp$linear_diagnosis <- "spline"

test_optim_1 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0 = 1e4, optfit = TRUE),
                            name = "double_log_spline_diagn")

test_optim_2 <- obj$enqueue(fitmod_csavr(brazil, eppmod = "rlogistic", B0 = 1e6, optfit = TRUE),
                            name = "rlog_spline_diagn")

test_optim_3 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logrw", B0 = 1e6, optfit = TRUE),
                            name = "rw_incid_spline_diagn")

test_optim_4 <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrw", B0 = 1e6, optfit = TRUE),
                            name = "rw_kappa_spline_diagn")

test_optim_5 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logspline", B0 = 1e6, optfit = TRUE),
                            name = "spline_incid_spline_diagn")

test_optim_6 <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrspline", B0 = 1e4, optfit = TRUE),
                            name = "spline_kappa_spline_diagn")


test_optim_1$status()
test_optim_2$status()
test_optim_3$status()
test_optim_4$status()
test_optim_5$status()
test_optim_6$status()

id_2 <- test_optim_2$id
id_3 <- test_optim_3$id
id_4 <- test_optim_4$id
id_5 <- test_optim_5$id

save(id_2,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/optim_poiss_2")
save(id_3,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/optim_poiss_3")
save(id_4,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/optim_poiss_4")
save(id_5,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/optim_poiss_5")

path_name <- "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/"
clust_ids <- list.files(path_name, full.names = T)
for(i in 1:3)
  load(clust_ids[[i+19]], verbose = T)

test_optim_2 <- obj$task_get(id_2)
test_optim_3 <- obj$task_get(id_3)
test_optim_4 <- obj$task_get(id_4)

optims <- list(test_optim_1$result(), test_optim_2$result(), test_optim_3$result(), 
               test_optim_4$result(), test_optim_5$result(), test_optim_6$result())

save(optims,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_31_7_2018/6_inicd_meths_spline_diag")

###########################################################################################
## Now for the knot linear runs ###########################################################
###########################################################################################
brazil$fp$neg_binom <- FALSE
brazil$fp$linear_diagnosis <- "knot_linear"

test_optim_1_knot <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0 = 1e4, optfit = TRUE),
                                 name = "double_log_knot_diagn")

test_optim_2_knot <- obj$enqueue(fitmod_csavr(brazil, eppmod = "rlogistic", B0 = 1e4, optfit = TRUE),
                                 name = "rlog_knot_diagn")

test_optim_3_knot <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logrw", B0 = 2e4, optfit = TRUE),
                                 name = "rw_incid_knot_diagn")

test_optim_4_knot <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrw", B0 = 2e4, optfit = TRUE),
                                 name = "rw_kappa_knot_diagn")

test_optim_5_knot <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logspline", B0 = 1e4, optfit = TRUE),
                                 name = "spline_incid_knot_diagn")

test_optim_6_knot <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrspline", B0 = 1e4, optfit = TRUE),
                                 name = "spline_kappa_knot_diagn")


test_optim_1_knot$status()
test_optim_2_knot$status()
test_optim_3_knot$status()
test_optim_4_knot$status()
test_optim_5_knot$status()
test_optim_6_knot$status()

optims_knot <- list(test_optim_1_knot$result(), test_optim_2_knot$result(), test_optim_3_knot$result(),
                    test_optim_4_knot$result(), test_optim_5_knot$result(), test_optim_6_knot$result())

save(optims_knot,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_31_7_2018/6_inicd_meths_KNOT_diag")

load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_31_7_2018/6_inicd_meths_KNOT_diag",
     verbose = T)

opto_outo <- plot_undiagnosed(optims_knot, model_labs = c("double_log","rlogistic","incid_rw","kappa_rw",
                                                          "incid_spline", "kappa_spline"))

opto_outo$combined_plot
opto_outo$rate_plot

#### RWs only 

knot_poiss_rws <- optims_knot[]


#############################################################################################
## Runnning the neg binom models on the cluster #############################################
#############################################################################################

brazil$fp$neg_binom <- TRUE
brazil$fp$linear_diagnosis <- "spline"

test_optim_1_neg <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0 = 1e4, optfit = TRUE),
                            name = "double_log_spline_diagn_neg")

test_optim_2_neg <- obj$enqueue(fitmod_csavr(brazil, eppmod = "rlogistic", B0 = 1e6, optfit = TRUE),
                            name = "rlog_spline_diagn_neg")

test_optim_3_neg <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logrw", B0 = 1e6, optfit = TRUE),
                            name = "rw_incid_spline_diagn_neg")

test_optim_4_neg <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrw", B0 = 1e6, optfit = TRUE),
                            name = "rw_kappa_spline_diagn_neg")

test_optim_5_neg <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logspline", B0 = 1e6, optfit = TRUE),
                            name = "spline_incid_spline_diagn_neg")

test_optim_6_neg <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrspline", B0 = 1e4, optfit = TRUE),
                            name = "spline_kappa_spline_diagn_neg")


test_optim_1_neg$status()
test_optim_2_neg$status()
test_optim_3_neg$status()
test_optim_4_neg$status()
test_optim_5_neg$status()
test_optim_6_neg$status()

id_2_neg <- test_optim_2_neg$id
id_3_neg <- test_optim_3_neg$id
id_4_neg <- test_optim_4_neg$id
id_5_neg <- test_optim_5_neg$id

save(id_2_neg,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/optim_NEG_2")
save(id_3_neg,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/optim_NEG_3")
save(id_4_neg,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/optim_NEG_4")
save(id_5_neg,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/optim_NEG_5")
path_name <- "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/"
clust_ids <- list.files(path_name, full.names = T)
for(i in 1:4)
  load(clust_ids[[i+15]], verbose = T)

test_optim_2_neg <- obj$task_get(id_2_neg)
test_optim_3_neg <- obj$task_get(id_3_neg)
test_optim_4_neg <- obj$task_get(id_4_neg)
test_optim_5_neg <- obj$task_get(id_5_neg)

optims_neg <- list(test_optim_1_neg$result(), test_optim_2_neg$result(), test_optim_3_neg$result(), 
               test_optim_4_neg$result(), test_optim_5_neg$result(), test_optim_6_neg$result())

save(optims_neg,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_31_7_2018/6_inicd_meths_spline_diag")

###########################################################################################
## Now for the knot linear runs ###########################################################
###########################################################################################
brazil$fp$neg_binom <- TRUE
brazil$fp$linear_diagnosis <- "knot_linear"

test_optim_1_knot_neg <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0 = 1e4, optfit = TRUE),
                                 name = "double_log_knot_neg_diagn")

test_optim_2_knot_neg <- obj$enqueue(fitmod_csavr(brazil, eppmod = "rlogistic", B0 = 1e4, optfit = TRUE),
                                 name = "rlog_knot_neg_diagn")

test_optim_3_knot_neg <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logrw", B0 = 2e4, optfit = TRUE),
                                 name = "rw_incid_knot_neg_diagn")

test_optim_4_knot_neg <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrw", B0 = 2e4, optfit = TRUE),
                                 name = "rw_kappa_knot_neg_diagn")

test_optim_5_knot_neg <- obj$enqueue(fitmod_csavr(brazil, incid_func = "incid_logspline", B0 = 1e4, optfit = TRUE),
                                 name = "spline_incid_knot_neg_diagn")

test_optim_6_knot_neg <- obj$enqueue(fitmod_csavr(brazil, eppmod = "logrspline", B0 = 1e4, optfit = TRUE),
                                 name = "spline_kappa_knot_neg_diagn")


test_optim_1_knot_neg$status()
test_optim_2_knot_neg$status()
test_optim_3_knot_neg$status()
test_optim_4_knot_neg$status()
test_optim_5_knot_neg$status()
test_optim_6_knot_neg$status()

optims_knot_neg <- list(test_optim_1_knot_neg$result(), test_optim_2_knot_neg$result(), test_optim_3_knot_neg$result(),
                    test_optim_4_knot_neg$result(), test_optim_5_knot_neg$result(), test_optim_6_knot_neg$result())

save(optims_knot_neg,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_31_7_2018/6_inicd_meths_knot_NEG_diag")


