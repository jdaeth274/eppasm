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

local_test_optim <- fitmod_csavr(brazil, incid_func ="incid_logrw", B0 = 1e4, optfit = TRUE)

test_optim_1_knot_neg <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0 = 1e4, optfit = TRUE),
                                     name = "double_log_knot_neg_diagn")

brazil$fp$neg_binom <- TRUE

local_test_optim_2 <- fitmod_csavr(brazil, incid_func ="incid_logrw", B0 = 1e4, optfit = TRUE)

###############################################################################
## Plotting functions #########################################################
###############################################################################

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
  if(length(optim_output) == 10){
    lengers <- 1
  }else{
      lengers <- length(optim_output)
    }
  
  for(i in 1:lengers){
    if(length(optim_output) == 10){
      
      list_version <- attributes(optim_output$mod)
      par_vals <- optim_output$par
      
    }else{
      list_version <- attributes(optim_output[[i]]$mod)
      par_vals <- optim_output[[i]]$par
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
      diag_rate <- par_vals[(length(par_vals) - 4) : length(par_vals)]
      knots <- c(1986, 1996, 2000, 2009, 2015)
      
      diagn_trend <- approx(knots, diag_rate, 1970:2021, rule = 2)$y
    }else{
      diag_rate <- par_vals[(length(par_vals) - 6) : length(par_vals)]
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
  
  if(length(optim_output) == 10){
    tot_deaths <- merge(tot_deaths, optim_output$likdat$aidsdeaths, all = TRUE)
    tot_diagnoses <- merge(tot_diagnoses, optim_output$likdat$diagnoses, all = TRUE)
  }else{
    tot_deaths <- merge(tot_deaths, optim_output[[1]]$likdat$aidsdeaths, all = TRUE)
    tot_diagnoses <- merge(tot_diagnoses, optim_output[[1]]$likdat$diagnoses, all = TRUE)
  }
  
  
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


###############################################################################
## 
outputs <- list(local_test_optim, local_test_optim_2)

test_output <- plot_undiagnosed(outputs, model_labs = c("Incid_logrw","neg_binom_log_rw_incid"))
test_output$combined_plot


brazil$fp$linear_diagnosis

brazil$fp$neg_binom <- TRUE

local_test_optim <- fitmod_csavr(brazil, incid_func ="incid_logrw", B0 = 1e4, optfit = TRUE)
