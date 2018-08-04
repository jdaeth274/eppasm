####################################################################################################
## Comparing the different modelling techniques and how they perform through time ##################
####################################################################################################

require(ggplot2)
devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")
require(reshape2)
require(ggpubr)

####################################################################################################
## Now lets load up the different datasets for this run to compare #################################
####################################################################################################

data_path_name <- "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/IMIS_2_8_2018/"
data_files <- list.files(data_path_name, full.names = T)
for(i in 1:length(data_files))
  load(data_files[[i]], verbose = T)

brazil_out1_cd4 <- tidy(results_list_neg_knot[[1]]) %>% data.frame(model = "Kappa RW", .)
brazil_out2_cd4 <- tidy(results_list_neg_knot[[2]]) %>% data.frame(model = "Incid RW", .)
brazil_out3_cd4 <- tidy(results_list_neg_knot[[3]]) %>% data.frame(model = "Kappa Spline", .)
brazil_out4_cd4 <- tidy(results_list_neg_knot[[4]]) %>% data.frame(model = "Incid Spline", .)
brazil_out_cd4 <- rbind(brazil_out1_cd4, brazil_out2_cd4, brazil_out3_cd4, brazil_out4_cd4)

brazil_out1_cd4_spline <- tidy(results_list_neg_spline[[1]]) %>% data.frame(model = "Kappa RW", .)
brazil_out2_cd4_spline <- tidy(results_list_neg_spline[[2]]) %>% data.frame(model = "Incid RW", .)
brazil_out3_cd4_spline <- tidy(results_list_neg_spline[[3]]) %>% data.frame(model = "Kappa Spline", .)
brazil_out4_cd4_spline <- tidy(results_list_neg_spline[[4]]) %>% data.frame(model = "Incid Spline", .)
brazil_out_cd4_spline <- rbind(brazil_out1_cd4_spline, brazil_out2_cd4_spline, brazil_out3_cd4_spline, brazil_out4_cd4_spline)

#####################################################################################################
## Now for our function to compare the different trends through time ################################
#####################################################################################################

rmse_per_time<-function(output_data,metric="prevalence",year_time_series=seq(1970,2021,1)){
  
  if(metric=="prevalence"){
    dataset_to_use <- output_data[output_data$outcome == "Prevalence (15-49)",]
  }
  
  if(metric=="incidence"){
    dataset_to_use <- output_data[output_data$outcome == "Incidence (15-49)",]
    
  }
  
  if(metric=="kappa"){
    dataset_to_use <- output_data[output_data$outcome == "log r(t)",]
  }
  
  
  time_to_test<-seq((year_time_series[1]-1970)+1,(year_time_series[length(year_time_series)]-2021)+52,1)
  
  tot_incid_error<-NULL
  tot_kappa_error<-NULL
  
  for (i in time_to_test[1]:time_to_test[length(time_to_test)]){
    
    year<- (i-1) + 1970
    
    year_data <- dataset_to_use[dataset_to_use$year == year, ]
    
    
    
    error_incid <- year_data[2, 4]  - year_data[4, 4]
    
    error_kappa <- year_data[1, 4] - year_data[3, 4]
    
    
    lower_bound_incid <- year_data[2, 7]  - year_data[4, 7]
    
    lower_bound_kappa <- year_data[1, 7]  - year_data[3, 7]
    
    upper_bound_incid <- year_data[2, 8]  - year_data[4, 8]
    
    upper_bound_kappa <- year_data[1, 8]  - year_data[3, 8]
    
    error_incid <- cbind.data.frame(error_incid,lower_bound_incid,upper_bound_incid, year)
    
    error_kappa <- cbind.data.frame(error_kappa, lower_bound_kappa, upper_bound_kappa, year)
    
    tot_incid_error <- rbind.data.frame(tot_incid_error, error_incid)
    
    tot_kappa_error <- rbind.data.frame(tot_kappa_error, error_kappa)
    
  }
  
  names(tot_incid_error) <- c("mean","low","high","year")
  names(tot_kappa_error) <- c("mean","low","high","year")
  
  return(list(error_incid = tot_incid_error, error_kappa = tot_kappa_error))
  
}

plotter_function_rmse<-function(list_of_rmse_results_dfs,plot_title,colour_by_sample_size=F){
  total_data<-list_of_rmse_results_dfs[[1]][[2]]
  
  for(i in 2:length(list_of_rmse_results_dfs)){
    total_data<-rbind.data.frame(total_data,list_of_rmse_results_dfs[[i]][[2]])
  }
  
  
  names_vector<-c(list_of_rmse_results_dfs[[1]][[2]]$type[1],
                  list_of_rmse_results_dfs[[2]][[2]]$type[1],
                  list_of_rmse_results_dfs[[3]][[2]]$type[1],
                  list_of_rmse_results_dfs[[4]][[2]]$type[1])
  colour_vector<-c("dodgerblue","red","blueviolet","green2")
  names(colour_vector) <- names_vector
  if(colour_by_sample_size == T){
    order_vector <- c(list_of_rmse_results_dfs[[1]][[2]]$order_type[1],
                      list_of_rmse_results_dfs[[5]][[2]]$order_type[1])
  }
  mean_plot<-ggplot(data = total_data,aes(x=time,y=RMSE,group=type))+geom_line(aes(colour=type,linetype=order_type),size=1.2)+
    scale_colour_discrete(name = "Sample size", breaks = c(names_vector), labels = c("100","500","1000","5000"))+
    scale_linetype_discrete(guide=FALSE) +
    labs(x="Time",y="RMSE",title=plot_title)
  
  if(colour_by_sample_size==T){
    mean_plot<-ggplot(data = total_data,aes(x=time,y=RMSE,group=type))+
      geom_line(aes(colour=sample_size,linetype=order_type),size=1.2)+
      #scale_colour_manual("Fitting method",values = colour_vector)+
      labs(x="time",y="RMSE",title=plot_title)+
      scale_colour_discrete(name = "Sample size",
                            breaks = c("100","500","1000","5000"), labels = c("100","500","1000","5000"))+
      scale_linetype_discrete("Penalty order",
                              breaks = c(order_vector), labels = c("First order","Second order"))
    
    
  }
  
  error_plot<-ggplot(data = total_data,aes(x=time,y=median,group=type))+geom_line(aes(colour=type,linetype=order_type),size=1.2)+
    geom_ribbon(aes(x=time,ymin=low,ymax=high,fill=type),alpha=0.36)+
    #scale_colour_manual("Fitting method",values = colour_vector)+
    #scale_fill_manual("Fitting method",values = colour_vector)+
    labs(x="time",y="error",title=plot_title)
  # scale_linetype_discrete("Sample size",
  #                       breaks = c(names_vector), labels = c("100","500","1000","5000"),
  #                       
  # scale_colour_discrete(name = "Sample size",
  #                     breaks = c(names_vector), labels = c("100","500","1000","5000"))
  #                     
  
  return(list(mean_plot=mean_plot,error_plot=error_plot))
  
  
}

##################################################################################################
## Run the functions #############################################################################
##################################################################################################

knot_diagnosis_diff <- rmse_per_time(output_data = brazil_out_cd4)
