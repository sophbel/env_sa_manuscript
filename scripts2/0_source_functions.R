library(patchwork)
library(INLA)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(purrr)
library(gridExtra)
library(tidyverse)
library(raster)
library(dlnm)
library(lubridate)
library(scales)
library(cowplot)
library(grDevices)
library(magrittr)
library(rlang)

### install packages that are missing
# options(repos = c(CRAN = "https://cloud.r-project.org/"))
# 
# # List of required packages
# required_packages <- c("patchwork", "INLA", "data.table","ggplot2", "dplyr", "tidyr","sf", "purrr", "gridExtra" 
#                        ,"tidyverse", "raster", "argparser")  # Add your packages here
# 
# # Function to check and install packages
# install_if_missing <- function(package) {
#   if (!requireNamespace(package, quietly = TRUE)) {
#     install.packages(package)
#   }
# }
# 
# # Check and install each package
# sapply(required_packages, install_if_missing)
# 
# # Load the packages
# invisible(sapply(required_packages, library, character.only = TRUE))


# rw1 model
hyper2.rw = list(prec = list(prior='pc.prec', param=c(0.3, 0.01))) # medium
# bym2 model
# probability of SD of theta1 > 1 = 0.01
hyper.bym2 = list(theta1 = list(prior="pc.prec", param=c(1, 0.01)),
                  theta2 = list(prior="pc", param=c(0.5, 0.5)))
'%notin%'<-Negate('%in%')

# Define Model--------
inla.mod <- function(form, fam, df, nthreads, config) {
  # Record the start time
  start_time <- Sys.time()
  config = config
  # Your existing inla.mod code
  mod <- INLA::inla(
    as.formula(form$formula),
    family = fam,
    offset=log(population_size/100000),
    control.inla = list(strategy = 'adaptive'),
    control.compute = list(dic = TRUE, config = config, waic=TRUE, cpo=TRUE),
    control.predictor = list(link = 1, compute = TRUE),
    inla.setOption(num.threads = nthreads),
    verbose = FALSE,
    data = df
  )
  
  # add the covariate 
  mod$cov <-form$covs
  # Record the end time
  end_time <- Sys.time()
  
  # Calculate the elapsed time
  elapsed_time <- end_time - start_time
  
  # Print the elapsed time for each model
  cat("Model runtime:", elapsed_time, "\n")
  
  return(mod)
}


### define inla model for a multivariate model -----
inla.mod.multivariate <- function(form, fam, df, stk, nthreads, config) {
  # Record the start time
  start_time <- Sys.time()
  
  if(fam == "binomial"){
    mod <- INLA::inla(
      formula = form$formula,
      family = fam,  # One per outcome
      Ntrials = data_long$month_seq,
      data = INLA::inla.stack.data(stk),  # Extract data from stack
      control.predictor = list(A = INLA::inla.stack.A(stk), compute = TRUE),
      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = config),
      inla.setOption(num.threads = nthreads)
    )
  }else{
  mod <- INLA::inla(
    formula = form$formula,
    family = fam,  # One per outcome
    offset = log(data_long$population_size/100000),
    data = INLA::inla.stack.data(stk),  # Extract data from stack
    control.predictor = list(A = INLA::inla.stack.A(stk), compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = config),
    inla.setOption(num.threads = nthreads)
  )
  }
  
  # Add covariate info
  mod$cov <- form$cov
  
  # Record end time
  end_time <- Sys.time()
  
  # Print runtime
  cat("Model runtime:", end_time - start_time, "\n")
  
  return(mod)
}


# DEFINE INLA MODEL WITH A SEQUENCING OFFSET --------
inla.mod.nooffset <- function(form, fam, df, nthreads, config) {
  # Record the start time
  start_time <- Sys.time()
  config = config
  # Your existing inla.mod code
  mod <- INLA::inla(
    as.formula(form$formula),
    family = fam,
    Ntrials = df$month_seq,
    # offset=combined_offset,
    control.inla = list(strategy = 'adaptive'),
    control.compute = list(dic = TRUE, config = config,waic=TRUE,cpo=TRUE),#, openmp.strategy="omp"
    control.predictor = list(link = 1, compute = TRUE),
    inla.setOption(num.threads = nthreads),
    verbose = FALSE,
    data = df
  )
  
  # add the covariate 
  mod$cov <-form$covs
  # Record the end time
  end_time <- Sys.time()
  
  # Calculate the elapsed time
  elapsed_time <- end_time - start_time
  
  # Print the elapsed time for each model
  cat("Model runtime:", elapsed_time, "\n")
  
  return(mod)
}
#################################################################################
## plots
#################################################################################

## plot Fit
plot.fit <- function(mod, data_mat, title=""){
  covs<-mod$cov
  mod[["summary.fitted.values"]] %>%
    bind_cols(data_mat %>% dplyr::select(date, disease)) %>%
    group_by(date) %>% 
    summarise(cases=sum(disease, na.rm=TRUE),
              mean=sum(mean),
              lower=sum(`0.025quant`),
              upper=sum(`0.975quant`), .groups="drop") %>%
    mutate(date=as.Date(date)) %>% 
    ggplot(aes(x=date)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="#D37295", alpha=0.4)+
    geom_line(aes(y=cases, colour="Observed"),alpha=0.3) +
    geom_line(aes(y=mean, colour="Fitted")) +
    scale_x_date(date_breaks="1 year", date_labels="%Y") +
    labs(x="Years", y="Counts", title = title, subtitle = covs, colour="") +
    scale_color_manual(values = c("Observed"="#499894", "Fitted"="#D37295")) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none",
          plot.title=element_text(hjust=0.5), panel.border = element_blank())
}

## plot effects
### Plot random effects-----
plot_random_effects<-function(base,univ, title_a = "", title_b = "", return = "both"){
  ## extract month effect
  # month_eff_base<-base_mod_list[[1]]$summary.random$id_m
  # month_eff_mod<-univ_mod_list[[1]]$summary.random$id_m
  first_year = 2005
  title_a <- title_a
  title_b <- title_b
  
  month_eff_base<-base$summary.random$id_m
  month_eff_mod<-univ$summary.random$id_m
  

  colnames(month_eff_base)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld")
  cov<-c(univ$cov)
  colnames(month_eff_mod)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld")
  
  ## extract year effect
  year_eff_base<-base$summary.random$id_y
  year_eff_mod<-univ$summary.random$id_y
  
  year_eff_base<-base$summary.random$id_y
  year_eff_mod<-univ$summary.random$id_y
  N_years <- nrow(year_eff_base)
  
  colnames(year_eff_base)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld")
  colnames(year_eff_mod)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld")
  ## number of months and years
  month_n<-nrow(month_eff_base)
  year_n_base<-nrow(year_eff_base)
  year_n_mod<-nrow(year_eff_mod)
  year_eff_base$ID <- first_year + seq(0,(N_years-1),)
  year_eff_mod$ID <- first_year + seq(0,(N_years-1),)
  
  ## add region column
  # reg_month<-c(rep(1,month_n),rep(2,month_n),rep(3,month_n))
  # reg_year_base<-c(rep(1,year_n_base),rep(2,year_n_base),rep(3,year_n_base))
  # reg_year_mod<-c(rep(1,year_n_mod),rep(2,year_n_mod),rep(3,year_n_mod))
  
  # month_eff_base$region<-factor(reg_month)
  # month_eff_mod$region<-factor(reg_month)
  # year_eff_base$region<-factor(reg_year_base)
  # year_eff_mod$region<-factor(reg_year_mod)
  ## linetypes
  lines<-c(title_a="dashed",title_b="solid")
  
  ## variable
  covs<-c(univ$data$cov)
  ## 
  a<-ggplot()+
    
    ## with region
    # geom_line(data=month_eff_base,aes(x=ID,y=median,group=region,color=region,linetype="Base Model"))+
    # geom_line(data=month_eff_mod,aes(x=ID,y=median,group=region,color=region,linetype="Indicator"))+
    
    ##without region
    geom_line(data=month_eff_base,aes(x=ID,y=median,linetype=title_a),color="#311B83")+
    geom_line(data=month_eff_mod,aes(x=ID,y=median,linetype=title_b))+
    geom_hline(yintercept=0,linetype="dashed",color="darkred",alpha=0.8)+
    theme_classic()+
    labs(x = "Month",
         y = "Median",
         linetype = "") +
    # ylim(-1,1)+
    # ylim(-0.4,0.4)+
    ggtitle("Monthly Random Effect") +
    labs(subtitle=covs)
  b<-ggplot()+
    geom_hline(yintercept=0,linetype="dashed",color="darkred",alpha=0.8)+
    geom_line(data=year_eff_base,aes(x=ID,y=median,linetype=title_a),color="#311B83")+
    geom_line(data=year_eff_mod,aes(x=ID,y=median,linetype=title_b))+
    theme_classic()+
    labs(x = "Year",
         y = "Median",
         linetype = "") +
    # ylim(-2,2)+
    ggtitle("Yearly Random Effect")+    
    labs(subtitle=covs)
  
  # return(a)
  if(return=="both"){
  return(grid.arrange(a,b, ncol=2, nrow=1))
  }
  if(return=="year"){
    return(b)
  }
  if(return=="month"){
    return(a)
  }
}


## plot random effect with provincial replicate for the seasonal effect
### Plot random effects-----
plot_random_effect_provRep<-function(base,univ, title_a = "", title_b = "", type = "annual"){
  
  title_a <- title_a
  title_b <- title_b
  
  if(type=="seasonal"){
  month_eff_base<-base$summary.random$id_m
  month_eff_mod<-univ$summary.random$id_m
  
  colnames(month_eff_base)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld")
  cov<-c(univ$cov)
  colnames(month_eff_mod)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld")
  ## number of months and years
  month_n<-nrow(month_eff_base)

  ## add region column
  reg_month <- rep(c("Eastern Cape", "Free State", "Gauteng", 
                     "KwaZulu-Natal", "Limpopo", "Mpumalanga", 
                     "North West", "Northern Cape", "Western Cape"), each = month_n)
  month_eff_base$region<-factor(reg_month)

  ## linetypes
  lines<-c(title_a="dashed",title_b="solid")
  
  ## variable
  covs<-c(univ$cov)
  ## 
  a<-ggplot()+
    # geom_line(data=month_eff_base,aes(x=ID,y=median,group=region,color=region,linetype="Base Model"))+
    # geom_line(data=month_eff_mod,aes(x=ID,y=median,group=region,color=region,linetype="Indicator"))+
    geom_line(data=month_eff_base,aes(x=ID,y=median,linetype=title_a, group=region))+
    geom_line(data=month_eff_mod,aes(x=ID,y=median,linetype=title_b))+
    geom_hline(yintercept=0,color="darkred",alpha=0.8)+
    # scale_linetype_manual(values = lines)+
    theme_classic()+
    theme( axis.text = element_text(size=15), axis.title = element_text(size=15))+
    labs(x = "Month",
         y = "Median",
         linetype = "",
         color="Region") +
    ylim(-1,1)+
    # ylim(-0.4,0.4)+
    # ggtitle("Monthly Random Effect") +
    labs(subtitle=covs)
  }
  if(type=="annual"){
    month_eff_base<-base$summary.random$id_y
    month_eff_mod<-univ$summary.random$id_y
    
    colnames(month_eff_base)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld")
    cov<-c(univ$cov)
    colnames(month_eff_mod)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld")
    ## number of months and years
    month_n<-length(unique(month_eff_base$ID))
    
    ## add region column
    reg_month <- rep(c("Eastern Cape", "Free State", "Gauteng", 
                       "KwaZulu-Natal", "Limpopo", "Mpumalanga", 
                       "North West", "Northern Cape", "Western Cape"), each = month_n)
    month_eff_base$region<-factor(reg_month)
    month_eff_mod$region<-factor(reg_month)
    ## linetypes
    lines<-c(title_a="dashed",title_b="solid")
    
    ## variable
    covs<-c(univ$cov)
    ## 
    a<-ggplot()+
      # geom_line(data=month_eff_base,aes(x=ID,y=median,group=region,color=region,linetype="Base Model"))+
      # geom_line(data=month_eff_mod,aes(x=ID,y=median,group=region,color=region,linetype="Indicator"))+
      geom_line(data=month_eff_base,aes(x=ID,y=median,linetype=title_a, group=region))+
      geom_line(data=month_eff_mod,aes(x=ID,y=median,linetype=title_b, group=region, color =region))+
      geom_hline(yintercept=0,color="darkred",alpha=0.8)+
      # scale_linetype_manual(values = lines)+
      theme_classic()+
      theme( axis.text = element_text(size=15), axis.title = element_text(size=15))+
      labs(x = "Year",
           y = "Median",
           linetype = "",
           color="Region") +
      ylim(-1,1)+
      # ylim(-0.4,0.4)+
      # ggtitle("Monthly Random Effect") +
      labs(subtitle=covs)+
      facet_wrap(region ~.)
  }
  
  return(a)
}

## plot nonlinear effects
plot.nonlinear.ef <- function(df,cov){
  cov<-cov
  df %>%
    ggplot() +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="red",alpha=0.8)+
    geom_line(aes(y=mean, x=ID), colour="#006666") +
    geom_ribbon(aes(ymin=lowerCI, ymax=upperCI, x=ID),
                alpha=0.4, fill="#17becf") +
    labs(x=cov, y="Effect") +
    # ggtitle(paste0(spvec[sp],", ",cov))+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12), plot.title=element_text(hjust=0.5))
}
## plot spatial effects
plot_spatial_effects<-function(mod,shp, N, structured = TRUE, title_a = "Spatial Random Effects"){
  mod<-mod
  shp_a1<-shp
  if(is.null(shp_a1$areaid)){
    if(N==52){
    shp_a1$areaid <- shp_a1$GID_2
    }else{
      shp_a1$areaid <- shp_a1$GID_1
      
    }
  }
  shp_a1 <- shp_a1[order(shp_a1$areaid), ]
  cov<-mod$cov
  unstructured_effects <- mod$summary.random$id_u$mean[1:N]
  structured_effects <- mod$summary.random$id_u$mean[(N+1):(N*2)]
  shp_a1$structured_effects<-structured_effects
  shp_a1$unstructured_effects<-unstructured_effects
  
  a<-ggplot(data = shp_a1) +
    geom_sf(aes(fill = unstructured_effects),color=NA) +
    # scale_fill_viridis_c() +
    scale_fill_gradient2(low="blue",mid="white",high="red",midpoint = 0)+
    theme_minimal() +
    # theme(legend.position = "none")+
    # labs(title = paste0("Spatial Random Effects: ",spvec[sp],", ",cov),
    labs(title = title_a,
         fill = "Random Effect",
         subtitle = "Full Effects")
  b<-ggplot(data = shp_a1) +
    geom_sf(aes(fill = structured_effects),color=NA) +
    scale_fill_gradient2(low="blue",mid="white",high="red",midpoint = 0)+
    theme_minimal() +
    labs(title = "",
         fill = "Random Effect",
         subtitle = "Structured Effects")
  # return(a)
  if(structured == TRUE){
    return(a+b)
  }else{
    return(a)
  }
}


## plot spatial effects Difference
plot_spatial_effects_diff<-function(mod, base_mod, shp, N, title_a = "", title_b = "", structured = FALSE){
  mod<-mod
  base_mod<-base_mod
  shp_a1<-shp
  cov<-mod$data$cov
  cov <- gsub("_lag0","",cov)
  unstructured_effects <- mod$summary.random$id_u$mean[1:N]
  unstructured_effects_base <- base_mod$summary.random$id_u$mean[1:N]
  unstructured_diff <- unstructured_effects_base-unstructured_effects
  
  structured_effects <- mod$summary.random$id_u$mean[(N+1):(N*2)]
  structured_effects_base <- base_mod$summary.random$id_u$mean[(N+1):(N*2)]
  structured_diff <- structured_effects_base -structured_effects
  
  shp_a1$structured_effects<-structured_diff
  shp_a1$unstructured_effects<-unstructured_diff
  
  if(structured == FALSE){
    a<-ggplot(data = shp_a1) +
      geom_sf(aes(fill = unstructured_effects),color=NA) +
      # scale_fill_viridis_c() +
      scale_fill_gradient2(low="red",mid="white",high="blue",midpoint = 0, limits = c(-0.5,0.5))+
      theme_minimal() +
      theme(legend.position = "right")+
      labs(title = cov,
           fill = "")
           # caption = paste0("unstrucutured effects: ",title_b,"-", title_a))
  }else{
    a<-ggplot(data = shp_a1) +
      geom_sf(aes(fill = structured_effects),color=NA) +
      # scale_fill_viridis_c() +
      scale_fill_gradient2(low="red",mid="white",high="blue",midpoint = 0 , limits = c(-0.5,0.5))+
      theme_minimal() +
      theme(legend.position = "right")+
      labs(title = cov,
           fill = "")
           # caption = paste0("structured effects: ", title_b,"-", title_a))
  }
  # labs(title = cov,
  #        fill = "Spatial Effect Difference",
  #        caption = "Blue: Model Decreased Effect; Red: Model Increased Effect")
  return(a)
}

#### plot difference to the best WAIC
plot_diff_waic<-function(summary_lags,title=""){
      title=title
      library(viridisLite)
      # Generate a viridis color scale
      viridis_colors <- viridis(256)
      # Assuming yellow is in the upper part, we can take colors up to 75% of the scale
      new_viridis_colors <- viridis_colors[quantile(1:192)]  # 75% of 256 is 192
      # Define a function to create a truncated viridis palette
      scale_color_truncated_viridis_d <- function(...) {
        discrete_scale("colour", "viridis", palette = truncated_viridis, ...)
      }
      ### colored by WAIC
      if(length(table(summary_lags$variable))>0){
      summary_lags$cov2 <- summary_lags$variable}else{
      summary_lags$cov2<-gsub("_", " ", summary_lags$cov)}
      # summary_lags<-summary_lags[-which(summary_lags$cov%in%c(socio_covs,"tasrange","sfcWind","tasmean3","tasmax3","tasmin3")),]
      qs<-quantile(summary_lags$diff_waic_lags)
      # create a new table to specify those which are red
      tmp2<-summary_lags
      tmp2$diff_waic_lags<-ifelse(tmp2$diff_waic_lags>0,NA,0)
      
      plot2<-ggplot(subset(summary_lags,summary_lags$diff_waic_lags>0))+
        facet_grid(cov2~.,scales="free", labeller = label_wrap_gen())+
        # facet_grid(.~cov2,scales="free", labeller = label_wrap_gen())+
        
        geom_point(aes(color=diff_waic_lags,x=lag,y=diff_waic_lags),size=3)+
        scale_color_viridis_b(breaks=qs[2:length(qs)],name="Difference to Best WAIC",labels=c("Better Model","","Worse Model",""), option="magma",begin=0.1,end=0.7)+
        # geom_point(data=subset(tmp2,tmp2$diff_waic_lags==0), aes(color="red",x=lag,y=diff_waic_lags),size=3)+
        geom_point(data=tmp2,aes(x=lag,y=diff_waic_lags),size=6,color="black",fill="red",shape="triangle")+
        theme_bw()+
        theme(strip.text.y=element_text(angle=90,size=9),
              strip.text.y.right = element_text(angle=360,size=11),
              strip.background = element_rect(fill="white"),
              axis.text = element_text(size=12),
              legend.title = element_text(size=12))+
        # ylim(-10,654)+
        ggtitle(title)+
        ylab("Difference to the best lag (WAIC)")+
        xlab("Lag (months)")
      return(plot2)

}

### plot the WAIC
plot_waic<-function(summary_lags){
  library(viridisLite)
  # Generate a viridis color scale
  viridis_colors <- viridis(256)
  # Assuming yellow is in the upper part, we can take colors up to 75% of the scale
  new_viridis_colors <- viridis_colors[quantile(1:192)]  # 75% of 256 is 192
  # Define a function to create a truncated viridis palette
  scale_color_truncated_viridis_d <- function(...) {
    discrete_scale("colour", "viridis", palette = truncated_viridis, ...)
  }
  ### colored by WAIC
  summary_lags$cov2<-gsub("_", " ", summary_lags$cov)
  qs<-quantile(summary_lags$waic)
  # create a new table to specify those which are red
  tmp2<-summary_lags
  tmp2$waic<-ifelse(tmp2$diff_waic_lags>0,NA,tmp2$waic)
  
  plot2<-ggplot(subset(summary_lags,summary_lags$diff_waic_lags>0))+
    facet_grid(.~cov2,scales="free", labeller = label_wrap_gen())+
    geom_point(aes(color=waic,x=lag,y=waic),size=3)+
    # scale_color_viridis_b(breaks=qs[2:length(qs)],name="WAIC",labels=c("Better Model","","Worse Model",""), option="magma",begin=0.1,end=0.7)+
    geom_point(data=tmp2,aes(x=lag,y=waic),size=3,color="black",fill="red",shape="triangle")+
    theme_bw()+
    theme(strip.text.y=element_text(angle=90,size=9),
          strip.text.y.right = element_text(angle=360,size=11),
          strip.background = element_rect(fill="white"),
          axis.text = element_text(size=12),
          legend.title = element_blank())+
    # ylim(-10,654)+
    ylab("WAIC")+
    xlab("Lag (months)")
  return(plot2)
  
}

### evaluation
# bayesian r squared function
rsq <- function(model, null,num_outcomes){
  dev_model <- model$deviance_dic 
  dev_null <- null[["dic"]][["deviance.mean"]]
  if(num_outcomes!=1){
    n <- model$N/num_outcomes
  }else{
  n <- model$N
  }
  rsq <- 1 - exp((-2 / n) * ((dev_model / -2) - (dev_null / -2)))
  return(round(rsq, 5))
}

mae.func <- function(mod, data) {
  # type options 
  # mean, 0.5quant, mode (see summary(mod[["summary.fitted.values"]]))
  ## mean absolute error
  # Find all of your absolute errors, xi – x.
  errors.mean <- mod[["summary.fitted.values"]][["mean"]] - data[["disease"]]
  # Add them all up.
  ae <- sum(abs(errors.mean),na.rm=T)
  # Divide by the number of errors.
  mean.ae <- ae / length(errors.mean)
  
  ## median absolute error
  # Find all of your absolute errors, xi – x.
  errors.med <- mod[["summary.fitted.values"]][["0.5quant"]] - data[["disease"]]
  # Add them all up.
  ae.med <- sum(abs(errors.med),na.rm=T)
  # Divide by the number of errors.
  med.ae <- ae.med / length(errors.med)
  
  mae_vec <- c(mean.ae, med.ae)
  names(mae_vec) <- c("mae", "med_ae")
  return(mae_vec)
}

eval.mod <- function(mod, data) {
  print("Running eval.mod")
  
  fixed_mat <- mod$summary.fixed 
  if((nrow(fixed_mat) == 1) ){
    fixed_mat <- mod$summary.fixed
  }else{
    fixed_mat <- fixed_mat[-1,]}
  
  # extract the mean and median absolute error
  maes <- mae.func(mod, data)
  
  # extract cov names
  if(length(mod$cov)==1){
    cov<-mod$cov}else{
      cov <- paste(mod$cov,collapse="__")
    }
  
  # deviance summary from dic
  
  # collect all info in dataframe
  df <- data.frame(
    cov = cov, 
    waic = mod$waic$waic,
    mae = maes[1],
    # cpo=mod$cpo$cpo
    cpo = -mean(log(mod[["cpo"]][["cpo"]]), na.rm=TRUE),
    dic = mod$dic$dic,
    deviance_dic = mod$dic$deviance.mean,
    N = nrow(mod$summary.fitted.values)
  )
  
  return(df)
}

# Function 1: Overall Model Summary
eval.binomial.multivariate <- function(mod, data_tmp) {
  print("Running eval_model_summary")
  
  covariates <- mod$cov
  summary_df <- data.frame(covariate = covariates)
  
  # Extract global metrics
  summary_df$waic <- mod$waic$waic
  summary_df$cpo <- -mean(log(mod$cpo$cpo), na.rm = TRUE)
  summary_df$dic <- mod$dic$dic
  summary_df$deviance_dic <- mod$dic$deviance.mean
  if((is.null(data_tmp$clusters))){
  summary_df$N <- nrow(data_tmp)/length(table(data_tmp$outcome))
  summary_df$N_outcomes <- length(table(data_tmp$outcome))
  }else{
    summary_df$N <- nrow(data_tmp)/length(table(data_tmp$clusters))
    summary_df$N_outcomes <- length(table(data_tmp$clusters))
  }
  
  return(summary_df)
}
# Function 2: Extract per-GPSC effects and MAE
extract_outcome_effects <- function(mod, data_tmp, family = "binomial") {
  print("Running extract_gpsc_effects")
  
  outcomes <- unique(data_tmp$outcome)
  covs <- mod$cov
  
  # Initialize result storage
  gpsc_results <- data.frame()
  
  for (outcome in outcomes) {
    # Extract fixed effects for the current outcome
    fixed_mat <- mod$summary.fixed[grep(outcome, rownames(mod$summary.fixed)), ]
    colnames(fixed_mat) <- c("mean", "sd", "lowerCI", "median", "upperCI", "mode", "kld")
    
    ## Extract outcome-specific data
    outcome_idx <- which(data_tmp$outcome == outcome)
    outcome_data <- data_tmp[outcome_idx, ]    
    # Extract fitted probabilities and calculate predicted counts
    fitted_probs <- mod$summary.fitted.values$mean[outcome_idx]
    predicted_counts <- fitted_probs * outcome_data$month_seq
    
    # Calculate MAE
    observed_counts <- outcome_data$count
    mae <- mean(abs(predicted_counts - observed_counts), na.rm = TRUE)
    
    # Store results
    gpsc_row <- data.frame(
      GPSC = outcome,
      effect_median = fixed_mat$median,
      effect_lowerCI = fixed_mat$lowerCI,
      effect_upperCI = fixed_mat$upperCI,
      mae = mae,
      cov = covs
    )
    
    gpsc_results <- rbind(gpsc_results, gpsc_row)
  }
  
  return(gpsc_results)
}


### for a low number of outcomes
eval.mod.multivariate <- function(mod, data_tmp, effect=TRUE, family="nbinomial") {
  print("Running eval.mod")
  
  # Get unique outcomes and covariates
  outcomes <- unique(data_tmp$outcome)
  covariates <- mod$cov
  
  # Initialize an empty data frame for storing results
  result_df <- data.frame(covariate = covariates)
  
  # Loop over each outcome to extract effects and metrics
  for (outcome in outcomes) {
    
    # Extract fixed effects for the current outcome
    fixed_mat <- mod$summary.fixed[grep(outcome, rownames(mod$summary.fixed)), ]
    colnames(fixed_mat) <- c("mean", "sd", "lowerCI", "median", "upperCI", "mode", "kld")
    
    ## extract outcome data
    outcome_data <- data_tmp[which(data_tmp$outcome==outcome),]
    
    ## Calculate the Mean Absolute Error (MAE) for each
    if(family=="binomial"){
      print("WRONG FUNCTION: RUN 'eval.binomial.multivariate() and extract_outcome_effects()")
      }else{
    errors.mean <- mod[["summary.fitted.values"]][["mean"]][which(data_tmp$outcome == outcome)] -
      outcome_data[["count"]][which(outcome_data$outcome == outcome)]
    mae <- mean(abs(errors.mean), na.rm = TRUE)
    result_df[[paste0("mae_", outcome)]] <- mae
    }
    
    # Store extracted metrics and effects in the result data frame
    if(effect==TRUE){
      result_df[[paste0("effect_med_", outcome)]] <- fixed_mat$median
      result_df[[paste0("effect_lowerCI_", outcome)]] <- fixed_mat$lowerCI
      result_df[[paste0("effect_upperCI_", outcome)]] <- fixed_mat$upperCI
      # result_df[[paste0("mae_", outcome)]] <- mae
    }
    # Add metrics only once (same across all covariates)
    if (!"waic" %in% colnames(result_df)) {
      result_df$waic <- mod$waic$waic
      # result_df$mae <- mae
      result_df$cpo <- -mean(log(mod[["cpo"]][["cpo"]]), na.rm = TRUE)
      result_df$dic <- mod$dic$dic
      result_df$deviance_dic <- mod$dic$deviance.mean
      result_df$N <- length(errors.mean)
    }
  }
  
  return(result_df)
}

eval.ranks<- function(mod_out){
  ## create a ranks core for the mae's and dic where 1 is the best model (lowest number) for each. 
  # mod_out is the output from function eval.mod and includes the column names "dic", "mean_ae", and "med_ae"
  ## create a ranks core for the mae's and dic where 1 is the best model (lowest number) for each. 
  # mod_out is the output from function eval.mod and includes the column names "dic", "mean_ae", and "med_ae"
  mod_out<-mod_out[with(mod_out, order(mae)),]
  mod_out$mae_score<-1:nrow(mod_out)
  # mod_out<-mod_out[with(mod_out, order(med_ae)),]
  # mod_out$med_ae_score<-1:nrow(mod_out)
  mod_out<-mod_out[with(mod_out,order(cpo)),]
  mod_out$cpo_score<-1:nrow(mod_out)
  mod_out<-mod_out[with(mod_out,order(waic)),]
  mod_out$waic_score<-1:nrow(mod_out)
  
  
  tmpscore<-mod_out$mae_score+mod_out$waic_score+mod_out$cpo_score
  mod_out<-mod_out[order(tmpscore),]
  mod_out$score<-1:nrow(mod_out)
  return(mod_out)
}


####################################################
### SET IDS
####################################################
## set function to create the id_indexes
re_id2 <- function(data, type){
  type <- type
  if(type =="week"){
    data$date<- data$week
   }
  if(type=="month"){
    data$date <- data$month
  }
  ## year and month
  data$year <- year(data$date)
  data$id_y <- as.numeric(factor(year(data$date), levels = seq(min(year(data$date)), max(year(data$date)))))
  data$month <- month(data$date)
  data$id_m <- as.numeric(factor(month(data$date), levels = seq(min(month(data$date)), max(month(data$date)))))
  
  ## check monthly
  tbl <- table(data$id_m, data$month)
  all_off_diag_zero <- all(tbl[row(tbl) != col(tbl)] == 0)
  if (all_off_diag_zero) {
    print("Only diagonal has numbers monthly")
  } else {
    print("There are non-zero values outside the diagonal")
  }
  ## check  yearly
  tbl <- table(data$id_y, data$year)
  all_off_diag_zero <- all(tbl[row(tbl) != col(tbl)] == 0)
  if (all_off_diag_zero) {
    print("Only diagonal has numbers yearly")
  } else {
    print("There are non-zero values outside the diagonal")
  }
  return(data)
}

################################################################################
#### Lag Dataframes
###############################################################################
# create lags for 1-3 months for all covariates and for each admin area separately
lag_dataframe<- function(datmat, cov_names, lag){
  print("lag data function")
  lag <- lag
  print(paste0("lag of: ",lag))
  cov_names=cov_names
  print(cov_names)
  if("date"%notin%colnames(datmat)){
    print("change week to date")
    datmat$date<-datmat$week
  }
  # Sort the data by id and time to ensure the correct order (to be on the safe side)
  datmat <- datmat %>%
    arrange(id_u, date, year, month)
  
  
  # Split the dataset by area
  data_split <- datmat %>%
    group_split(id_u)
  
  # Function to create lagged variables for a specified lag and variable
  create_lagged_data <- function(df, variable, lag) {
    mutate(df, !!paste0(variable, "_lag", lag) := lag(!!sym(variable), lag))
  }
  
  # Apply the lag function to each split dataset for each variable and lags specified
  if(lag==6){
  lagged_data <- map(data_split, ~{
    data_for_variable <- .x
    for (variable in cov_names) {
      data_for_variable <- create_lagged_data(data_for_variable, variable, 1)
      data_for_variable <- create_lagged_data(data_for_variable, variable, 2)
      data_for_variable <- create_lagged_data(data_for_variable, variable, 3)
      data_for_variable <- create_lagged_data(data_for_variable, variable, 4)
      data_for_variable <- create_lagged_data(data_for_variable, variable, 5)
      data_for_variable <- create_lagged_data(data_for_variable, variable, 6)
      
    }
    data_for_variable  # Return the dataframe with lagged variables for all specified variables
  } )
  }
  if(lag==12){
    lagged_data <- map(data_split, ~{
      data_for_variable <- .x
      for (variable in cov_names) {
        data_for_variable <- create_lagged_data(data_for_variable, variable, 1)
        data_for_variable <- create_lagged_data(data_for_variable, variable, 2)
        data_for_variable <- create_lagged_data(data_for_variable, variable, 3)
        data_for_variable <- create_lagged_data(data_for_variable, variable, 4)
        data_for_variable <- create_lagged_data(data_for_variable, variable, 5)
        data_for_variable <- create_lagged_data(data_for_variable, variable, 6)
        data_for_variable <- create_lagged_data(data_for_variable, variable, 7)
        data_for_variable <- create_lagged_data(data_for_variable, variable, 8)
        data_for_variable <- create_lagged_data(data_for_variable, variable, 9)
        data_for_variable <- create_lagged_data(data_for_variable, variable, 10)
        data_for_variable <- create_lagged_data(data_for_variable, variable, 11)
        data_for_variable <- create_lagged_data(data_for_variable, variable, 12)
        
      }
      data_for_variable  # Return the dataframe with lagged variables for all specified variables
    } )
  }
  
  # Combine the lagged datasets back together
  data_mat <- bind_rows(lagged_data, .id = "id_u")
  data_mat$id_u<-as.numeric(data_mat$id_u) # reformat id_u after rowbind
  return(data_mat)
}


# Construct formulas -----
construct_formula <- function(vars, mod) {
  formulas_list <- lapply(vars, function(var) {
    if (grepl("_nl$", var)) {
      var_base <- sub("_nl$", "", var)  # Remove the "_nl" suffix
      transformed_var <- paste0("f(inla.group(", var_base, ", n = 8, method='cut'), model='rw2')")
      var_name <- var_base
    } else {
      transformed_var <- var
      var_name <- var
    }
    
    if(mod=="mod7"){
    formula <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                             'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                             'f(id_y, model = "iid", replicate = id_prov)',
                             'vaccination_period','population_density',
                  transformed_var
                  ),"disease")
    }
    
    list(formula = formula, covs = var_name)
  })
  
}


############### PCV ART COMPARISON FUNCTIONS ###################################
############### construct a formula to test perturbations with province replication and without #########

# Function to build a base formula list
perturbation_formula_list <- function(response_var = "vts", replicate = FALSE) {
  
  if(replicate == FALSE){
  # Define the base components used in all models
  base_components <- c(
    'f(id_u, model = "bym2", graph = g, scale.model = TRUE, adjust.for.con.comp = TRUE)',
    'f(id_m, model = "rw2", cyclic = TRUE, scale.model = TRUE, constr = TRUE)',
    'f(id_y, model = "iid")')
    # Define the covariate sets and their descriptions
    covariates <- list(
      list(var = "ART_coverage",
           desc = "spatial, seasonal, and annual accounting for ART coverage"),
      list(var = "PCV_coverage",
           desc = "spatial, seasonal, and annual accounting for PCV coverage"),
      list(var = "vaccination_period",
           desc = "spatial, seasonal, and annual accounting for vaccination_period"),
      list(var = "post_vaccination_2009",
           desc = "spatial, seasonal, and annual accounting for 2009 vaccine"),
      list(var = "ART_coverage_national",
           desc = "spatial, seasonal, and annual accounting for ART Nationally"),
      list(var = c("vaccination_period", "ART_coverage"),
           desc = "spatial, seasonal, and annual accounting for vaccination_period and ART coverage")
    )
  
  }
  if(replicate ==TRUE){
    # Define the base components used in all models
    base_components <- c(
      'f(id_u, model = "bym2", graph = g, scale.model = TRUE, adjust.for.con.comp = TRUE)',
      'f(id_m, model = "rw2", cyclic = TRUE, scale.model = TRUE, constr = TRUE)',
      'f(id_y, model = "iid", replicate = id_prov)'
    )
    # Define the covariate sets and their descriptions
    covariates <- list(
      list(var = "ART_coverage",
           desc = "spatial, seasonal, and annual (replicate) accounting for ART coverage"),
      list(var = "PCV_coverage",
           desc = "spatial, seasonal, and annual (replicate) accounting for PCV coverage"),
      list(var = "vaccination_period",
           desc = "spatial, seasonal, and annual (replicate) accounting for vaccination_period"),
      list(var = "post_vaccination_2009",
           desc = "spatial, seasonal, and annual (replicate) accounting for 2009 vaccine"),
      list(var = "ART_coverage_national",
           desc = "spatial, seasonal, and annual (replicate) accounting for ART Nationally"),
      list(var = c("vaccination_period", "ART_coverage"),
           desc = "spatial, seasonal, and annual (replicate) accounting for vaccination_period and ART coverage")
    )
  }
  
  # Initialize the output list
  out_list <- list()
  
  # Loop over each covariate set and build the model
  for (i in seq_along(covariates)) {
    base_form <- list()
    form <- reformulate(
      c(1, base_components, covariates[[i]]$var),
      response = response_var
    )
    base_form$formula <- as.formula(form)
    base_form$covs <- paste0(response_var, " ", covariates[[i]]$desc)
    out_list[[i]] <- base_form
  }
    # assign(list_name, out_list, envir = .GlobalEnv)
  
  # Also return it for convenience
  return(out_list)
}
### unscale scaled variables
unscale<-function(list,mean_1,sd_1){
  new_ID<-list$ID
  new_ID<-(list$ID*sd_1)+mean_1
  if(length(new_ID)==1){
    new_ID
  }else{
    new_ID
  }
  return(new_ID)
}


##### MODEL EVALUATION #####
evaluate_model_list <- function(model_list,
                                data,
                                base_intercept,
                                num_outcomes = 1) {
  # 1️ Evaluate all models in the list
  mod_out <- lapply(model_list, function(mod) eval.mod(mod, data))
  
  # 2️ Combine all model evaluation metrics into one data frame
  mod_out <- do.call(rbind, mod_out)
  
  # 3️ Evaluate and rank the models across performance metrics
  mod_out <- eval.ranks(mod_out)
  
  # 4️ Clean and standardize model and covariate names
  mod_out$model <- rownames(mod_out)
  mod_out$cov <- gsub(" ", "_", mod_out$cov)
  mod_out$cov <- gsub(",", "", mod_out$cov)
  
  # 5️ Compute R² for each model
  mod_out$rsq <- apply(
    mod_out,
    1,
    function(x) rsq(
      mod_out[mod_out$cov == x["cov"], ],
      base_intercept,
      num_outcomes = num_outcomes
    )
  )
  
  # 6️ Extract numeric model index and base name
  mod_out$num <- gsub("mae", "", mod_out$model)
  mod_out$num[which(mod_out$num == "")] <- 0
  mod_out$num <- as.numeric(mod_out$num)
  mod_out$base <- paste0("base", "_none")
  
  #  Return the clean evaluated data frame
  return(mod_out)
}


############### END PCV ART COMPARISON FUNCTIONS ########3



### extract lag and create lag column
extract_best_lag <- function(data, lagN=6, type){
  # This function extracts lags from the 'data' dataframe and creates a new lag column.
  # It calculates the WAIC difference relative to the best value for each 'cov' group.
  
  # Arguments:
  # - data: A data frame containing lagged covariate values with a column named 'cov'.
  # - lagN: Maximum lag to consider (default is 6).
  # - type: Specifies the processing type. Options are:
  #     * "cov" - applies 'eval.ranks' to evaluate ranks within each lag group.
  #     * "all" - calculates WAIC differences without applying 'eval.ranks'.
  # extract the lags and create new column
  fin<-NULL
  for(i in 0:lagN){
    uno<-data[grep(paste0("lag",i),data$cov),]
    if(nrow(uno)==0){next}
    uno$cov<-unlist(strsplit(uno$cov,paste0("_lag",i)))
    uno$lag<-i
    fin<-rbind(fin,uno)
  }
  tmp1<-fin

  
  ## calculate difference to best
  covs<-unique(tmp1$cov)
  summary_lags<-NULL
  out_lags_sum<-NULL
  for(j in 1:length(covs)){
    df2<-subset(tmp1,tmp1$cov==covs[j])
    if(nrow(df2)>0){
      if(type=="cov"){
      df2<-eval.ranks(df2)
      }
      df2$diff_waic_lags<-df2$waic-min(df2$waic)
      out_lags_sum<-rbind(out_lags_sum,df2)
      
    }
    
    # summary_lags<-rbind(summary_lags,out_lags_sum)
    summary_lags<-out_lags_sum

  }
  return(summary_lags)
}

#### save model outputs
extract_model_results <- function(univ_model_list, df, unscaled_data, type = "model") {
  model_out <- list()
  
  for (k in seq_along(univ_model_list)) {
    if(type=="model"){
    cov <- univ_model_list[[k]]$cov
    }
    if(type=="summary"){
      cov <- univ_model_list[[k]]$data$cov
      }
    
    if(length(strsplit(cov,"_")[[1]])>2){
      cov_scaled <- sapply(strsplit(cov, "_"), function(x) paste0(x[3],"_",x[4]))
      }
    
    
    # Goodness of fit
    print("running eval.mod")
    gof <- eval.mod(univ_model_list[[k]], df)
    # gof$rsq <- rsq(gof,gpsc_seasonal_list[[k]],1)
    gof$rsq <- rsq(gof,int_mod,1)
    
    # Fixed effects
    fixed <- univ_model_list[[k]]$summary.fixed
    
    # Random effects
    mod_spat <- univ_model_list[[k]]$summary.random$id_u
    mod_spat$model_n <- k
    
    mod_y <- univ_model_list[[k]]$summary.random$id_y
    mod_y$model_n <- k
    
    mod_m <- univ_model_list[[k]]$summary.random$id_m
    mod_m$model_n <- k
    
    # Summary 
    # mod hyper
    mod_sum <- summary(univ_model_list[[k]])
    
    # Save mean and sd for each variable
    if(length(strsplit(cov,"_")[[1]])>2){
      mean_val <- mean(unscaled_data[[cov_scaled]], na.rm = TRUE)
      sd_val <- sd(unscaled_data[[cov_scaled]], na.rm = TRUE)
    }else{
      mean_val <- mean(unscaled_data[[cov]], na.rm = TRUE)
      sd_val <- sd(unscaled_data[[cov]], na.rm = TRUE)
    }
    
    # Nonlinear run that
    if (length(univ_model_list[[k]]$summary.random) > 3) {
      mod_var <- univ_model_list[[k]]$summary.random[[4]]
      mod_var$model_n <- k
      
      model_out[[k]] <- list(
        data = list(cov = cov, mean = mean_val, sd = sd_val),
        gof = gof,
        summary.fixed = fixed,
        summary.random = list(id_u = mod_spat, id_y = mod_y, id_m = mod_m, mod_var = mod_var, mod_sum = mod_sum)
      )
    } else {
      model_out[[k]] <- list(
        data = list(cov = cov, mean = mean_val, sd = sd_val),
        gof = gof,
        summary.fixed = fixed,
        summary.random = list(id_u = mod_spat, id_y = mod_y, id_m = mod_m, mod_sum = mod_sum)
      )
    }
  }
  
  return(model_out)
}


############## EXTRACT DATA FROM CROSSPREDITCT FOR GPSC INTERACTION IN DLNM####
# Create a function to extract predictions into a tidy format
extract_cp_gpsc_data <- function(cp, level_label, gpsc_name, cov_oi) {
  
  mat_cp_lag <- NULL
  mat_cp <- matrix(nrow=c(nrow(cp$matfit)), ncol = 13)
  for(k in 1:ncol(cp$matfit)){
    mat_cp[,1] <- cp$matfit[,k]
    mat_cp[,2] <- cp$matlow[,k]
    mat_cp[,3] <- cp$mathigh[,k]
    mat_cp[,4] <- rownames(cp$matfit)
    mat_cp[,5] <- gpsc_name
    mat_cp[,6] <- colnames(cp$matfit)[k]
    mat_cp[,7] <- k-1
    mat_cp[,8] <- level_label
    mat_cp[,9] <- cov_oi
    mat_cp[,10] <- cp$allfit
    mat_cp[,11] <- cp$alllow
    mat_cp[,12] <- cp$allhigh
    mat_cp[,13] <- cp$predvar
    mat_cp_lag <- rbind(mat_cp_lag, mat_cp)
  }
  mat_cp_lag <- as.data.frame(mat_cp_lag)
  colnames(mat_cp_lag) <- c("fit","lowerCI","upperCI","var", "GPSC" ,"lag", "lag_num", "interaction_level", "covariate", "cumulative_fit", "cum_lowerCI", "cum_upperCI","predvar")
  mat_cp_lag[,c(1:4,7,10:13)] <- sapply(mat_cp_lag[,c(1:4,7,10:13)],  as.numeric)
  return(mat_cp_lag)
}
