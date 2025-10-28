
####LOAD DATA & LIBRARIES #####################################################
path_to_package <- "/home/sbelman/Documents/Extra_Projects/IDExtremes/GHRmodel/ghrmodel_0.0.0.9000.tar.gz"
################################################################################
#### PURPOSE ##########
################################################################################
## THIS SCRIPT REPLACES 05_1_ WITHOUT INTERACTIONS BUT RUNNING SEPARATE MODELS 
## FOR SEPARATE OUTCOMES WITH PM2.5, PM10, AND HURS. 
## IT INCLUDES THE BASE SCRIPT WHEREBY A DLNM IS RUN FOR EACH VARIABLE INDEPENDENTLY 
## NOT RUNNING THIS WITH INTERACTIONS BUT JUST SIMPLE OUTCOMES

install.packages(path_to_package , 
                 repos = NULL, type = "source", INSTALL_opts = c("--no-multiarch", "--no-test-load"))

library(ghrmodel)
source("/home/sbelman/Documents/env_sa_manuscript/scripts2/0_source_functions.R")

### set resolution
time = "weekly"
space = "adm2"
outcome_type = "province" ## options are "serotype" "demographic" "province"
precov = TRUE

## load spatial data
if(space == "adm1"){
  shp<-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_1.shp")
  ## read adjacency matrix
  g <- inla.read.graph(filename = "/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/sa_adjacency_map_adm1.adj")
}
if(space == "adm2"){
  shp<-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_2.shp")
  ## read adjacency matrix
  g <- inla.read.graph(filename = "/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/sa_adjacency_map.adj")
}

# load  data depending on aggregations
# data <-fread(file=paste0("/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly_lag_sc.csv"))
# data_unscaled <- fread(file=paste0("/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly_lag.csv"))

if(time == "weekly" & space == "adm1"){
  data <-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly_lag_sc.csv")
  data_unscaled <- fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly_lag.csv")
}
if(time == "weekly" & space == "adm2"){
  data <-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm2_weekly_lag_sc.csv")
  data_unscaled <- fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm2_weekly_lag.csv")
}
if(time == "monthly" & space == "adm1"){
  data <-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_monthly_lag_sc.csv")
  data_unscaled <- fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_monthly_lag.csv")
}
if(time == "monthly" & space == "adm2"){
  data <-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm2_monthly_lag_sc.csv")
  data_unscaled <- fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm2_monthly_lag.csv")
}

### set variables and rename appropriately
if("date"%notin%colnames(data)){
  print("Rename Week to Date")
  df<- data
  df <- df %>%
    mutate(date = as.Date(week))}else{
      print("Date already there")
      df<-data
      df$date<-as.Date(data$date)
    }
df <- df %>%
  mutate(
    post_vaccination_2009 = ifelse(date >= as.Date("2009-01-01"), 1, 0),
    post_vaccination_2011 = ifelse(date >= as.Date("2011-01-01"), 1, 0),
    post_covid_2020 = ifelse(date >= as.Date("2020-01-01"), 1, 0),
    vaccination_period = ifelse(date < as.Date("2009-01-01"), 1, ifelse(date >= as.Date("2011-01-01"), 3, 2)),
    time_since_vaccination_2009 = as.numeric(date - as.Date("2009-01-01")) * post_vaccination_2009,
    time_since_vaccination_2011 = as.numeric(date - as.Date("2011-01-01")) * post_vaccination_2011,
    time_since_covid_2020 = as.numeric(date - as.Date("2020-01-01")) * post_covid_2020,
    female = female_count,
    male = male_count
  ) %>%
  mutate(
    across(
      starts_with("GPSC") & ends_with("_count"),
      ~ ifelse(. > 0, 1, 0),
      .names = "{.col}_present"
    )
  )

colnames(df)[grep("present",colnames(df))] <- gsub("_count","", colnames(df)[grep("present",colnames(df))])

### include province as factors for replications
df$id_prov <- as.numeric(factor(df$NAME_1, levels = c("Eastern_Cape", "Free_State", "Gauteng", 
                                                      "KwaZulu-Natal", "Limpopo", "Mpumalanga", 
                                                      "North_West", "Northern_Cape", "Western_Cape")))

df$vaccination_period <- as.factor(df$vaccination_period)

df$week <- as.Date(df$week)  # strips data.table class and metadata
df <- df %>%
  group_by(week, NAME_1) %>%
  summarise(week_seq_prov = sum(sequenced), .groups = 'drop') %>%
  mutate(
    week = as.Date(week),
    NAME_1 = as.character(NAME_1))%>%
  left_join(df, by = c("week","NAME_1")) 

## subset by only pre covid
if(precov==TRUE){
  df <- subset(df,df$date < as.Date("2020-01-01"))
  endyear = 2019
}else{
  endyear = 2023
}


hurs_grp <- inla.group(df$hurs_lag0, n = 5)  # You can adjust 'n' (number of groups) depending on granularity
df$hurs_grp <- hurs_grp

### add environmental covariates
all <- colnames(df)
all_gpscs <- all
all <- grep("lag0",all, value = TRUE)
cov_names <- grep("hurs|pm2p5|pm10|so2", all, value = TRUE)
cov_names_labels <- gsub("_lag0", "", cov_names)

mod_all <- dlnm_results <- NULL
outcomes_demo <- c("mening_count","bact_count","other_count","female_count","age_lt6","age_15t64","age_18t64", "age_gt65","pcv13","pcv7","nvt")
## seros with >2000
data2<- fread(file="/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_point_base.csv",quote=FALSE, header = TRUE)
dtsero <- data.table(table(data2$serotype))[data.table(table(data2$serotype))$N>2000]$V1
dtsero <- dtsero[dtsero!=""]

if(outcome_type=="demographic"){
  outcomes <- outcomes_demo
}
if(outcome_type == "serotype"){
if(space=="adm1"){
  outcomes_sero <- paste0("Sero",dtsero,"_count")
  # outcomes_sero <- grep("Sero",all_gpscs, value=TRUE)
  outcomes <- c(outcomes_sero)
  outcomes <- outcomes[outcomes!="Sero6A_count"] ## will not converge
}
if(space=="adm2"){
  outcomes_sero <- paste0("Sero",dtsero,"_count")
  outcomes <- c(outcomes_sero)
  outcomes <- outcomes[outcomes!="Sero6B_count"] ## will not converge
}
}
############################### RESULTS WITH DIFFERENT OUTCOMES  ###############
for(o in 1:length(outcomes)){
  print(outcomes[o])
  print(paste0("outcome: ", o,"/",length(outcomes)))
  df$outcome <- NULL
  df$outcome <- df[[outcomes[o]]]
  
  
  ################################################################################
  ###################### RUN BASE MODELS #########################################
  ################################################################################
  ## run an intercept model
  forms <- NULL
  forms$formula <- as.formula(outcome ~1)
  forms$cov <- 'int'
  int_mod <- inla.mod(form = forms, fam = "nbinomial", df , nthreads = 4, config = FALSE)
  ## run a random effects only model
  forms <- NULL
  forms$formula <- as.formula(paste0("outcome ~ ",
                                     'f(id_u, model = "bym2", graph = g, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) + ',  # Spatial effect per outcome
                                     'f(id_m, model = "rw2", cyclic = TRUE, scale.model = TRUE, constr = TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) + ',# seasonal effect per outcome
                                     'f(id_y, model = "iid",replicate = id_prov,  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +', #interannual effect per outcome
                                     'vaccination_period +', 'population_density'
                                     ))
  forms$cov <- 're'
  re_mod <- inla.mod(form = forms, fam = "nbinomial", df, nthreads = 4, config = FALSE)
  
  saveRDS(int_mod, file = paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/suboutcomes/intercept_model_",outcomes[o],"_",time,"_",space,"_",endyear,".rds"))
  saveRDS(re_mod, file = paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/suboutcomes/re_only_model_",outcomes[o],"_",time,"_",space,"_",endyear,".rds"))
  ###################### END RUN BASE MODELS #########################################
  
  
############################### START LOOP THROUGH COVARIATES  #################
model_out <- cp_list  <- plot_heatmap <- plot_slices <- list()
for(c in 1:length(cov_names)) {
  print(paste0("Running Cov: ", c, "/", length(cov_names)))
  cov_oi <- cov_names[c]
  print(cov_oi)
  
  ############################### PREPARE CROSSBASIS  ##########################
  # Creating a crossbasis with a vector and allow the function do the lags including the group for me
  if(time=="weekly"){
    max_lag <- 8
    lag_knots <- c(2,4) # Log-spaced knots
    lag_knots <- logknots(max_lag, 2)
  }else{
    max_lag <- 3
    lag_knots <- c(1.5,2.5) # Log-spaced knots
    
  }
  
  if(cov_oi%in%c("prlrsum_lag0","prlrmean_lag0")){
    exp_knots <- logknots(max(df[[cov_oi]]),2)
  }else{
    exp_knots <- quantile(df[[cov_oi]], probs = c(0.33, 0.66), na.rm = TRUE)
  }
  
  if(cov_oi%in%c("hurs_lag0")){
    cb <- crossbasis(
      x = df[[cov_oi]],
      lag = max_lag,
      argvar = list(fun = "ns", knots = exp_knots),
      arglag = list(fun = "ns", df = 3),
      # arglag = list(fun = "ns", knots = lag_knots),
      group = df$id_u
    )
  }else{
    cb <- crossbasis(
      x = df[[cov_oi]],
      lag = max_lag,
      argvar = list(fun = "ns", knots = exp_knots),
      arglag = list(fun = "ns", df = 2),
      # arglag = list(fun = "ns", knots = lag_knots),
      group = df$id_u
    )
  }
  
  
  # create data frame with crossbasis
  df_pre_bound <- df
  cb_df <- cbind(df_pre_bound, cb)
  
  # write formula
  cb_form<- list()
  cb_form$formula <-
    as.formula(
      paste(
        "outcome ~",
        ### only the crossbasis
        paste(colnames(cb), collapse = " + "),  # Include all crossbasis variables
        "+ f(id_u, model = 'bym2', graph = g, scale.model = T, adjust.for.con.comp = TRUE,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
        "+ f(id_m, model = 'rw2', cyclic = T, scale.model = T, constr = T,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
        "+ f(id_y, model = 'iid', replicate = id_prov,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
        "+ vaccination_period",
        "+ population_density"
      )
    )
  cb_form$covs<-paste0("crossbasis_", cov_oi, "_none")
  
  
  ######## RUN MODELS WITH CROSSBASIS ##########################################
  print("running INLA")
  mod <- INLA::inla(
    formula = cb_form$formula,
    family = "nbinomial",  # One per outcome
    offset=log(population_size/100000),
    control.inla = list(strategy = 'adaptive'),
    control.compute = list(dic = TRUE, waic=TRUE, cpo=TRUE, config = FALSE, return.marginals = TRUE),
    control.predictor = list(link = 1, compute = TRUE),
    control.fixed = list(correlation.matrix = TRUE),
    inla.setOption(num.threads = 4),
    verbose = FALSE,
    data = cb_df)
  mod$cov <- cov_names[c]
  
  ######## CROSSPREDICTION AND PLOT ############################################
  print("extract covariance")
  ### set center based on South Africa 
  if(cov_oi%in%c("pm2p5_lag0","pm10_lag0","so2_lag0")){
    if(cov_oi=="pm2p5_lag0"){
      cen = 40
    }
    if(cov_oi=="pm10_lag0"){
      cen = 75
    }
    # if(cov_oi=="o3_lag0"){
    #   cen = 100
    # }
    if(cov_oi=="so2_lag0"){
      cen = 125
    }
  }else{
    ### set center at the median for the rest
    cen = ((range(data[[cov_oi]],na.rm = T)[2]-range(data[[cov_oi]],na.rm = T)[1])/2)+range(data[[cov_oi]],na.rm = T)[1]
  }
  ### extract covariance and variance
  original_coefs <- mod$summary.fixed$mean[c(1:ncol(cb)+1)]  # Extract fixed-effect coefficients
  original_vcov <- mod$misc$lincomb.derived.covariance[1:ncol(cb)+1,1:ncol(cb)+1] # Extract variance-covariance
  
  print("run crosspred no interaction")
  cp <- crosspred(basis = cb, cen=cen, coef = original_coefs, vcov = original_vcov) ## have to center it somewhere because its not linear so its not just at time 0.
  
  ## extract mean and sd
  sd_1 <- sd(data_unscaled[[cov_names[c]]], na.rm = T)
  mean_1 <- mean(data_unscaled[[cov_names[c]]], na.rm=T)
  ## rescale if necessary
  if (grepl("^(hurs|absh|prlrsum|prlrmean|prlrmax)(_lag[0-12])?$", cov_names[c])) {
    print(paste("Rescaling:",cov_names[c]))
    new_predvar<-cp$predvar
    new_predvar<-(cp$predvar*sd_1)+mean_1
    cp$predvar <- new_predvar
    rownames(cp$matfit) <- new_predvar
    
    df_all <- extract_cp_gpsc_data(cp,"none",outcomes[o], cov_names[c])
    dlnm_results <- rbind(dlnm_results,df_all)
    
  }else{
    print(paste("No need to rescale", cov_names[c]))
    df_all <- extract_cp_gpsc_data(cp,"none",outcomes[o], cov_names[c])
    dlnm_results <- rbind(dlnm_results,df_all)
    write.table(df_all, file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/suboutcomes/dlnm_results_",cov_oi,"_",outcomes[o],"_",time,"_",space,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
  }

  ########################### EXTRACT MODEL RESULTS TO SAVE ############
  cov<-mod$cov
  #gof
  gof <- eval.mod(mod,df)
  gof$rsq <- rsq(gof, int_mod, 1)
  # fixed effects
  fixed<-mod$summary.fixed
  #random effects
  mod_spat<-mod$summary.random$id_u
  mod_y<-mod$summary.random$id_y
  mod_m<-mod$summary.random$id_m
  # mod hyper
  mod_sum <- summary(mod)
  # save mean and sd for each variable
  mean <- mean(data_unscaled[[cov]],na.rm=T)
  sd <- sd(data_unscaled[[cov]],na.rm=T)
  
  cp$cov <- cov_names[c]
  cp_list[[c]] <- cp
  model_out[[c]]<-list(data=list(cov=cov,mean=mean,sd=sd),gof=gof, fixed=fixed,summary.random=list(id_u=mod_spat,id_y=mod_y,id_m=mod_m,  mod_sum=mod_sum), cp = cp) 
} ### end of the loop through covariates
############## SUMMARY OF COVARIATES ###################################
mod_sum <- lapply(model_out, function(x) x$gof)
mod_sum <- do.call(rbind, mod_sum)
mod_sum2 <- mod_sum[,c("cov","waic","mae","cpo","rsq")]
mod_sum2$cov <- gsub("_lag0","",mod_sum2$cov)
mod_sum2$outcome <- outcomes[o]
mod_all <- rbind(mod_all, mod_sum2)
############## SAVE FILES ##############################################
saveRDS(model_out,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/suboutcomes/model_out_summary_list_",outcomes[o],"_",time,"_",space,"_",endyear,".rds"))
saveRDS(cp_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/suboutcomes/crosspred_list_",outcomes[o],"_",time,"_",space,"_dlnm_",endyear,".rds"))
write.table(mod_sum2, file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/suboutcomes/mod_gof_dlnm_",outcomes[o],"_",time,"_",space,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")


} ## end loop through outcomes

write.table(dlnm_results, file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/suboutcomes/dlnm_subsetoutcomes_results_fits_",time,"_",space,"_",outcome_type,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(mod_all, file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/suboutcomes/mod_gof_dlnm_alloutcomes_",time,"_",space,"_",outcome_type,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")






# ################################################################################
# ####### NO INTERACTION EFFECT PLOTS (spatiotemporal resolutions) ############
# ################################################################################
# cov_names <- unique(dlnm_results$cov)
# cov_names <- c("tasmax","hurs","absh","pm2p5","pm10")
# ax_labs <- c("Max. Temperature","Relative Humidity (%)","Absolute Humidity", 
#              expression(paste('Concentration (', mu, 'g/m'^3, ')')),
#              expression(paste('Concentration (', mu, 'g/m'^3, ')')))
# dlnm_results$cov <- gsub("_lag0","",dlnm_results$covariate)
# ## select vectors
# dis_vec <- c("mening_count","bact_count")
# age_vec <- c("age_lt5","age_gt65","age_gt80","female_count")
# sero_vec <- grep("Sero",unique(dlnm_results$GPSC),value=TRUE)
# pcv_vec <- c("pcv13","pcv7","nvt")
# ### plot all covs and 6 weeks of lags
# dlnm_results_sub <- subset(dlnm_results, dlnm_results$cov=="pm2p5")
# plotcov_list <- list()
# # for(c in 1:length(cov_names)){
#   tmp <- subset(dlnm_results_sub, dlnm_results_sub$cov%in%cov_names[c] & dlnm_results_sub$lag_num %in% c(0:6) & dlnm_results_sub$GPSC%in%sero_vec )
#   tmp$rr <- exp(tmp$fit)
#   tmp$lowerCI_rr <- exp(tmp$lowerCI)
#   tmp$upperCI_rr <- exp(tmp$upperCI)
#   tmp$lowerCI_rr[which(tmp$lowerCI_rr<0.8)] <- 0.8
#   tmp$upperCI_rr[which(tmp$upperCI_rr>1.1)] <- 1.1
#   tmp$lag_week <- paste0("Week",tmp$lag_num)
#   tmp$lag_week <- factor(tmp$lag_week, levels = paste0("Week",seq(0,6,1)), labels = paste0("Week",seq(0,6,1)))
#   plotcov <- ggplot(tmp)+
#     geom_line(aes(x=var,y=rr,group=lag_week,group=GPSC))+
#     geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
#     geom_ribbon(aes(x=var, ymin=lowerCI_rr,ymax=upperCI_rr, group=lag_week,group=GPSC),alpha=0.3)+
#     theme_bw()+
#     xlab("Var")+
#     ylab("Relative Risk")+
#     # scale_y_continuous(trans="log10")+
#     scale_y_continuous(trans="log10", limits = c(0.8,1.11), breaks = c(0.8, 1, 1.1))+
#     ggtitle(cov_names[c])+
#     facet_grid(GPSC~lag_week, ncol=7)+
#     xlab(ax_labs[c])+
#     theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13))
#   plotcov_list[[c]] <- plotcov
# # }
# 
# 
# plotcov_list[[1]]
# plotcov_list[[2]]
# plotcov_list[[3]]
# plotcov_list[[4]]
# plotcov_list[[5]]
# 
# 
# #################### DLNM ######################################################
# p1 <- ggplot(dlnm_results_sub)+
#   geom_hline(yintercept=1, linetype="dashed",color="red")+
#   geom_line(aes(x = predvar, y = exp(cumulative_fit), group=cov))+
#   geom_ribbon(aes(x = predvar, ymin = exp(cum_lowerCI), ymax= exp(cum_upperCI), group=cov), alpha=0.5)+
#   theme_bw()+
#   ylab("Relative Risk")+
#   ggtitle(outcomes[1])+
#   ylim(0.1,2)+
#   xlab("Variable")+
#   facet_wrap(cov~., scales="free")
# dlnm_results_sub <- subset(dlnm_results, dlnm_results$GPSC == outcomes[2])
# p2 <- ggplot(dlnm_results_sub)+
#   geom_hline(yintercept=1, linetype="dashed",color="red")+
#   geom_line(aes(x = predvar, y = exp(cumulative_fit), group=cov))+
#   geom_ribbon(aes(x = predvar, ymin = exp(cum_lowerCI), ymax= exp(cum_upperCI), group=cov), alpha=0.5)+
#   theme_bw()+
#   ylab("Relative Risk")+
#   ggtitle(outcomes[2])+
#   ylim(0.1,2)+
#   xlab("Variable")+
#   facet_wrap(cov~., scales="free")
# p1/p2
