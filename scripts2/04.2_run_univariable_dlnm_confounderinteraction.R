################################################################################
#### PURPOSE ##########
################################################################################
## THIS SCRIPT IS TESTING THE INCORPORATING OF A CONFOUNDING INTERACTION WITH TEMPERATURE AND PM2.5
## WE ALSO MANUALLY SUBSTITUTE INCLUDING HURS AND ABSH (5 CUTS SPECIFICED AS ABSH_GRP AND HURS_GRP)
## ADDITIONALLY, TESTING INCORPORATING TEMPERATURE AS AN ADDITIONAL CROSS BASIS. ## 

## allows models from 2005-2019 and 2005-2023
## also can run with and without the interaction of the environmental effects and at all spatial and temporal levels

setwd("/home/sbelman/Documents/env_sa_manuscript/")
####LOAD DATA & LIBRARIES #####################################################
source("scripts2/0_source_functions.R")
### set if interaction is true or not
interaction = TRUE
### set resolution
time = "weekly"
space = "adm2"
precov = TRUE
permute = FALSE
## load spatial data
if(space == "adm1"){
  shp<-st_read("input_datasets/shps/gadm41_namematch_ZAF_1.shp")
  ## read adjacency matrix
  g <- inla.read.graph(filename = "input_datasets/shps/sa_adjacency_map_adm1.adj")
}
if(space == "adm2"){
  shp<-st_read("input_datasets/shps/gadm41_namematch_ZAF_2.shp")
  ## read adjacency matrix
  g <- inla.read.graph(filename = "input_datasets/shps/sa_adjacency_map.adj")
}

# load  data depending on aggregations
# data <-fread(file=paste0("dataframes/sa_adm1_weekly_lag_sc.csv"))
# data_unscaled <- fread(file=paste0("dataframes/sa_adm1_weekly_lag.csv"))

if(time == "weekly" & space == "adm1"){
  data <-fread(file="dataframes/sa_adm1_weekly_lag_sc.csv")
  data_unscaled <- fread(file="dataframes/sa_adm1_weekly_lag.csv")
}
if(time == "weekly" & space == "adm2"){
  data <-fread(file="dataframes/sa_adm2_weekly_lag_sc.csv")
  data_unscaled <- fread(file="dataframes/sa_adm2_weekly_lag.csv")
}
if(time == "monthly" & space == "adm1"){
  data <-fread(file="dataframes/sa_adm1_monthly_lag_sc.csv")
  data_unscaled <- fread(file="dataframes/sa_adm1_monthly_lag.csv")
}
if(time == "monthly" & space == "adm2"){
  data <-fread(file="dataframes/sa_adm2_monthly_lag_sc.csv")
  data_unscaled <- fread(file="dataframes/sa_adm2_monthly_lag.csv")
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
    vaccination_period = ifelse(date < as.Date("2009-01-01"), 1, ifelse(date >= as.Date("2011-01-01"), 3, 2)),
    post_covid_2020 = ifelse(date >= as.Date("2020-01-01"), 1, 0),
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

## include the proportion of each GPSC per the number sequenced each month
# df <- df %>%
#   group_by(year, month, NAME_1) %>%
#   summarise(month_seq_prov = sum(sequenced), .groups = 'drop') %>%
#   left_join(df, by = c("year","month","NAME_1")) 
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
  data_unscaled <- subset(data_unscaled,data_unscaled$date < as.Date("2020-01-01"))
  endyear = 2019
}else{
  endyear = 2023
}

## make a group for humidity
# Create grouped version
hurs_grp <- inla.group(df$hurs_lag3, n = 5)  # You can adjust 'n' (number of groups) depending on granularity
df$hurs_grp <- hurs_grp
absh_grp <- inla.group(df$absh_lag3, n = 5)  # You can adjust 'n' (number of groups) depending on granularity
df$absh_grp <- absh_grp


### add environmental covariates
all <- colnames(df)
all_gpscs <- all
all <- grep("lag0",all, value = TRUE)
if(interaction == TRUE){
  # cov_names <- grep("tasmax|hurs|absh|pm2p5|pm10|o3|so2", all, value = TRUE)
  cov_names <- grep("pm2p5|pm10|so2", all, value = TRUE)
}else{
# cov_names <- grep("hurs|pm2p5|pm10", all, value = TRUE)
cov_names <- grep("pm2p5|pm10|so2", all, value = TRUE)
}
cov_names_labels <- gsub("_lag0", "", cov_names)

# 
# ### select some serotypes to include
# data2<- fread(file="input_datasets/disease/SA_disease_point_base.csv",quote=FALSE, header = TRUE)
# dtsero <- data.table(table(data2$serotype))[data.table(table(data2$serotype))$N>800]
# pcv_vec <- c("4","6B","9V","14","18C","19F","23F","1","3","5","6A","7F","19A")
# dtsero[which(dtsero$V1%notin%pcv_vec)]

### select which environmental factors will be included
# dtgpsc <- data.table(table(data2$GPSC))[which(data.table(table(data2$GPSC))$N>100 | data.table(table(data2$GPSC))$V1%in%c("8","41"))]
# dtgpsc[order(-N)]$V1
# gpsc_vec <- paste0("GPSC",dtgpsc$V1,"_count")
envint_vec <- c("tas")
envint_vec <- paste0(envint_vec,"_lag3")

##### LOAD INTERCEPT MODELS FOR R2 CALCULATION ##################################
int_mod <- readRDS(file=paste0("models/base_models/base_model_main_",time,"_intercept_",space,"_",endyear,".rds"))
re_mod <- readRDS(file=paste0("models/base_models/base_model_",time,"_20092011_popdens_",space,"_",endyear,".rds"))

##### SET UP LOOPS FOR MODELS    ###############################################
## define loop length
    if(interaction==TRUE){
      envint_vec_sub <- envint_vec
      mod_sum_all <- envint_results <- NULL
    }else{
      envint_vec_sub <- "none"
      dlnm_results <- NULL
    }

### initiate envint loop or loop through single "none" vector for no interaction
for(gp in 1:length(envint_vec_sub)){

  print(paste0("Running Models for: ", envint_vec_sub[gp]))
    ### set the interaction variable
    interact_var = envint_vec_sub[gp]
        ### create empty lists
        model_out <- cp_list  <- plot_heatmap <- plot_slices <- list()
        for(c in 1:length(cov_names)) {
          print(paste0("Running Cov: ", c, "/", length(cov_names)))
          cov_oi <- cov_names[c]
          print(cov_oi)
          
        ############################### PREPARE CROSSBASIS  ############################
          # Creating a crossbasis with a vector and allow the function do the lags including the group for me
          if(time=="weekly"){
            max_lag <- 8
            lag_knots <- c(2,4) # Log-spaced knots
            lag_knots <- logknots(max_lag, 2)
          }else{
            max_lag <- 3
            lag_knots <- c(1.5,2.5) # Log-spaced knots
          }
          
          # if(cov_oi%in%c("prlrsum_lag0","prlrmean_lag0")){
          #   exp_knots <- logknots(max(df[[cov_oi]]),2)
          # }else{
            # exp_knots <- quantile(df[[cov_oi]], probs = c(0.33, 0.66), na.rm = TRUE)
          # }
          
          exp_knots <- quantile(df[[cov_oi]], probs = c(0.33, 0.66), na.rm = TRUE)
            
          
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
          
          ####### set up crossbasis with and without interaction
          if (interaction == TRUE) {
            ## save df
            df_pre_bound <- df
            # ## put the environmental weekly over the amount sequenced per month so that the GPSCs are comparable
            # if (grepl("^GPSC", interact_var)) {
            #   prop <- ifelse(df_pre_bound$week_seq_prov == 0, 0,
            #                  df_pre_bound[[interact_var]] / df_pre_bound$week_seq_prov)
            #   
            # } else {
            #   stop(paste("Unknown interaction variable:", interact_var))
            # }
            interactnums <- df_pre_bound[[interact_var]] 
            
            ## for non gpsc interaction
            if(grepl("^(hurs|absh|prlrsum|prlrmean|prlrmax)(_lag[0-8])?$", interact_var)){
              scaled_interactnums <- scale(interactnums, center = TRUE, scale = FALSE)
              interact_var_scaled <- as.numeric(scaled_interactnums)
              sd_1 <- sd(data_unscaled[[interact_var]], na.rm = T)
              mean_1 <- mean(data_unscaled[[interact_var]], na.rm=T)
              # if scaling
              # original_interactnums <- interact_var_scaled * sd_1 + mean_1
              original_interactnums <- data_unscaled[[interact_var]]
            }else{
              # just scale and don't remove zeroes as we did with the GPSCs
              scaled_interactnums <- scale(interactnums, center = TRUE, scale = TRUE)
              interact_var_scaled <- as.numeric(scaled_interactnums)
              sd_1 <- sd(data_unscaled[[interact_var]], na.rm = T)
              mean_1 <- mean(data_unscaled[[interact_var]], na.rm=T)
              # if only centering
              original_interactnums <- interact_var_scaled*sd_1  + mean_1
            }

            ## define z given scaled interaction
            z_low <- quantile(interact_var_scaled, 0.15, na.rm = T)
            z_med <- quantile(interact_var_scaled, 0.5, na.rm = T)
            z_high <- quantile(interact_var_scaled, 0.85, na.rm = T)
            
            z_low_true <- quantile(original_interactnums , 0.15, na.rm = T)
            z_med_true <- quantile(original_interactnums , 0.5, na.rm = T)
            z_high_true <- quantile(original_interactnums , 0.85, na.rm = T)
            
            # Step 3: apply interaction as usual
            cbinteract <- cb * interact_var_scaled
            # cbinteract[is.na(cbinteract)] <- 0
            
            # Name the columns to track the interaction
            colnames(cbinteract) <- paste0(colnames(cb), ":", interact_var)
            
            # Combine everything into the same dataset
            cb_df <- cbind(df_pre_bound, cb, cbinteract)
          }
            if(interaction==FALSE){
            # create data frame with crossbasis
            df_pre_bound <- df
            cb_df <- cbind(df_pre_bound, cb)
            }
        ###################### PREPARE FORMULA #########################################
          if(interaction==TRUE){
          # write formula
          cb_form<- list()
          cb_form$formula <-
            as.formula(
              paste(
                "disease ~",
                paste(c(
                  colnames(cb),                 # main effect of crossbasis - this creates instability in the model due to colinearity
                  interact_var,                 # main effect of envint proportion
                  colnames(cbinteract)),         # interaction terms
                  collapse = " + "),
                ### only the crossbasis
                # paste(colnames(cb), collapse = " + "),  # Include all crossbasis variables
                "+ f(id_u, model = 'bym2', graph = g, scale.model = T, adjust.for.con.comp = TRUE, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
                "+ f(id_m, model = 'rw2', cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
                "+ f(id_y, model = 'iid', replicate = id_prov,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
                "+ vaccination_period",
                "+ population_density",
                paste0("+ f(absh_grp, model = 'rw2', scale.model = T, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))" )
              )
            )
          cb_form$covs<-paste0("crossbasis_", cov_oi, "_interaction_", gsub("_count","",interact_var))
        }else{
          # write formula
          cb_form<- list()
          cb_form$formula <-
            as.formula(
              paste(
                "disease ~",
                ### only the crossbasis
                paste(colnames(cb), collapse = " + "),  # Include all crossbasis variables
                "+ f(id_u, model = 'bym2', graph = g, scale.model = T, adjust.for.con.comp = TRUE,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
                "+ f(id_m, model = 'rw2', cyclic = T, scale.model = T, constr = T,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
                "+ f(id_y, model = 'iid', replicate = id_prov,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
                "+ vaccination_period",
                "+ population_density",
                paste0("+ f(absh_grp, model = 'rw2', scale.model = T, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))" )
              )
            )
          cb_form$covs<-paste0("crossbasis_", cov_oi, "_none")
        }
        
        ######## RUN MODELS WITH CROSSBASIS #############################################
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
          
          if(cov_names[c] %in% c("hurs_lag0","pm2p5_lag0","pm10_lag0")){
            # saveRDS(mod, file = paste0("models/dlnms/modifiers/dlnm_model_univariable_",gsub("_lag0","",cov_names[c]),"_",space,"_",time,".rds"))
          }
        
        ######## CROSSPREDICTION AND PLOT ##############################################
          print("extract covariance")
          ### set center
          if(cov_oi%in%c("pm2p5_lag0","pm10_lag0","o3_lag0","so2_lag0")){
            if(cov_oi=="pm2p5_lag0"){ cen = 40 }
            if(cov_oi=="pm10_lag0"){ cen = 75 }
            # if(cov_oi=="o3_lag0"){ cen = 100 }
            if(cov_oi=="so2_lag0"){ cen = 125 }
          }else{
            ### set center at the median for the rest
            cen = ((range(data[[cov_oi]],na.rm = T)[2]-range(data[[cov_oi]],na.rm = T)[1])/2)+range(data[[cov_oi]],na.rm = T)[1]
          }          ### extract covariance and variance
          original_coefs <- mod$summary.fixed$mean[c(1:ncol(cb)+1)]  # Extract fixed-effect coefficients
          original_vcov <- mod$misc$lincomb.derived.covariance[1:ncol(cb)+1,1:ncol(cb)+1] # Extract variance-covariance
       
          ### extract interacting covariance and variance
          if(interaction==TRUE){
            original_coefs_interact <- mod$summary.fixed$mean[grep(paste0(":",interact_var),rownames(mod$summary.fixed))]
            original_vcov_interact <- mod$misc$lincomb.derived.covariance[grep(paste0(":",interact_var),rownames(mod$summary.fixed)),
                                                                          grep(paste0(":",interact_var),rownames(mod$summary.fixed))] 
            
            #### for cross variance matrix between both
            coef_names<- colnames(mod$misc$lincomb.derived.covariance)
            main_effects <-coef_names[c(1:ncol(cb)+1)]
            interaction_terms <- grep(":", coef_names, value = TRUE)
            # Extract the covariance submatrix between main effects and interaction terms
            main_effects_idx <- which(coef_names %in% main_effects)  # Get indices for main effects
            interaction_idx <- which(coef_names %in% interaction_terms)  # Get indices for interaction terms
            # Now extract the cross-covariance terms (off-diagonal elements)
            original_vcov_cross <-  mod$misc$lincomb.derived.covariance[main_effects_idx, interaction_idx]
            
            ## crosspred with interacting cross basis
            print("run crosspred for interaction")
            
            ### MOVED THIS TO ABOVE
            ## define z given scaled interaction
            # z_low <- quantile(scaled_interactnums, 0.25)
            # z_med <- quantile(scaled_interactnums, 0.5)
            # z_high <- quantile(scaled_interactnums, 0.75)
            # 
            # z_low_true <- quantile(original_interactnums , 0.2)
            # z_med_true <- quantile(original_interactnums , 0.5)
            # z_high_true <- quantile(original_interactnums , 0.75)
            
            ## interaction crossprediciton
            coef_low <- original_coefs + z_low * original_coefs_interact
            coef_med <- original_coefs + z_med * original_coefs_interact
            coef_high <- original_coefs + z_high * original_coefs_interact
            
            # Calculate the variance-covariance matrices for these coefficients
            vcov_low <- original_vcov + 
              (z_low^2) * original_vcov_interact + 
              2 * z_low * original_vcov_cross
            
            vcov_med <- original_vcov + 
              (z_med^2) * original_vcov_interact + 
              2 * z_med * original_vcov_cross
            
            vcov_high <- original_vcov + 
              (z_high^2) * original_vcov_interact + 
              2 * z_high * original_vcov_cross
  
            ## crosspred for each low and high proportions
            # cp_low <- crosspred(basis = cbinteract, cen = cen,
            #                     coef = coef_low,
            #                     vcov = vcov_low)
            # cp_med <- crosspred(basis = cbinteract, cen = cen,
            #                      coef = coef_med,
            #                      vcov = vcov_med)
            # cp_high <- crosspred(basis = cbinteract, cen = cen,
            #                      coef = coef_high,
            #                      vcov = vcov_high)
            
            cp_low <- crosspred(basis = cb, cen = cen,
                                coef = coef_low,
                                vcov = vcov_low)
            cp_med <- crosspred(basis = cb, cen = cen,
                                coef = coef_med,
                                vcov = vcov_med)
            cp_high <- crosspred(basis = cb, cen = cen,
                                 coef = coef_high,
                                 vcov = vcov_high)

            }else{
            print("run crosspred no interaction")
            cp <- crosspred(basis = cb, cen=cen, coef = original_coefs, vcov = original_vcov) ## have to center it somewhere because its not linear so its not just at time 0.
          }
          ## extract mean and sd
          sd_1 <- sd(data_unscaled[[cov_names[c]]], na.rm = T)
          mean_1 <- mean(data_unscaled[[cov_names[c]]], na.rm=T)
          
          ##################### RESCALE CROSS PREDICTION FROM SCALED VARIABLES #########
          if(interaction==TRUE){
            if (grepl("^(hurs|absh|prlrsum|prlrmean|prlrmax)(_lag[0-12])?$", cov_names[c])) {
              print(paste("Rescaling:",cov_names[c]))
              # new_predvar<-cp$predvar
              # new_predvar<-(cp$predvar*sd_1)+mean_1
              # cp$predvar <- new_predvar
              # rownames(cp$matfit) <- new_predvar
              
              ## low proportion predictions
              new_predvar<-cp_low$predvar
              new_predvar<-(cp_low$predvar*sd_1)+mean_1
              cp_low$predvar <- new_predvar
              rownames(cp_low$matfit) <- new_predvar
              
              ## high proportion predictions
              new_predvar<-cp_med$predvar
              new_predvar<-(cp_med$predvar*sd_1)+mean_1
              cp_med$predvar <- new_predvar
              rownames(cp_med$matfit) <- new_predvar
              
              ## high proportion predictions
              new_predvar<-cp_high$predvar
              new_predvar<-(cp_high$predvar*sd_1)+mean_1
              cp_high$predvar <- new_predvar
              rownames(cp_high$matfit) <- new_predvar
              
          
              # Extract risks
              envint_name <- gsub("_lag0","",envint_vec_sub[gp])
              df_low <- extract_cp_gpsc_data(cp_low, "low", envint_name, cov_names[c])
              df_low$envmodifier_val <- z_low_true
              df_med <- extract_cp_gpsc_data(cp_med, "med", envint_name, cov_names[c])
              df_med$envmodifier_val <- z_med_true
              df_high <- extract_cp_gpsc_data(cp_high, "high", envint_name, cov_names[c])
              df_high$envmodifier_val <- z_high_true
              
              envint_results <- rbind(envint_results, df_low,df_med, df_high)
              

            }else{
              print(paste("No need to rescale", cov_names[c]))
              # Extract risks
              envint_name <- gsub("_count","",envint_vec[gp])
              df_low <- extract_cp_env_data(cp_low, "low", envint_name , cov_names[c])
              df_low$envmodifier_val <- z_low_true
              
              df_med <- extract_cp_env_data(cp_med, "med", envint_name, cov_names[c])
              df_med$envmodifier_val <- z_med_true
              
              df_high <- extract_cp_env_data(cp_high, "high", envint_name, cov_names[c])
              df_high$envmodifier_val <- z_high_true
              
              envint_results <- rbind(envint_results, df_low, df_med, df_high)
              
            }
          }else{
            if (grepl("^(hurs|absh|prlrsum|prlrmean|prlrmax)(_lag[0-8])?$", cov_names[c])) {
              print(paste("Rescaling:",cov_names[c]))
              new_predvar<-cp$predvar
              new_predvar<-(cp$predvar*sd_1)+mean_1
              cp$predvar <- new_predvar
              rownames(cp$matfit) <- new_predvar
              
              df_all <- extract_cp_env_data(cp,"none","no_interaction", cov_names[c])
              dlnm_results <- rbind(dlnm_results,df_all)
              
            }else{
              print(paste("No need to rescale", cov_names[c]))
              df_all <- extract_cp_gpsc_data(cp,"none","no_interaction", cov_names[c])
              dlnm_results <- rbind(dlnm_results,df_all)
            }
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
          ### save all
          if(interaction==TRUE){
            cp_high$cov <- cp_low$cov <- cp_med$cov <- cov_names[c]
            # cp_list[[c]] <- list(cp = cp, cp_low = cp_low, cp_med = cp_med ,cp_high = cp_high)
            cp_list[[c]] <- list(cp_low = cp_low, cp_med = cp_med ,cp_high = cp_high)
            model_out[[c]]<-list(data=list(cov=cov,mean=mean,sd=sd),gof=gof, fixed=fixed,summary.random=list(id_u=mod_spat,id_y=mod_y,id_m=mod_m,  mod_sum=mod_sum), cp = list(cp_low = cp_low, cp_med = cp_med, cp_high = cp_high))
          }else{
            cp$cov <- cov_names[c]
            cp_list[[c]] <- cp
            model_out[[c]]<-list(data=list(cov=cov,mean=mean,sd=sd),gof=gof, fixed=fixed,summary.random=list(id_u=mod_spat,id_y=mod_y,id_m=mod_m,  mod_sum=mod_sum), cp = cp) 
          }
        } ### end of the loop through covariates
        
        ############## SUMMARY OF COVARIATES ###################################
        mod_sum <- lapply(model_out, function(x) x$gof)
        mod_sum <- do.call(rbind, mod_sum)
        mod_sum2 <- mod_sum[,c("cov","waic","mae","cpo","rsq")]
        mod_sum2$cov <- gsub("_lag0","",mod_sum2$cov)
        mod_sum2$interact <- gsub("_count","",interact_var)
        
        if(interaction==TRUE){
          mod_sum_all <- rbind(mod_sum_all,mod_sum2)
        }
        
        if(permute == TRUE){
          
        }else{
        ############## SAVE FILES ##############################################
        # saveRDS(model_out,file=paste0("models/dlnms/modifiers/model_out_summary_list_",time,"_",space,"_dlnm_",interact_var,"_",max_lag,"_",endyear,".rds"))
        # saveRDS(cp_list,file=paste0("models/dlnms/modifiers/crosspred_list_",time,"_",space,"_dlnm_",interact_var,"_",maxlag,"_",endyear,".rds"))
        write.table(mod_sum2, file=paste0("models/dlnms/modifiers/mod_gof_dlnm_",time,"_",space,"_",interact_var,"_",max_lag,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
}
  } ### end of loop through modifier


if(interaction==TRUE){
write.table(envint_results, file=paste0("models/dlnms/modifiers/envmod_results_fits_",time,"_",space,"_allmodifiers_propprov_",max_lag,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = ",")
write.table(mod_sum_all, file=paste0("models/dlnms/modifiers/mod_gof_dlnm",time,"_",space,"_allmodifiers_propprov_",max_lag,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = ",")
}else{
  write.table(dlnm_results, file=paste0("models/dlnms/modifiers/nointeraction_results_fits_",time,"_",space,"_",max_lag,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = ",")
}

################################################################################
time = "weekly"
space = "adm1"
max_lag = 8
endyear = 2019
## PLOTS ##
## pm2.5 effect with temperature variables (min and max)
dlnm_results <- fread(file=paste0("models/dlnms/modifiers/nointeraction_results_fits_",time,"_",space,"_",max_lag,"_",endyear,".csv"))
colnames(dlnm_results)[which(colnames(dlnm_results)=="GPSC")] <- "envmodifier"
dlnm_results$envmodifier_val <- NA
envint_results <- fread(file=paste0("models/dlnms/modifiers/envmod_results_fits_",time,"_",space,"_allmodifiers_propprov_",max_lag,"_",endyear,".csv"))
### bind them as if they're the same
outdf <- rbind(dlnm_results,envint_results)
### join them so I can plot on the same graph
nointdt <- dlnm_results[,c("fit","lowerCI","upperCI","lag","lag_num","cumulative_fit","cum_lowerCI","cum_upperCI","predvar")]
colnames(nointdt) <- c("fit_noint","lowerCI_noint","upperCI_noint","lag","lag_num","cumulative_fit_noint","cum_lowerCI_noint","cum_upperCI_noint","predvar")
envint_none <- left_join(envint_results,nointdt, by = c("lag","lag_num","predvar"))
# outdf <- envint_results[which(envint_results$envmodifier%in%c("tasmin_lag0","tasmax_lag0","tas_lag0","so2_lag0","hurs_lag0","absh_lag0")),]
envint_none$interaction_level <- factor(envint_none$interaction_level, levels = c("low","med","high"), labels =  c("low","med","high"))

## modifier names
mod_vec <- unique(envint_none$envmodifier)
mod_vecname <- gsub("_lag0","",mod_vec)
envint_none$envmodifier <- factor(envint_none$envmodifier, levels = mod_vec, labels = mod_vecname)
mod_lagsplot <-ggplot(envint_none[which(envint_none$interaction_level%in%c("low","high"))]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  geom_line(aes(x=var, y=fit, group = interaction(interaction_level, envmodifier, lag_num), color = interaction_level)) +
  geom_ribbon(aes(x=var, ymin = lowerCI, ymax = upperCI, group = interaction(interaction_level, envmodifier, lag_num), fill = interaction_level), alpha = 0.2) +
  geom_line(aes(x=var, y=fit_noint, group = interaction(interaction_level, envmodifier, lag_num)), color = "black", linetype = "dashed") +
  theme_bw() +
  xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
  ylab("Effect") + 
  labs(fill = "effect modifier\nproportion", color = "effect modifier\nproportion")+
  theme(axis.text.x = element_text(angle=45, hjust =1))+
  facet_grid(envmodifier ~ lag_num, scales ="free")



outdfcum <- envint_none[which(envint_none$lag_num ==3),]
cummod_plot <- ggplot(outdfcum) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  geom_line(aes(x=var, y=cumulative_fit, group = interaction(interaction_level, envmodifier), color = interaction_level)) +
  geom_ribbon(aes(x=var, ymin = cum_lowerCI, ymax = cum_upperCI, group = interaction(interaction_level, envmodifier), fill = interaction_level), alpha = 0.2) +
  geom_line(aes(x=var, y=cumulative_fit_noint, group = interaction(interaction_level, envmodifier, lag_num)), color = "black", linetype = "dashed") +
  theme_bw() +
  xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
  ylab("Effect") + 
  labs(fill = "effect modifier\nproportion", color = "effect modifier\nproportion")+
  facet_wrap(envmodifier ~ .)


envmodifier_table <- envint_none %>%
  dplyr::select(envmodifier, interaction_level, envmodifier_val) %>%
  distinct() %>%
  arrange(envmodifier, interaction_level)
envmodifier_table2 <- envmodifier_table %>%   # replace df with your object name
  pivot_wider(
    names_from = interaction_level,
    values_from = envmodifier_val
  ) %>%
  mutate(across(c(low, med, high), ~ signif(.x, 2))) %>%
  arrange(envmodifier)


write.table(envmodifier_table2, file = "figures/revision_figs/modifier_proportions.csv", sep = ",",quote = FALSE, row.names = FALSE, col.names = TRUE)
##### save sensitivity plots
pdf("figures/revision_figs/cum_meteomodplot.pdf", width = 8, height =5)
print(cummod_plot)
ggsave("figures/revision_figs/cum_meteomodplot.png", width = 8, height =5)
dev.off()
pdf("figures/revision_figs/weekly_meteomodplot.pdf", width = 10, height =6)
print(mod_lagsplot)
ggsave("figures/revision_figs/weekly_meteomodplot.png", width = 10, height =6)
dev.off()


outdfcum_rr <- outdfcum |>
  dplyr::mutate(
    rr_fit   = exp(cumulative_fit),
    rr_low   = exp(cum_lowerCI),
    rr_high  = exp(cum_upperCI),
    rr_fit_noint   = exp(cumulative_fit_noint),
    rr_low_noint  = exp(cum_lowerCI_noint),
    rr_high_noint  = exp(cum_upperCI_noint)
  )

ggplot(outdfcum_rr) +
  geom_line(aes(x=var, y=rr_fit, group = interaction(interaction_level, envmodifier), color = interaction_level)) +
  geom_ribbon(aes(x=var, ymin = rr_low, ymax = rr_high, group = interaction(interaction_level, envmodifier), fill = interaction_level), alpha = 0.2) +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")+
  facet_wrap(envmodifier ~ .) +
  ylab("relative risk") +
  xlab("pm2.5 concentration")



###### run simple linear interaction inla model ######

### humidity ###
form<- list()
form$formula <-
  as.formula(
    paste(
      "disease ~",
      "pm2p5_lag3*absh_lag0",  # Include all crossbasis variables
      "+ f(id_u, model = 'bym2', graph = g, scale.model = T, adjust.for.con.comp = TRUE,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
      "+ f(id_m, model = 'rw2', cyclic = T, scale.model = T, constr = T,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
      "+ f(id_y, model = 'iid', replicate = id_prov,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
      "+ vaccination_period",
      "+ population_density"
    )
  )
form$covs<-"pm2.5_hurs"
mod <- INLA::inla(
  formula = form$formula,
  family = "nbinomial",  # One per outcome
  offset=log(population_size/100000),
  control.inla = list(strategy = 'adaptive'),
  control.compute = list(dic = TRUE, waic=TRUE, cpo=TRUE, config = FALSE, return.marginals = TRUE),
  control.predictor = list(link = 1, compute = TRUE),
  # control.fixed = list(correlation.matrix = TRUE),
  inla.setOption(num.threads = 4),
  verbose = FALSE,
  data = df)
mod$cov <- "int"

mod$summary.fixed
beta_pm  <- mod$summary.fixed["pm2p5_lag3","mean"]
beta_int <- mod$summary.fixed["pm2p5_lag3:absh_lag0","mean"]

low_h  <- quantile(df$absh_lag0, 0.1)
high_h <- quantile(df$absh_lag0, 0.9)

effect_low  <- beta_pm + beta_int * low_h
effect_high <- beta_pm + beta_int * high_h

c(low = effect_low, high = effect_high)
### %stronger effect at high humidity as compared to low humidity
effect_high/effect_low



### temperature ##
form<- list()
form$formula <-
  as.formula(
    paste(
      "disease ~",
      "pm2p5_lag3*tasmax_lag0",  # Include all crossbasis variables
      "+ f(id_u, model = 'bym2', graph = g, scale.model = T, adjust.for.con.comp = TRUE,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
      "+ f(id_m, model = 'rw2', cyclic = T, scale.model = T, constr = T,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
      "+ f(id_y, model = 'iid', replicate = id_prov,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
      "+ vaccination_period",
      "+ population_density"
    )
  )
form$covs<-"pm2.5_tasmin"


print("running INLA")
mod <- INLA::inla(
  formula = form$formula,
  family = "nbinomial",  # One per outcome
  offset=log(population_size/100000),
  control.inla = list(strategy = 'adaptive'),
  control.compute = list(dic = TRUE, waic=TRUE, cpo=TRUE, config = FALSE, return.marginals = TRUE),
  control.predictor = list(link = 1, compute = TRUE),
  # control.fixed = list(correlation.matrix = TRUE),
  inla.setOption(num.threads = 4),
  verbose = FALSE,
  data = df)
mod$cov <- "int"

mod$summary.fixed
beta_pm  <- mod$summary.fixed["pm2p5_lag3","mean"]
beta_int <- mod$summary.fixed["pm2p5_lag3:tasmax_lag0","mean"]

low_h  <- quantile(df$tasmax_lag0, 0.1)
high_h <- quantile(df$tasmax_lag0, 0.9)

effect_low  <- beta_pm + beta_int * low_h
effect_high <- beta_pm + beta_int * high_h

c(low = effect_low, high = effect_high)

### percent stronger effect of PM2.5 at higher temperatures as compared to lower ones
effect_high/effect_low



