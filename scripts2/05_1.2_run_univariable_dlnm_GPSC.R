################################################################################
#### PURPOSE ##########
################################################################################
## THIS SCRIPT REPLACES 05_1_ FOR THE GPSC INTERACTIONS AT A WEEKLY ADMIN 1 OR ADMIN 2LEVEL
## IT INCLUDES THE BASE SCRIPT WHEREBY A DLNM IS RUN FOR EACH VARIABLE INDEPENDENTLY 
## AND THEN ITERATIVELY INCLUDING EACH GPSC TO DETERMINE IF THE GPSC HAS A MODIFYING 
## EFFECT ON THE ENVIRONMENTAL VARIABLE. INCLUDING THE GPSC PROPORTIONS VIA A BINOMIAL
## PROBABILITY MODEL. 
## We are running this at weekly administrative region 1 WITH THE INTERACTION to
## maximize the power of the temporal effect but not reducing the power too much.


####LOAD DATA & LIBRARIES #####################################################
path_to_package <- "/home/sbelman/Documents/Extra_Projects/IDExtremes/GHRmodel/ghrmodel_0.0.0.9000.tar.gz"

install.packages(path_to_package , 
                 repos = NULL, type = "source", INSTALL_opts = c("--no-multiarch", "--no-test-load"))
library(ghrmodel)
source("/home/sbelman/Documents/BRD/scripts/0_source_functions.R")
### set if interaction is true or not
interaction = FALSE
### set resolution
time = "weekly"
space = "adm2"

## load spatial data
if(space == "adm1"){
  shp<-st_read("/home/sbelman/Documents/BRD/SouthAfrica/shps/gadm41_namematch_ZAF_1.shp")
  ## read adjacency matrix
  g <- inla.read.graph(filename = "/home/sbelman/Documents/BRD/SouthAfrica/shps/sa_adjacency_map_adm1.adj")
}
if(space == "adm2"){
  shp<-st_read("/home/sbelman/Documents/BRD/SouthAfrica/shps/gadm41_namematch_ZAF_2.shp")
  ## read adjacency matrix
  g <- inla.read.graph(filename = "/home/sbelman/Documents/BRD/SouthAfrica/shps/sa_adjacency_map.adj")
}

# load  data depending on aggregations
# data <-fread(file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm1_weekly_lag_sc.csv"))
# data_unscaled <- fread(file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm1_weekly_lag.csv"))

if(time == "weekly" & space == "adm1"){
  data <-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm1_weekly_lag_sc.csv")
  data_unscaled <- fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm1_weekly_lag.csv")
}
if(time == "weekly" & space == "adm2"){
  data <-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_weekly_lag_sc.csv")
  data_unscaled <- fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_weekly_lag.csv")
}
if(time == "monthly" & space == "adm1"){
  data <-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm1_monthly_lag_sc.csv")
  data_unscaled <- fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm1_monthly_lag.csv")
}
if(time == "monthly" & space == "adm2"){
  data <-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_monthly_lag_sc.csv")
  data_unscaled <- fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_monthly_lag.csv")
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

### subset by time we have sequenced
df <- subset(df,df$date <= as.Date("2021-01-01"))

### add environmental covariates
all <- colnames(df)
all_gpscs <- all
all <- grep("lag0",all, value = TRUE)
if(interaction == TRUE){
  # cov_names <- grep("tasmax|hurs|absh|pm2p5|pm10|o3|so2", all, value = TRUE)
  cov_names <- grep("hurs|pm2p5|pm10", all, value = TRUE)
}else{
# cov_names <- grep("tasmax|tasmin|hurs|absh|prlrsum|prlrmean|sfcWind|pm2p5|pm10|o3|so2|spei3|spei6|spi3|spi6", all, value = TRUE)
cov_names <- grep("hurs|pm2p5|pm10", all, value = TRUE)
}
cov_names_labels <- gsub("_lag0", "", cov_names)

### select some serotypes to include
data2<- fread(file="/home/sbelman/Documents/BRD/SouthAfrica/disease/SA_disease_point_base.csv",quote=FALSE, header = TRUE)
dtsero <- data.table(table(data2$serotype))[data.table(table(data2$serotype))$N>800]
pcv_vec <- c("4","6B","9V","14","18C","19F","23F","1","3","5","6A","7F","19A")
dtsero[which(dtsero$V1%notin%pcv_vec)]

### select which GPSCs will be includes
# gpsc_vec <- grep("GPSC", all_gpscs, value = TRUE)
# gpsc_vec <- grep("count", gpsc_vec, value = TRUE) ## if including the proportions
dtgpsc <- data.table(table(data2$GPSC))[which(data.table(table(data2$GPSC))$N>100 | data.table(table(data2$GPSC))$V1%in%c("8","41"))]
dtgpsc[order(-N)]$V1
gpsc_vec <- paste0("GPSC",dtgpsc$V1,"_count")

##### VISUALIZE DATA ###########################################################
## visualize the weekly admin 1 gpsc data
# gpsc_annual_plot <- gpsc_season_plot <- list()
# for(gp in 1:length(gpsc_vec)){
# gpsc_season_plot[[gp]] <- df %>%
#   group_by(month)%>%
#   summarise(count = sum(!!sym(gpsc_vec[gp])),
#             month_tot = sum(week_seq_prov), .groups = 'drop') %>%
#   ungroup() %>%
#   ggplot()+
#     geom_line(aes(x=month,y=count/month_tot))+
#     ggtitle(gpsc_vec[gp])+
#     theme(legend.position = "none")
# 
# gpsc_annual_plot[[gp]] <- df %>%
#   group_by(date)%>%
#   summarise(count = sum(!!sym(gpsc_vec[gp])),
#             month_tot = sum(week_seq_prov), .groups = 'drop') %>%
#   ungroup() %>%
#   ggplot()+
#   geom_line(aes(x=date,y=count/month_tot))+
#   ggtitle(gpsc_vec[gp])+
#   theme(legend.position = "none")
# }
# grid.arrange(grobs=gpsc_season_plot)
# grid.arrange(grobs=gpsc_annual_plot)

##### LOAD INTERCEPT MODELS FOR R2 CALCULATION ##################################
int_mod <- readRDS(file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/base_models/base_model_main_",time,"_intercept_",space,".rds"))
re_mod <- readRDS(file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/base_models/base_model_",time,"_20092011_popdens_",space,".rds"))

##### SET UP LOOPS FOR MODELS    ###############################################
## define loop length
    if(interaction==TRUE){
      gpsc_vec_sub <- gpsc_vec
      rr_ratio_all <- mod_sum_all <- gpsc_results <- NULL
    }else{
      gpsc_vec_sub <- "none"
      dlnm_results <- NULL
    }

### initiate GPSC loop or loop through single "none" vector for no interaction
for(gp in 1:length(gpsc_vec_sub)){
  print(paste0("Running Models for: ",gpsc_vec_sub[gp]))
    ### set the interaction variable
    interact_var = gpsc_vec_sub[gp]
        ### create empty lists
        model_out <- cp_list  <- plot_heatmap <- plot_slices <- list()
        for(c in 1:length(cov_names)) {
          print(paste0("Running Cov: ", c, "/", length(cov_names)))
          cov_oi <- cov_names[c]
          print(cov_oi)
          
        ############################### PREPARE CROSSBASIS  ############################
          # Creating a crossbasis with a vector and allow the function do the lags including the group for me
          if(time=="weekly"){
            max_lag <- 6
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
          
          ####### set up crossbasis with and without interaction
          if (interaction == TRUE) {
            ## save df
            df_pre_bound <- df
            ## put the GPSC weekly over the amount sequenced per month so that the GPSCs are comparable
            if (grepl("^GPSC", interact_var)) {
              prop <- ifelse(df_pre_bound$week_seq_prov == 0, 0,
                             df_pre_bound[[interact_var]] / df_pre_bound$week_seq_prov)
              
            } else {
              stop(paste("Unknown interaction variable:", interact_var))
            }
            
            # Create the interacted crossbasis
            # Step 1: scale only non-zero props
            prop_nonzero <- prop[prop > 0]
            scaled_nonzero <- scale(prop_nonzero, center = TRUE, scale = TRUE)
            
            # Step 2: put scaled values back into a vector matching full dataset length
            interact_var_scaled <- rep(0, length(prop))  # default to 0
            interact_var_scaled[prop > 0] <- scaled_nonzero
            
            # Step 3: apply interaction as usual
            cbinteract <- cb * interact_var_scaled
            cbinteract[is.na(cbinteract)] <- 0
            
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
                  interact_var,                 # main effect of GPSC proportion
                  colnames(cbinteract)),         # interaction terms
                  collapse = " + "),
                ### only the crossbasis
                # paste(colnames(cb), collapse = " + "),  # Include all crossbasis variables
                "+ f(id_u, model = 'bym2', graph = g, scale.model = T, adjust.for.con.comp = TRUE, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
                "+ f(id_m, model = 'rw2', cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
                "+ f(id_y, model = 'iid', replicate = id_prov,  hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
                "+ post_vaccination_2009",
                "+ post_vaccination_2011",
                "+ population_density"
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
                "+ post_vaccination_2009",
                "+ post_vaccination_2011",
                "+ population_density"
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
            # saveRDS(mod, file = paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/models/dlnm_model_univariable_",gsub("_lag0","",cov_names[c]),"_",space,"_",time,".rds"))
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
            
            ## define z given scaled interaction
            z_low <- quantile(scaled_nonzero, 0.1)
            z_med <- quantile(scaled_nonzero, 0.5)
            z_high <- quantile(scaled_nonzero, 0.8)
            
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
            cp_low <- crosspred(basis = cbinteract, cen = cen,
                                coef = coef_low,
                                vcov = vcov_low)
            cp_med <- crosspred(basis = cbinteract, cen = cen,
                                 coef = coef_med,
                                 vcov = vcov_med)
            cp_high <- crosspred(basis = cbinteract, cen = cen,
                                 coef = coef_high,
                                 vcov = vcov_high)
            
            
            # ## Relative Risk incorporating the variance covariacne
            # X_high <- cp_high$coefficients
            # X_low <- cp_low$coefficients
            # vcov_model <- original_vcov_interact
            # X_diff <- X_high - X_low
            # SE_diff <- sqrt(t(X_diff) %*% vcov_model %*% X_diff)
            # fit_high <- cp_high$allfit
            # fit_low  <- cp_low$allfit
            # RR_ratio <- exp(fit_high - fit_low)
            # diff <- fit_high - fit_low
            # RR_ratio_lower <- exp(diff - 1.96 * SE_diff)
            # RR_ratio_upper <- exp(diff + 1.96 * SE_diff)
            # rr_ratiodf <- data.frame(cbind(names(RR_ratio),RR_ratio, RR_ratio_lower, RR_ratio_upper))
            # rr_ratiodf$gpsc <- gpsc_vec_sub[gp]
            # rr_ratiodf$cov <- cov_names[c]
            # colnames(rr_ratiodf) <- c("var", "rr_ratio","rr_lowerCI","rr_upperCI", "gpsc","cov")
            # rr_ratiodf[,2:4] <- sapply(rr_ratiodf[,2:4], as.numeric)
            # 
            }else{
            print("run crosspred no interaction")
            cp <- crosspred(basis = cb, cen=cen, coef = original_coefs, vcov = original_vcov) ## have to center it somewhere because its not linear so its not just at time 0.
          }
          ## extract mean and sd
          sd_1 <- sd(data_unscaled[[cov_names[c]]], na.rm = T)
          mean_1 <- mean(data_unscaled[[cov_names[c]]], na.rm=T)
          
          ##################### RESCALE CROSS PREDICTION FROM SCALED VARIABLES #########
          if(interaction==TRUE){
            if (grepl("^(hurs|absh|prlrsum|prlrmean|prlrmax)(_lag[0-6])?$", cov_names[c])) {
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
              gpsc_name <- gsub("_count","",gpsc_vec_sub[gp])
              df_low <- extract_cp_gpsc_data(cp_low, "low", gpsc_name, cov_names[c])
              df_med <- extract_cp_gpsc_data(cp_med, "med", gpsc_name, cov_names[c])
              df_high <- extract_cp_gpsc_data(cp_high, "high", gpsc_name, cov_names[c])
              gpsc_results <- rbind(gpsc_results, df_low,df_med, df_high)
              
              # relative risk ratio results
              ## rename RR
              # rr_ratiodf$var <- new_predvar
              # rr_ratio_all <- rbind(rr_ratio_all, rr_ratiodf)
              
            }else{
              print(paste("No need to rescale", cov_names[c]))
              # Extract risks
              gpsc_name <- gsub("_count","",gpsc_vec[gp])
              df_low <- extract_cp_gpsc_data(cp_low, "low", gpsc_name, cov_names[c])
              df_med <- extract_cp_gpsc_data(cp_med, "med", gpsc_name, cov_names[c])
              df_high <- extract_cp_gpsc_data(cp_high, "high", gpsc_name, cov_names[c])
        
              gpsc_results <- rbind(gpsc_results, df_low, df_med, df_high)
              
              # relative risk ratio results
              # rr_ratio_all <- rbind(rr_ratio_all, rr_ratiodf)
              
            }
          }else{
            if (grepl("^(hurs|absh|prlrsum|prlrmean|prlrmax)(_lag[0-6])?$", cov_names[c])) {
              print(paste("Rescaling:",cov_names[c]))
              new_predvar<-cp$predvar
              new_predvar<-(cp$predvar*sd_1)+mean_1
              cp$predvar <- new_predvar
              rownames(cp$matfit) <- new_predvar
              
              df_all <- extract_cp_gpsc_data(cp,"none","no_interaction", cov_names[c])
              dlnm_results <- rbind(dlnm_results,df_all)
              
            }else{
              print(paste("No need to rescale", cov_names[c]))
              df_all <- extract_cp_gpsc_data(cp,"none","no_interaction", cov_names[c])
              dlnm_results <- rbind(dlnm_results,df_all)
            }
          }
          
          # plot_crosspred_coef(cp, type = "slices", lag = 0:8)
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
        
        ############## SAVE FILES ##############################################
        saveRDS(model_out,file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/model_out_summary_list_",time,"_",space,"_dlnm_",interact_var,"_",max_lag,"week_mixeddf.rds"))
        # saveRDS(cp_list,file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/crosspred_list_",time,"_",space,"_dlnm_",interact_var,"_",maxlag,"week.rds"))
        write.table(mod_sum2, file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/mod_gof_dlnm_",time,"_",space,"_",interact_var,"_",max_lag,"week_mixeddf.csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")

  } ### end of loop through gpscs


if(interaction==TRUE){
write.table(rr_ratio_all, file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/rr_ratio_all_",time,"_",space,"_allGPSCs_propprov_",max_lag,"week_mixeddf.csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(gpsc_results, file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/gpsc_results_fits_",time,"_",space,"_allGPSCs_propprov_",max_lag,"week_mixeddf.csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(mod_sum_all, file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/mod_gof_dlnm",time,"_",space,"_allGPSCs_propprov_",max_lag,"week_mixeddf.csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
}else{
  write.table(dlnm_results, file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/nointeraction_results_fits_",time,"_",space,"_",max_lag,"week_mixeddf.csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
}


# ################################################################################





