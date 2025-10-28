################################################################################
#### PURPOSE ##########
################################################################################
## IT INCLUDES THE BASE SCRIPT WHEREBY A DLNM IS RUN FOR EACH VARIABLE INDEPENDENTLY 
## AND THEN ITERATIVELY INCLUDING EACH GPSC TO DETERMINE IF THE GPSC HAS A MODIFYING 
## EFFECT ON THE ENVIRONMENTAL VARIABLE. INCLUDING THE GPSC PROPORTIONS VIA A BINOMIAL
## PROBABILITY MODEL. 
## We are running this at weekly administrative region 2 to
## maximize the power of the temporal effect but not reducing the power too much.


####LOAD DATA & LIBRARIES #####################################################
path_to_package <- "/home/sbelman/Documents/Extra_Projects/IDExtremes/GHRmodel/ghrmodel_0.0.0.9000.tar.gz"

install.packages(path_to_package , 
                 repos = NULL, type = "source", INSTALL_opts = c("--no-multiarch", "--no-test-load"))
library(ghrmodel)
source("/home/sbelman/Documents/env_sa_manuscript/scripts2/0_source_functions.R")
### set if interaction is true or not
interaction = FALSE
### set resolution
time = "weekly"
space = "adm1"
### set the time period 2005 - 2019 is precov and 2005-2023 is not precov
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
  endyear = 2019
}else{
  endyear = 2023
}
### add environmental covariates
all <- colnames(df)
all_gpscs <- all
all <- grep("lag0",all, value = TRUE)
if(interaction == TRUE){
  # cov_names <- grep("tasmax|hurs|absh|pm2p5|pm10|o3|so2", all, value = TRUE)
  cov_names <- grep("hurs|absh|pm2p5|pm10|prlrsum", all, value = TRUE)
  cov_names2 <- grep("hurs|absh|tasmax|pm2p5|pm10|so2|prlrsum", all, value = TRUE)
  
}else{
  if(space == "adm2"){
  cov_names <- grep("hurs|absh|pm2p5|pm10|prlrsum", all, value = TRUE)
  cov_names2 <- grep("hurs|absh|tasmax|tasmin|pm2p5|pm10|so2|prlrsum", all, value = TRUE)
  }else{
    cov_names <- grep("hurs|absh|pm2p5|pm10", all, value = TRUE)
    cov_names2 <- grep("hurs|absh|tasmax|tasmin|pm2p5|pm10|so2", all, value = TRUE)
  }
  }
cov_names_labels <- gsub("_lag0", "", cov_names)
### select which GPSCs will be includes
gpsc_vec <- grep("GPSC", all_gpscs, value = TRUE)
gpsc_vec <- grep("count", gpsc_vec, value = TRUE) ## if including the proportions 

### select some serotypes to include
data2<- fread(file="/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_point_base.csv",quote=FALSE, header = TRUE)
dtsero <- data.table(table(data2$serotype))[data.table(table(data2$serotype))$N>800]
pcv_vec <- c("4","6B","9V","14","18C","19F","23F","1","3","5","6A","7F","19A")
dtsero[which(dtsero$V1%notin%pcv_vec)]

## make a group for humidity
# Create grouped version
hurs_grp <- inla.group(df$hurs_lag0, n = 5)  # You can adjust 'n' (number of groups) depending on granularity
df$hurs_grp <- hurs_grp
absh_grp <- inla.group(df$absh_lag0, n = 5)  # You can adjust 'n' (number of groups) depending on granularity
df$absh_grp <- absh_grp


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
int_mod <- readRDS(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_main_",time,"_intercept_",space,"_",endyear,".rds"))
re_mod <- readRDS(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_",time,"_20092011_popdens_",space,"_",endyear,".rds"))


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
          model_out[[c]] <- list()
          print(paste0("Running Cov: ", c, "/", length(cov_names)))
          for(c2 in 1:length(cov_names2)){
            if(cov_names2[c2]==cov_names[c]){next}
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
          
          if(cov_oi%in%c("prlrsum_lag0","prlrmean_lag0")){
            exp_knots <- logknots(max(df[[cov_oi]]),2)
          }else{
            exp_knots <- quantile(df[[cov_oi]], probs = c(0.33, 0.66), na.rm = TRUE)
          }
          
          if(cov_oi%in%c("hurs_lag0","absh_lag0")){
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
            if(cov_names2[c2]%in%c("hurs_lag0")){
              # write formula for hurs
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
                    "+ vaccination_period",
                    "+ population_density",
                    paste0("+ f(hurs_grp, model = 'rw2', scale.model = T, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))" )
                  )
                )
            }
            if(cov_names2[c2]%in%c("absh_lag0")){
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
                    "+ vaccination_period",
                    "+ population_density",
                    paste0("+ f(absh_grp, model = 'rw2', scale.model = T, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))" )
                  )
                )
            }## end loop not hurs
        if(cov_names2[c2]%notin%c("hurs_lag0","absh_lag0")){
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
                "+ vaccination_period",
                "+ population_density",
                paste("+",cov_names2[c2])
              )
            )
            }## end loop not hurs
            cb_form$covs<-paste0("crossbasis_", cov_oi, "_interaction_", gsub("_count","",interact_var),"_",cov_names2[c2])
            
        }else{
          if(cov_names2[c2]=="hurs_lag0"){
            # write formula for hurs          
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
                paste0("+ f(hurs_grp, model = 'rw2', scale.model = T, hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))" )
              )
            )
          }
          if(cov_names2[c2]=="absh_lag0"){
            # write formula for hurs          
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
          }
          if(cov_names2[c2]%notin%c("hurs_lag0","absh_lag0")){
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
                  paste("+",cov_names2[c2])
                )
              )
          }
          cb_form$covs<-paste0("crossbasis_", cov_oi, "_none_",cov_names2[c2])
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
          mod$cov2 <- cov_names2[c2]
        
        ######## CROSSPREDICTION AND PLOT ##############################################
          print("extract covariance")
          ### set center
          if(cov_oi%in%c("pm2p5_lag0","pm10_lag0","o3_lag0","so2_lag0")){
            if(cov_oi=="pm2p5_lag0"){ cen = 40 }
            if(cov_oi=="pm10_lag0"){ cen = 75 }
            # if(cov_oi=="o3_lag0"){ cen = 100 } ## removing this and centering at the median as it fluctuates day and night and follows an 8 hour cycle
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
            
            
            ### Relative Risk incorporating the variance covariacne
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
              gpsc_name <- gsub("_count","",gpsc_vec_sub[gp])
              df_low <- extract_cp_gpsc_data(cp_low, "low", gpsc_name, cov_names[c])
              df_med <- extract_cp_gpsc_data(cp_med, "med", gpsc_name, cov_names[c])
              df_high <- extract_cp_gpsc_data(cp_high, "high", gpsc_name, cov_names[c])
              gpsc_results <- rbind(gpsc_results, df_low,df_med, df_high)
              
              # relative risk ratio results
              ## rename RR
              rr_ratiodf$var <- new_predvar
              rr_ratio_all <- rbind(rr_ratio_all, rr_ratiodf)
              
            }else{
              print(paste("No need to rescale", cov_names[c]))
              # Extract risks
              gpsc_name <- gsub("_count","",gpsc_vec[gp])
              df_low <- extract_cp_gpsc_data(cp_low, "low", gpsc_name, cov_names[c])
              df_med <- extract_cp_gpsc_data(cp_med, "med", gpsc_name, cov_names[c])
              df_high <- extract_cp_gpsc_data(cp_high, "high", gpsc_name, cov_names[c])
        
              gpsc_results <- rbind(gpsc_results, df_low, df_med, df_high)
              
              # relative risk ratio results
              rr_ratio_all <- rbind(rr_ratio_all, rr_ratiodf)
              
            }
          }else{
            if (grepl("^(hurs|absh|prlrsum|prlrmean|prlrmax)(_lag[0-12])?$", cov_names[c])) {
              print(paste("Rescaling:",cov_names[c]))
              new_predvar<-cp$predvar
              new_predvar<-(cp$predvar*sd_1)+mean_1
              cp$predvar <- new_predvar
              rownames(cp$matfit) <- new_predvar
              
              df_all <- extract_cp_gpsc_data(cp,"none","no_interaction", cov_names[c])
              df_all$cov2 <- cov_names2[c2]
              dlnm_results <- rbind(dlnm_results,df_all)
              
            }else{
              print(paste("No need to rescale", cov_names[c]))
              df_all <- extract_cp_gpsc_data(cp,"none","no_interaction", cov_names[c])
              df_all$cov2 <- cov_names2[c2]
              dlnm_results <- rbind(dlnm_results,df_all)
            }
          }
 
          
          
          ########################### EXTRACT MODEL RESULTS TO SAVE ############
          cov<-mod$cov
          cov2 <- mod$cov2
          #gof
          gof <- eval.mod(mod,df)
          gof$rsq <- rsq(gof, int_mod, 1)
          # fixed effects
          fixed<-mod$summary.fixed
          #random effects
          mod_spat<-mod$summary.random$id_u
          mod_y<-mod$summary.random$id_y
          mod_m<-mod$summary.random$id_m
          if(cov_names2[c2]=="hurs_lag0"){
            mod_hurs <- mod$summary.random$hurs_grp
          }
          # mod hyper
          mod_sum <- summary(mod)
          # save mean and sd for each variable
          mean <- mean(data_unscaled[[cov]],na.rm=T)
          sd <- sd(data_unscaled[[cov]],na.rm=T)
          ### save all
          if(interaction==TRUE){
            cp_high$cov <- cp_low$cov <- cp_med$cov <- cov_names[c]
            cp_list[[c]] <- list(cp = cp, cp_low = cp_low, cp_med = cp_med ,cp_high = cp_high)
            cp_list[[c]] <- list(cp_low = cp_low, cp_med = cp_med ,cp_high = cp_high, rr_ratio = rr_ratiodf)
            model_out[[c]][[c2]]<-list(data=list(cov=cov,cov2 = cov2, mean=mean,sd=sd),gof=gof, fixed=fixed,summary.random=list(id_u=mod_spat,id_y=mod_y,id_m=mod_m,  mod_sum=mod_sum), cp = list(cp_low = cp_low, cp_med = cp_med, cp_high = cp_high, rr_ratio = rr_ratiodf))
          }else{
            cp$cov <- cov_names[c]
            cp_list[[c]] <- cp
            if(cov_names2[c2]=="hurs_lag0"){
              model_out[[c]][[c2]]<-list(data=list(cov=cov, cov2 = cov2, mean=mean,sd=sd),gof=gof, fixed=fixed,summary.random=list(id_u=mod_spat,id_y=mod_y,id_m=mod_m, id_hurs=mod_hurs, mod_sum=mod_sum), cp = cp) 
            }else{
            model_out[[c]][[c2]]<-list(data=list(cov=cov, cov2 = cov2, mean=mean,sd=sd),gof=gof, fixed=fixed,summary.random=list(id_u=mod_spat,id_y=mod_y,id_m=mod_m, mod_sum=mod_sum), cp = cp) 
            }
            }
        } ### end of the loop through covariates
        } ## end of second covariate bivariable interaction loop
        
        ############## SUMMARY OF COVARIATES ###################################
        mod_sum <- do.call(rbind, lapply(model_out, function(lvl1) {
          do.call(rbind, lapply(lvl1, function(x) {
            # Extract gof and add cov2 as a new column
            gof_df <- x$gof
            gof_df$cov2 <- x$data$cov2
            return(gof_df)
          }))
        }))
        mod_sum2 <- mod_sum[,c("cov","waic","mae","cpo","rsq","cov2")]
        mod_sum2$cov <- gsub("_lag0","",mod_sum2$cov)
        mod_sum2$interact <- gsub("_count","",interact_var)
        mod_sum2$cov2 <- gsub("_lag0","",mod_sum2$cov2)
        
        if(interaction==TRUE){
          mod_sum_all <- rbind(mod_sum_all,mod_sum2)
        }
        
        ############## SAVE FILES ##############################################
        saveRDS(model_out,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/bivariable/model_out_summary_list_",time,"_",space,"_dlnm_",interact_var,"_bivariable_",max_lag,"_",endyear,".rds"))
        # saveRDS(cp_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/bivariable/crosspred_list_",time,"_",space,"_dlnm_",interact_var,"_bivariable_",max_lag,"_",endyear,".rds"))
        write.table(mod_sum2, file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/bivariable/mod_gof_dlnm_",time,"_",space,"_",interact_var,"_bivariable_",max_lag,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")

  } ### end of loop through gpscs


if(interaction==TRUE){
write.table(rr_ratio_all, file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/bivariable/rr_ratio_all_",time,"_",space,"_allGPSCs_propprov_bivariable_",max_lag,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(gpsc_results, file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/bivariable/gpsc_results_fits_",time,"_",space,"_allGPSCs_propprov_bivariable_",max_lag,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(mod_sum_all, file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/bivariable/mod_gof_dlnm",time,"_",space,"_allGPSCs_propprov_bivariable_",max_lag,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
}else{
  write.table(dlnm_results, file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/bivariable/nointeraction_results_fits_",time,"_",space,"_bivariable_",max_lag,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
}
 
# 
# ################################################################################
mod_sum2[order(mod_sum2$waic),]
# ##  VISUALIZE BIVARIABLE OUTPUTS 
st_vec <- c("weekly_adm1","weekly_adm2")

### read in bivariable
stlist <- list()
for(s in 1:length(st_vec)){
  m1 <- fread(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/bivariable/nointeraction_results_fits_",st_vec[s],"_bivariable_",max_lag,"_",endyear,".csv"))
  m1$space_time <- st_vec[s]
  stlist[[s]] <- m1
}
fit_list_biv <- rbindlist(stlist)
fit_list_biv$cov <- gsub("_lag0","",fit_list_biv$covariate)
fit_list_biv$cov2 <- gsub("_lag0","",fit_list_biv$cov2)

### read in univariable
stlist <- list()
for(s in 1:length(st_vec)){
  m1 <- fread(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/bivariable/nointeraction_results_fits_",st_vec[s],"_",max_lag,"_",endyear,".csv"))
  m1$space_time <- st_vec[s]
  stlist[[s]] <- m1
}
fit_list <- rbindlist(stlist)
fit_list$cov <- gsub("_lag0","",fit_list$covariate)
fit_list$cov2 <- "none"
fit_list <- rbind(fit_list, fit_list_biv)
fit_list$lag_week <- paste0("Week",fit_list$lag_num)
fit_list$lag_week <- factor(fit_list$lag_week, levels = paste0("Week",seq(0,12,1)), labels = paste0("Week",seq(0,12,1)))
# convert to relative risks
fit_list$rr <- exp(fit_list$fit)
fit_list$lowerCI_rr <- exp(fit_list$lowerCI)
fit_list$upperCI_rr <- exp(fit_list$upperCI)
fit_list$cumrr <- exp(fit_list$cumulative_fit)
fit_list$cumlowerCI_rr <- exp(fit_list$cum_lowerCI)
fit_list$cumupperCI_rr <- exp(fit_list$cum_upperCI)

### Relative humidity with pm10 and pm2.5 and without any ----------------------
tmp <- subset(fit_list, fit_list$cov%in%c("hurs") & fit_list$cov2%in%c("pm2p5","pm10","tasmax","tasmin","none") & fit_list$lag_num %in% c(0:12) & fit_list$space_time=="weekly_adm2")
tmp$cov2 <- factor(tmp$cov2, levels = c("pm2p5","pm10","so2","tasmax","tasmin","none"))
color_labels <- c("pm2p5" = "red", "pm10" = "orange","so2"= "pink", "hurs" = "darkblue","tasmax" = "purple","tasmin" = "brown", "none" = "darkgreen")
linetype_labels <-c("pm2p5" = "solid", "pm10" = "solid","so2" = "solid", "hurs" = "solid","tasmax" = "solid", "tasmin" = "solid", "none" = "longdash")
plot_hurs <- ggplot(tmp)+
  geom_line(aes(x=var,y=rr,group=interaction(lag_week,space_time, cov2), color=cov2, linetype = cov2),alpha=0.7)+
  geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
  geom_ribbon(aes(x=var, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,space_time,cov2), fill=cov2),alpha=0.04)+
  scale_color_manual(values = color_labels)+
  scale_fill_manual(values = color_labels)+
  scale_linetype_manual(values = linetype_labels)+
  theme_bw()+
  xlab("Var")+
  ylab("Relative Risk")+
  # scale_y_continuous(trans="log10")+
  scale_y_continuous(trans="log10", limits = c(0.85,1.15), breaks = c(0.85, 1, 1.15))+
  ggtitle("Relative Humidity")+
  labs(linetype = "additive variable", color = "additive variable", fill = "additive variable")+
  facet_wrap(.~lag_week, nrow = 1)+
  xlab("Percent")+
  theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13))
plot_hurs

tmp <- subset(fit_list, fit_list$cov%in%c("hurs") & fit_list$cov2%in%c("pm2p5","pm10","so2","tasmax","tasmin","none")& fit_list$space_time=="weekly_adm2" & fit_list$lag_num==2)
tmp$cov2 <- factor(tmp$cov2, levels = c("pm2p5","pm10","so2","tasmax","tasmin","none"))
plot_hurs_l2 <- ggplot(tmp)+
  geom_line(aes(x=var,y=cumrr,group=interaction(lag_week,space_time, cov2), color=cov2, linetype = cov2),alpha=0.5)+
  geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
  geom_ribbon(aes(x=var, ymin=cumlowerCI_rr,ymax=cumupperCI_rr, group=interaction(lag_week,space_time,cov2), fill=cov2),alpha=0.04)+
  scale_color_manual(values = color_labels)+
  scale_fill_manual(values = color_labels)+
  scale_linetype_manual(values = linetype_labels)+
  theme_bw()+
  xlab("Var")+
  ylab("Relative Risk")+
  # scale_y_continuous(trans="log10")+
  # scale_y_continuous(trans="log10", limits = c(0.9,1.1), breaks = c(0.9, 1, 1.1))+
  ggtitle("Relative Humidity")+
  labs(linetype = "additive variable", color = "additive variable", fill = "additive variable")+
  # facet_wrap(.~lag_week, nrow = 1)+
  xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
  theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
        strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13),
        legend.position = "none")

### Relative humidity with pm10 and pm2.5 and without any ----------------------
tmp <- subset(fit_list, fit_list$cov%in%c("pm2p5") & fit_list$cov2%in%c("hurs","tasmax","tasmin","none") & fit_list$lag_num %in% c(0:12) & fit_list$space_time=="weekly_adm2")
tmp$cov2 <- factor(tmp$cov2, levels = c("hurs","tasmax","tasmin","none"))
plot_pm2p5 <- ggplot(tmp)+
  geom_line(aes(x=var,y=rr,group=interaction(lag_week,space_time, cov2), color=cov2, linetype = cov2),alpha=0.7)+
  geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
  geom_ribbon(aes(x=var, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,space_time,cov2), fill=cov2),alpha=0.05)+
  theme_bw()+
  xlab("Var")+
  ylab("Relative Risk")+
  # scale_y_continuous(trans="log10")+
  scale_y_continuous(trans="log10", limits = c(0.85,1.15), breaks = c(0.85, 1, 1.15))+
  scale_color_manual(values = color_labels)+
  scale_fill_manual(values = color_labels)+
  scale_linetype_manual(values = linetype_labels)+
  ggtitle("PM 2.5")+
  labs(linetype = "additive variable", color = "additive variable", fill = "additive variable")+
  facet_wrap(.~lag_week, nrow = 1)+
  xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
  theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13))

tmp <- subset(fit_list, fit_list$cov%in%c("pm2p5") & fit_list$cov2%in%c("hurs","tasmax","tasmin","none") & fit_list$lag_num %in% c(0:12) & fit_list$space_time=="weekly_adm2" & fit_list$lag_num==2)
tmp$cov2 <- factor(tmp$cov2, levels = c("hurs","tasmax","tasmin","none"))
plot_pm2p5_l2 <- ggplot(tmp)+
  geom_line(aes(x=var,y=rr,group=interaction(lag_week,space_time, cov2), color=cov2, linetype = cov2),alpha=0.7)+
  geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
  geom_ribbon(aes(x=var, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,space_time,cov2), fill=cov2),alpha=0.05)+
  theme_bw()+
  xlab("Var")+
  ylab("Relative Risk")+
  # scale_y_continuous(trans="log10")+
  scale_y_continuous(trans="log10", limits = c(0.9,1.1), breaks = c(0.9, 1, 1.1))+
  scale_color_manual(values = color_labels)+
  scale_fill_manual(values = color_labels)+
  scale_linetype_manual(values = linetype_labels)+
  ggtitle("PM 2.5")+
  labs(linetype = "additive variable", color = "additive variable", fill = "additive variable")+
  # facet_wrap(.~lag_week, nrow = 1)+
  # xlab(ax_labs[c])+
  theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
        strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13),
        legend.position = "none")

### Relative humidity with pm10 and pm2.5 and without any ----------------------
tmp <- subset(fit_list, fit_list$cov%in%c("pm10") & fit_list$cov2%in%c("hurs","tasmax","tasmin","none") & fit_list$lag_num %in% c(0:12) & fit_list$space_time=="weekly_adm2")
tmp$cov2 <- factor(tmp$cov2, levels = c("hurs","tasmax","tasmin","none"))
plot_pm10 <- ggplot(tmp)+
  geom_line(aes(x=var,y=rr,group=interaction(lag_week,space_time, cov2), color=cov2, linetype=cov2),alpha=0.7)+
  geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
  geom_ribbon(aes(x=var, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,space_time,cov2), fill=cov2),alpha=0.05)+
  theme_bw()+
  xlab("Var")+
  ylab("Relative Risk")+
  # scale_y_continuous(trans="log10")+
  scale_y_continuous(trans="log10", limits = c(0.85,1.15), breaks = c(0.85, 1, 1.15))+
  scale_color_manual(values = color_labels)+
  scale_fill_manual(values = color_labels)+
  scale_linetype_manual(values = linetype_labels)+
  ggtitle("PM 10")+
  labs(linetype = "additive variable", color = "additive variable", fill = "additive variable")+
  facet_wrap(.~lag_week, nrow = 1)+
  xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
  theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13))
tmp <- subset(fit_list, fit_list$cov%in%c("pm10") & fit_list$cov2%in%c("hurs","tasmax","tasmin","none") & fit_list$lag_num %in% c(0:12) & fit_list$space_time=="weekly_adm2" & fit_list$lag_num == 2)
tmp$cov2 <- factor(tmp$cov2, levels = c("hurs","tasmax","tasmin","none"))
plot_pm10_l2 <- ggplot(tmp)+
  geom_line(aes(x=var,y=rr,group=interaction(lag_week,space_time, cov2), color=cov2, linetype=cov2),alpha=0.7)+
  geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
  geom_ribbon(aes(x=var, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,space_time,cov2), fill=cov2),alpha=0.05)+
  theme_bw()+
  xlab("Var")+
  ylab("Relative Risk")+
  # scale_y_continuous(trans="log10")+
  scale_y_continuous(trans="log10", limits = c(0.9,1.1), breaks = c(0.9, 1, 1.1))+
  scale_color_manual(values = color_labels)+
  scale_fill_manual(values = color_labels)+
  scale_linetype_manual(values = linetype_labels)+
  ggtitle("PM 10")+
  labs(linetype = "additive variable", color = "additive variable", fill = "additive variable")+
  # facet_wrap(.~lag_week, nrow = 1)+
  # xlab(ax_labs[c])+
  theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
        strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13),
        legend.position = "none")

plot_hurs/plot_pm2p5/plot_pm10

plot_hurs_l2 + plot_pm2p5_l2 + plot_pm10_l2

