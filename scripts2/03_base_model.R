################################################################################
######## CREATE BASE MODEL SCRIPT
## Author: Sophie Belman
## Date: 13 October 2025
## Test random effects to estimate seasonality of pneumo in an INLA framework
## comparing monthly_district to weekly_district
################################################################################

################################################################################
####LOAD DATA##########
################################################################################
source("/home/sbelman/Documents/env_sa_manuscript/scripts2/0_source_functions.R")
# weekly=FALSE
precov=TRUE ### set whether the run ends in 2020 or 2023
time = "weekly"
space="adm2"
      ## load disease data
  data<-fread(file=paste0("/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_",space,"_",time,"_lag_sc.csv"))


## subset by only pre covid
if(precov==TRUE){
  data <- subset(data,data$date < as.Date("2020-01-01"))
  endyear = 2019
}else{
  endyear = 2023
}
      ## load shape
      # landscan_raster <- raster("/home/sbelman/Documents/BRD/SouthAfrica/sociodemographic/landscan_2022/landscan-global-2022.tif")
      ## read adjacency matrix
      if(space=="adm1"){
        g <- inla.read.graph(filename = "/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/sa_adjacency_map_adm1.adj")
        shp<-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_1.shp")
        
      }else{
        g <- inla.read.graph(filename = "/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/sa_adjacency_map.adj")
        shp<-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_2.shp")
      }
      ## load formulas
      base_form_list<-readRDS("/home/sbelman/Documents/env_sa_manuscript/formulas/base_form_list.rds")
      covs <- sapply(base_form_list, function(x) x$cov)
      #### Prepare effects to account for vaccination
      if("date"%notin%colnames(data)){
            df<- data
            df <- df %>%
              mutate(date = as.Date(week))}else{
                # mutate(date = as.Date(month))}else{
                  
                df<-data
              }
            # test interaction term
      df$date<-as.Date(df$date)
            df <- df %>%
              mutate(
                post_vaccination_2009 = ifelse(date >= as.Date("2009-01-01"), 1, 0),
                post_vaccination_2011 = ifelse(date >= as.Date("2011-01-01"), 1, 0),
                vaccination_period = ifelse(date < as.Date("2009-01-01"), 1, ifelse(date >= as.Date("2011-01-01"), 3, 2)),
                post_covid_2020 = ifelse(date >= as.Date("2020-01-01"), 1, 0),
                time_since_vaccination_2009 = as.numeric(date - as.Date("2009-01-01")) * post_vaccination_2009,
                time_since_vaccination_2011 = as.numeric(date - as.Date("2011-01-01")) * post_vaccination_2011,
                time_since_covid_2020 = as.numeric(date - as.Date("2020-01-01")) * post_covid_2020,
                PCV_coverage = PCV_coverage_3dose,
                ART_coverage = ART_coverage
              )
            
            ### include province as factors for replications
            df$id_prov <- as.numeric(factor(df$NAME_1, levels = c("Eastern_Cape", "Free_State", "Gauteng", 
                                                                  "KwaZulu-Natal", "Limpopo", "Mpumalanga", 
                                                                  "North_West", "Northern_Cape", "Western_Cape")))
            ### vaccination period as a factor
            df$vaccination_period <- as.factor(df$vaccination_period)
            df_all <- df
            # df <- subset(df,df$date <= as.Date("2021-01-01"))
            
            ### check for overdispersion
            mean_cases <- mean(df$disease, na.rm=T)
            var_cases <- var(df$disease, na.rm=T)
            var_cases/mean_cases
            
            # table(df$vaccination_period, df$year)
            
            
 ############################ run base formula used later ####
            base_form<-list()
            form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
                                  'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
                                  'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
                                  'vaccination_period', 'population_density'),
                                "disease")
            base_form$formula <- as.formula(form)
            base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for both vaccination binary, and population density")
        
            base_main <- inla.mod(base_form, fam = "nbinomial", df = df_all, nthreads=4, config=FALSE)
            saveRDS(base_main, file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_",time,"_20092011_popdens_",space,"_",endyear,".rds"))
            
          
################################################################################
####RUN MODELS#########
################################################################################
      base_mod_list<-list()
      threads=4
      for(i in 1:length(base_form_list)) {
          print("run nbinomial")
          base_mod_list[[i]] <- inla.mod(base_form_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
      }

      mod_out <- lapply(base_mod_list, function(mod) eval.mod(mod, df)) # output model evaluation metrics as a list
      mod_out <- do.call(rbind, mod_out) # collapse into df

      mod_out$model <- rownames(mod_out)
      mod_out$cov <- gsub(" ","_",mod_out$cov)
      mod_out$cov <- gsub(",","",mod_out$cov)


      #### test r2 for all mdoels
      mod_out$rsq <- apply(mod_out, 1, function(x) rsq(mod_out[mod_out$cov == x["cov"], ], base_mod_list[[1]], num_outcomes = 1))
      mod_out$num <- gsub("mae","",mod_out$model)
      mod_out$num[which(mod_out$num=="")] <- 0
      mod_out$num <- as.numeric(mod_out$num)
      mod_out$base <- paste0("base",mod_out$num)

      mod_out <- mod_out[,c("num","base","waic","mae","cpo","rsq","cov")]

     #### read and save the base model
#    saveRDS(base_mod_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_list_weekly_",endyear,".rds"))

     write.table(mod_out,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_summary_statisics_",time,"_",space,"_",endyear,".csv"),quote=FALSE,row.names=FALSE,col.names=TRUE, sep =",")
     # saveRDS(base_mod_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_list_",time,"_",space,"_",endyear,".rds"))
     saveRDS(base_mod_list[[1]],file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_main_",time,"_intercept_",space,"_",endyear,".rds"))
   
#     
##########################################################################################
####### TEST MIXED BASE MODELS WITH SOCIODEMOGRAPHIC VARS ##################      
########################################################################################## 
     ### load intercept
    # int_mod <- readRDS(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_main_",time,"_intercept_",space,"_",endyear,".rds"))
    #   base_test_list <- list()
    #   base_form<-list()
    #   form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'ART_coverage', 'vaccination_period', 'population_density'),
    #                       "disease")
    #   base_form$formula <- as.formula(form)
    #   base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for ART and both vaccines")
    #   base_test_list[[1]]<-base_form
    # 
    #   base_form<-list()
    #   form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'ART_coverage', 'post_vaccination_2009', 'population_density'),
    #                       "disease")
    #   base_form$formula <- as.formula(form)
    #   base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for ART and 2009 vaccine")
    #   base_test_list[[2]]<-base_form
    # 
    #   base_form<-list()
    #   form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'ART_coverage', 'post_vaccination_2009', 'flu_positivity', 'population_density'),
    #                       "disease")
    #   base_form$formula <- as.formula(form)
    #   base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for ART, 2009 vaccine, flu positivity")
    #   base_test_list[[3]]<-base_form
    # 
    # 
    #   base_form<-list()
    #   form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'ART_coverage', 'PCV_coverage', 'flu_positivity', 'population_density'),
    #                       "disease")
    #   base_form$formula <- as.formula(form)
    #   base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for ART, PCV_coverage, flu positivity")
    #   base_test_list[[4]]<-base_form
    # 
    #   base_form<-list()
    #   form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    #                         'vaccination_period', 'population_density'),
    #                       "disease")
    #   base_form$formula <- as.formula(form)
    #   base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for both vaccination factor, and population density")
    #   base_test_list[[5]]<-base_form
    # 
    #   base_test_models<-list()
    #   threads=4
    #   for(i in 1:length(base_test_list)) {
    #     # for(i in 5){
    #     print("run nbinomial")
    #     base_test_models[[i]] <- inla.mod(base_test_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
    #   }
    # 
    #   mod_out <- lapply(base_test_models, function(mod) eval.mod(mod, df)) # output model evaluation metrics as a list
    #   mod_out <- do.call(rbind, mod_out) # collapse into df
    # 
    #   mod_out$model <- rownames(mod_out)
    #   mod_out$cov <- gsub(" ","_",mod_out$cov)
    #   mod_out$cov <- gsub(",","",mod_out$cov)
    # 
    #   #### test r2 for all mdoels
    #   mod_out$rsq <- apply(mod_out, 1, function(x) rsq(mod_out[mod_out$cov == x["cov"], ],int_mod, num_outcomes = 1))
    #   mod_out$num <- 1:nrow(mod_out)
    #   mod_out$num <- as.numeric(mod_out$num)
    #   mod_out$base <- paste0("intervention_test")
    # 
    #   mod_out <- mod_out[,c("num","base","waic","mae","cpo","rsq","cov")]
    # 
    #   saveRDS(base_test_models[[5]], file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_",time,"_20092011_popdens_",space,"_",endyear,".rds"))
    #   write.table(mod_out, file = paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/mod_out_",time,"_intervention_test_",space,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "," )

      
##########################################################################################
####### VISUALIZE BASE MODELS ##################      
##########################################################################################  
      # library(patchwork)
      # res_vec <- c("monthly_adm1","monthly_adm2","weekly_adm1","weekly_adm2")
      # time_vec <- c("monthly","monthly","weekly","weekly")
      # space_vec <- c("adm1","adm2","adm1","adm2")
      # basestring_vec <- paste0(time_vec, "_20092011_popdens_", space_vec)
  #     
  #     ### read in a list of the base models
  #     mbase <- lapply(basestring_vec, function(re) readRDS(paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_", re, "_",endyear,".rds")))
  #     
  #     #### seasonal and annual random effects
  #     mod_idm_all <- mod_idy_all <- NULL
  #     for(re in 1:length(res_vec)){
  #       ## seasonal random effects
  #       mod_idm <- mbase[[re]]$summary.random$id_m
  #       mod_idm$model_type <- res_vec[re]
  #       mod_idm$month_name <- NULL
  #       month_names <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", 
  #                        "October", "November", "December")
  #       mod_idm$month_name <- rep(month_names, nrow(mod_idm)/12)
  #       mod_idm$month_name <- factor(mod_idm$month_name, levels = c(month_names))
  #       mod_idm_all <- rbind(mod_idm_all, mod_idm)
  #       
  #       ## annual random effects
  #       mod_idy <- mbase[[re]]$summary.random$id_y
  #       year_n <- nrow(mod_idy)/9 
  #       reg_month_year <- rep(c("Eastern_Cape", "Free_State", "Gauteng",
  #                                                      "KwaZulu-Natal", "Limpopo", "Mpumalanga",
  #                                                      "North_West", "Northern_Cape", "Western_Cape"), each = year_n)
  #       mod_idy$region <- factor(reg_month_year)
  #       mod_idy$model_type <- res_vec[re]
  #       mod_idy_all <- rbind(mod_idy_all, mod_idy)
  #     }
  #   
  #     ### plot spatial effects
  #     spat_list <- list()
  #     for(re in 1:length(res_vec)){
  #       if(space_vec[re]=="adm1"){
  #         shp<-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_1.shp")
  #       }
  #       if(space_vec[re]=="adm2"){
  #         shp<-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_2.shp")
  #       }
  #       spat_list[[re]] <- plot_spatial_effects(mbase[[re]], shp, nrow(shp), structured=FALSE, title_a = res_vec[re])
  #     }
  #     wrap_plots(spat_list)
  #     
  #   
  #     #### extract spatial effects and save
  #     re=4
  #     if(res_vec[[re]]=="weekly_adm1"){
  #       name_vec <- unique(data$adm1_name)
  #       n <- length(name_vec)
  #       full_spatial <- mbase[[3]]$summary.random$id_u[1:n,]
  #       full_spatial$adm2_name <- name_vec
  #       ## pull RR
  #       full_spatial_rr <- full_spatial
  #       full_spatial_rr[,2:6] <- exp(full_spatial[,2:6])
  #     }
  #     if(res_vec[[re]]=="weekly_adm2"){
  #       name_vec <- unique(data$adm2_name)
  #       name_grid <- unique(data[,c("NAME_1","adm2_name")])
  #       n <- length(name_vec)
  #       full_spatial <- mbase[[4]]$summary.random$id_u[1:n,]
  #       full_spatial$adm2_name <- name_vec
  #       full_spatial <- left_join(full_spatial,name_grid, by = "adm2_name")
  #       ## pull RR
  #       full_spatial_rr <- full_spatial
  #       full_spatial_rr[,2:6] <- exp(full_spatial[,2:6])
  # 
  #     }
  #     write.table(full_spatial_rr, paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/spatial_effects_",time,"_",space,"_",endyear,".csv"),
  #                 sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
  #     full_spatial_rr[order(full_spatial_rr$mean),]
  #     ## extract specific values
  #     westcape <- subset(full_spatial, full_spatial$NAME_1=="Western_Cape")
  #     
  # 
  #     
  #     ### subset the seasonal and annual by each in turn
  #     # Define ordered month names 
  #     colnames(mod_idm_all)[grep("quant",colnames(mod_idm_all))] <- c("lowerCI","median","upperCI")
  #     colnames(mod_idy_all)[grep("quant",colnames(mod_idy_all))] <- c("lowerCI","median","upperCI")
  #     
  #     mp <- subset(mod_idm_all, mod_idm_all$model_type == "weekly_adm2")
  #     mp$month_name <- factor(mp$month_name, levels = c(month_names))
  #     yp <- subset(mod_idy_all, mod_idy_all$model_type == "weekly_adm2")
  #     
  #     m <- ggplot(mp)+
  #       geom_line(aes(x = month_name, y = median, group=1))+
  #       geom_hline(yintercept=0, linetype="dashed", color="red")+
  #       geom_ribbon(aes(x = month_name, ymin = lowerCI, ymax = upperCI, group =1), alpha = 0.5)+
  #       theme_bw()+
  #       theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 13), axis.title = element_text(size = 13))+
  #       xlab("Month")+
  #       ylab("Seasonal Effect")
  #     
  #     yp$year <- 2004+ yp$ID
  #     y <- ggplot(yp)+
  #       geom_hline(yintercept=0, linetype="dashed", color="red")+
  #       geom_ribbon(aes(x = year, ymin = lowerCI, ymax = upperCI, group = region, fill = region), alpha = 0.3)+
  #       geom_line(aes(x = year, y = median, color = region, group=region))+
  #       theme_bw()+
  #       theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 13), axis.title = element_text(size = 13))+
  #       xlab("Year")+
  #       ylab("Annual Effect")+
  #       labs(fill="Province",color = "Province")+
  #       facet_wrap(region ~.)
  # 
  # pdf(paste0("/home/sbelman/Documents/BRD/SouthAfrica/manuscript/random_effects/res_weekly_adm2_",endyear,".pdf"), width = 13, height =4.5)
  #     print(      spat_list[[4]]+m+y)
  #     dev.off()
  # cbind(mp$month_name,exp(mp[,2:6]))
#       
#     
#     ## plot random effects
#     plot.fit(mod7, df)
#     plot_random_effects(mod7,mod3)
#     # plot_random_effect_provRep(mod6,mod3)
#     plot_spatial_effects(mod7, shp, 52)
# 
#     ## plot the annual random effect stratified by region for the best model (model 6) as compared to no replicate (model 3)
#     base <- mod12
#     test <- mod13
#     # load res
#     month_eff_base<-base$summary.random$id_m
#     month_eff_mod<-test$summary.random$id_m
#     year_eff_base <- base$summary.random$id_y
#     year_eff_mod <- test$summary.random$id_y
#     # name cols
#     colnames(year_eff_base)<-colnames(year_eff_mod)<-colnames(month_eff_mod)<-colnames(month_eff_base)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld")
#     cov<-c(test$cov)
#     ## number of months and years
#     # month_n<-nrow(month_eff_mod)/9
#     month_n<-nrow(month_eff_mod)/12
#     
#     year_n<-nrow(year_eff_mod)/9
#     
#     ## add region column
#     reg_month <- rep(c("Eastern Cape", "Free State", "Gauteng",
#                        "KwaZulu-Natal", "Limpopo", "Mpumalanga",
#                        "North West", "Northern Cape", "Western Cape"), each = month_n)
#     reg_month_year <- rep(c("Eastern Cape", "Free State", "Gauteng",
#                        "KwaZulu-Natal", "Limpopo", "Mpumalanga",
#                        "North West", "Northern Cape", "Western Cape"), each = year_n)
#     month_eff_mod$region<-factor(reg_month)
#     year_eff_mod$region <- factor(reg_month_year)
#     ## linetypes
#     lines<-c("Base Model"="dashed","Test"="solid")
#     ## variable
#     covs<-c(test$cov)
#     ## monthly effect (seasonal) [by province]
#     a<-ggplot()+
#       geom_line(data=month_eff_base,aes(x=ID,y=median,linetype="Base Model"))+
#       geom_line(data=month_eff_mod,aes(x=ID,y=median,linetype="Test"))+
#       geom_hline(yintercept=0,linetype="dashed",color="red",alpha=0.3)+
#       theme_classic()+
#       theme( axis.text = element_text(size=15), axis.title = element_text(size=15))+
#       labs(x = "Month",
#            y = "Median",
#            linetype = "",
#            color="Region") +
#       ylim(-1,1)+
#       labs(subtitle=covs)
#     ## regional effect (province)
#     b<-ggplot()+
#       geom_line(data=year_eff_base,aes(x=ID,y=median,linetype="Base Model"))+
#       geom_line(data=year_eff_mod,aes(x=ID,y=median,linetype="Test", color= region))+
#       # scale_linetype_manual(values = lines)+
#       geom_hline(yintercept=0,linetype="dashed",alpha=0.5,color="red")+
#       theme_classic()+
#       labs(x = "Year",
#            y = "Median",
#            linetype = "") +
#       # ylim(-2,2)+
#       ggtitle("Yearly Random Effect")+    
#       labs(subtitle=covs) +
#       facet_wrap(region~.)
#     a+b
#     ## disease cases
#     c <- ggplot(df) +
#       geom_line(aes(x=date,y=disease,group= GID_2,group=NAME_1))+
#       theme_bw()+
#       facet_wrap(NAME_1~., scales="free")
#     
#     ## disease plus random effect
#     b+c
# 
# base_mod<-base_mod_list[[3]]
# plot_random_effects(mod_list[[1]],base_mod)
# plot_spatial_effects(mod13, shp, 52)
# plot_spatial_effects(mod12, shp, 52)
