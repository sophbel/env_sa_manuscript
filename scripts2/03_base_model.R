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
time = "monthly"
space="adm1"
threads = 4
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
            
            base_intercept <- inla.mod(base_form_list[[1]], fam = "nbinomial", df = df, nthreads=4, config=FALSE)
            saveRDS(base_intercept, file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_main_",time,"_intercept_",space,"_",endyear,".rds"))
            
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
    # ## load intercept
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

      
