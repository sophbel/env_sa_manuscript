

################################################################################
#### LOAD DATA AND LIBRARIES
################################################################################
source("/home/sbelman/Documents/env_sa_manuscript/scripts2/0_source_functions.R")
# data<-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm2_weekly_lag_sc.csv")
data<-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm2_weekly_lag.csv")
shp<-st_read("/home/sbelman/Documents/BRD/SouthAfrica/shps/gadm41_namematch_ZAF_2.shp")
## read adjacency matrix
g <- inla.read.graph(filename = "/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/sa_adjacency_map.adj")
base_intercept <- readRDS(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_main_weekly_intercept_adm2_2019.rds"))
precov = TRUE
################################################################################
#### PREPARE DATA
################################################################################
df<- data
df <- df %>%
  mutate(date = as.Date(date))
## create interaction terms
df <- df %>%
  mutate(
    vaccination_period = ifelse(date < as.Date("2009-01-01"), 1, ifelse(date >= as.Date("2011-01-01"), 3, 2)),
    post_vaccination_2009 = ifelse(date >= as.Date("2009-01-01"), 1, 0),
    post_vaccination_2011 = ifelse(date >= as.Date("2011-01-01"), 1, 0),
    post_covid_2020 = ifelse(date >= as.Date("2020-01-01"), 1, 0),
    time_since_vaccination_2009 = as.numeric(date - as.Date("2009-01-01")) * post_vaccination_2009,
    time_since_vaccination_2011 = as.numeric(date - as.Date("2011-01-01")) * post_vaccination_2011,
    time_since_covid_2020 = as.numeric(date - as.Date("2020-01-01")) * post_covid_2020,
    # PCV_coverage = PCV_coverage_3dose*100,
    # ART_coverage = ART_coverage*100
    PCV_coverage = PCV_coverage_3dose,
    ART_coverage = ART_coverage
  )

## vectors with variable names
all_names<-colnames(df)
### environmental
env_vars<-c(grep("sp",all_names,value=TRUE),grep("tas",all_names,value=TRUE),
            grep("hurs",all_names,value=TRUE),
            grep("absh",all_names,value=TRUE),grep("prlr",all_names,value=TRUE),grep("pm",all_names,value=TRUE),grep("o3",all_names,value=TRUE))
gene_vars<-grep("monthly_count",all_names,value=TRUE)

### include province as factors for replications
df$id_prov <- as.numeric(factor(df$NAME_1))

### group the ART coverage nationally 
art_national <- df %>%
  group_by(date) %>%
  summarise(ART_coverage_national = mean(ART_coverage, na.rm=T))

df <- left_join(df, art_national, by="date")

### create vt column
df <- df %>%
  mutate(vts = pcv7 + pcv13)
df$vaccination_period <- as.factor(df$vaccination_period)

## subset by only pre covid
if(precov==TRUE){
  df <- subset(df,df$date < as.Date("2020-01-01"))
  endyear = 2019
}else{
  endyear = 2023
}

################################################################################
#### write a formula list ####
################################################################################
base_pcvartprovrep_formula_list <- base_pcvart_formula_list <- list()
base_form<-list()
form <- as.formula(disease ~ 1)
base_form$formula <- as.formula(form)
base_form$covs<-paste0("Intercept")
base_pcvart_formula_list[[1]]<-base_form

######### base models for comparison in this list #######
base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual")
base_pcvart_formula_list[[2]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication)")
base_pcvart_formula_list[[3]]<-base_form

######### include coverage of ART and PCV #########
base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")',
                      'ART_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for ART coverage")
base_pcvart_formula_list[[4]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")',
                      'PCV_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for PCV coverage")
base_pcvart_formula_list[[5]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")',
                      'vaccination_period'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for vaccination_period")
base_pcvart_formula_list[[6]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")',
                      'post_vaccination_2009'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for 2009 vaccine")
base_pcvart_formula_list[[7]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")',
                      'ART_coverage_national'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for ART coverage Nationally")
base_pcvart_formula_list[[8]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")',
                      'vaccination_period',
                      'ART_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for vaccine period and ART coverage")
base_pcvart_formula_list[[9]]<-base_form

### include province replication ###
base_pcvartprovrep_formula_list <- list()

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)',
                      'ART_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for ART coverage")
base_pcvartprovrep_formula_list[[1]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)',
                      'PCV_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for PCV coverage")
base_pcvartprovrep_formula_list[[2]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)',
                      'vaccination_period'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for vaccination_period")
base_pcvartprovrep_formula_list[[3]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)',
                      'post_vaccination_2009'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for 2009 vaccine")
base_pcvartprovrep_formula_list[[4]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)',
                      'ART_coverage_national'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for ART coverage Nationally")
base_pcvartprovrep_formula_list[[5]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)',
                      'vaccination_period',
                        'ART_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for vaccination_period and ART coverage")
base_pcvartprovrep_formula_list[[6]]<-base_form



########## use function from source functions here ############
##### include non-vaccine types only and no province replication##### ##### ##### 
base_pcvart_nvt_formula_list <- perturbation_formula_list(response_var = "nvt", replicate = FALSE)
##### include vaccine types only and no province replication##### ##### ##### 
base_pcvart_vts_formula_list <- perturbation_formula_list(response_var = "vts", replicate = FALSE)
##### adults only and no province replication##### ##### ##### 
base_pcvart_adults_formula_list <- perturbation_formula_list(response_var = "age_15t64", replicate = FALSE)
##### children only and no province replication##### ##### ##### 
base_pcvart_children_formula_list <- perturbation_formula_list(response_var = "age_lt6", replicate = FALSE)
##### include non-vaccine types only and with province replication##### ##### ##### 
base_pcvartprovrep_nvt_formula_list <- perturbation_formula_list(response_var = "nvt", replicate = TRUE)
##### include vaccine types only and with province replication##### ##### ##### 
base_pcvartprovrep_vts_formula_list <- perturbation_formula_list(response_var = "vts", replicate = TRUE)
##### adults only and with province replication##### ##### ##### 
base_pcvartprovrep_adults_formula_list <- perturbation_formula_list(response_var = "age_15t64", replicate = TRUE)
##### children only and with province replication##### ##### ##### 
base_pcvartprovrep_children_formula_list <- perturbation_formula_list(response_var = "age_lt6", replicate = TRUE)

################################################################################
#### RUN PCV ART MODELS
################################################################################
##################### with no province replication ###################
# ##### no province replication
# base_pcvart_model_list <-list()
# threads=4
# for(i in 1:length(base_pcvart_formula_list)) {
#   print(i)
#   print("run nbinomial")
#   base_pcvart_model_list[[i]] <- inla.mod(base_pcvart_formula_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
# }
# 
# saveRDS(base_pcvart_model_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvart_model_list_weekly_",endyear,".rds"))


# #################### with province replication ###################
# base_pcvartprovrep_model_list <-list()
# threads=4
# for(i in 1:length(base_pcvartprovrep_formula_list)) {
# 
#   print(i)
#   print("run nbinomial")
#   base_pcvartprovrep_model_list[[i]] <- inla.mod(base_pcvartprovrep_formula_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
# }
# saveRDS(base_pcvartprovrep_model_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/base_pcvartprovrep_model_list_weekly_",endyear,".rds"))

# base_pcvart_model_list <- readRDS(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvart_model_list_weekly_",endyear,".rds"))
mod_out <- evaluate_model_list(base_pcvart_model_list, df, base_intercept, num_outcomes = 1)
mod_out$outcome <- "disease"
base_pcvartprovrep_model_list <- readRDS(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/base_pcvartprovrep_model_list_weekly_",endyear,".rds"))
mod_out_pr <- evaluate_model_list(base_pcvartprovrep_model_list, df, base_intercept, num_outcomes = 1)
mod_out_pr$outcome <- "disease"
### save the base models
mod_out_bases <- rbind(mod_out, mod_out_pr)
write.table(mod_out_bases,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvart_summary_statisics_weekly_",endyear,".csv"),quote=FALSE,row.names=FALSE,col.names=TRUE)


##################### NO PROVINCE REPLICATION - SUBSET MODELS ###################
        ### with no province replication but just non-vaccine types
        base_pcvart_nvtmodel_list <-list()
        threads=4
        for(i in 1:length(base_pcvart_nvt_formula_list)) {
          print(i)
          print("run nbinomial")
          base_pcvart_nvtmodel_list[[i]] <- inla.mod(base_pcvart_nvt_formula_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
        }      
        saveRDS(base_pcvart_nvtmodel_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvart_nvtmodel_list_weekly_",endyear,".rds"))
        mod_out_nvt <- evaluate_model_list(base_pcvart_nvtmodel_list, df, base_intercept, num_outcomes = 1)
        mod_out_nvt$outcome <- "nvts"
        ### extract fixed effects and save
              ### fixed effect models no province replicates non-vaccine types alone
              modART_nvt <- base_pcvart_nvtmodel_list[[1]]
              modPCV_nvt <- base_pcvart_nvtmodel_list[[2]]
              mod_vaxperiod_nvt <- base_pcvart_nvtmodel_list[[3]]
              mod_2009_nvt <- base_pcvart_nvtmodel_list[[4]]
              mod_ARTnational_nvt <- base_pcvart_nvtmodel_list[[5]]
              mod_vaxART_nvt <- base_pcvart_nvtmodel_list[[6]]
              ### fixed table
              fixed_effects_nvts <- rbind(modART_nvt$summary.fixed[2,],
                                          modPCV_nvt$summary.fixed[2,],
                                          mod_vaxperiod_nvt$summary.fixed[2,],
                                          mod_2009_nvt$summary.fixed[2,],
                                          mod_ARTnational_nvt$summary.fixed[2,],
                                          mod_vaxART_nvt$summary.fixed[2:3,])
              fixed_effects_nvts$outcome <- "nvts"
        
        
        ### with no province replication but just vaccine types
        base_pcvart_vtsmodel_list <-list()
        threads=4
        for(i in 1:length(base_pcvart_vts_formula_list)) {
          
          print(i)
          print("run nbinomial")
          base_pcvart_vtsmodel_list[[i]] <- inla.mod(base_pcvart_vts_formula_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
        }      
        saveRDS(base_pcvart_vtsmodel_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvart_vtsmodel_list_weekly_",endyear,".rds"))
        mod_out_vts <- evaluate_model_list(base_pcvart_vtsmodel_list, df, base_intercept, num_outcomes = 1)
        mod_out_vts$outcome <- "vts"
        ### extract fixed effects and save 
              ### fixed effect models no province replicates vaccine types alone
              modART_vts <- base_pcvart_vtsmodel_list[[1]]
              modPCV_vts <- base_pcvart_vtsmodel_list[[2]]
              mod_vaxperiod_vts <- base_pcvart_vtsmodel_list[[3]]
              mod_2009_vts <- base_pcvart_vtsmodel_list[[4]]
              mod_ARTnational_vts <- base_pcvart_vtsmodel_list[[5]]
              mod_vaxART_vts <- base_pcvart_vtsmodel_list[[6]]
              ## fixed
              fixed_effects_vts <- rbind(modART_vts$summary.fixed[2,],
                                         modPCV_vts$summary.fixed[2,],
                                         mod_vaxperiod_vts$summary.fixed[2,],
                                         mod_2009_vts$summary.fixed[2,],
                                         mod_ARTnational_vts$summary.fixed[2,],
                                         mod_vaxART_vts$summary.fixed[2:3,])
              fixed_effects_vts$outcome <- "vts"
        
        ### with no province replication but just adults
              base_pcvart_adultmodel_list <-list()
        threads=4
        for(i in 1:length(base_pcvart_adults_formula_list)) {
          
          print(i)
          print("run nbinomial")
          base_pcvart_adultmodel_list[[i]] <- inla.mod(base_pcvart_adults_formula_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
        }      
        saveRDS(base_pcvart_adultmodel_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvart_adultmodel_list_weekly_",endyear,".rds"))
        mod_out_adults <- evaluate_model_list(base_pcvart_adultmodel_list, df, base_intercept, num_outcomes = 1)
        mod_out_adults$outcome <- "adults"
        ### extract fixed effects and save
            ### fixed effect models no province replicates adults
            modART_adults <- base_pcvart_adultmodel_list[[1]]
            modPCV_adults <- base_pcvart_adultmodel_list[[2]]
            mod_vaxperiod_adults <- base_pcvart_adultmodel_list[[3]]
            mod_2009_adults <- base_pcvart_adultmodel_list[[4]]
            mod_ARTnational_adults <- base_pcvart_adultmodel_list[[5]]
            mod_vaxART_adults <- base_pcvart_adultmodel_list[[6]]
            ### fixed table
            fixed_effects_adults <- rbind(modART_adults$summary.fixed[2,],
                                          modPCV_adults$summary.fixed[2,],
                                          mod_vaxperiod_adults$summary.fixed[2,],
                                          mod_2009_adults$summary.fixed[2,],
                                          mod_ARTnational_adults$summary.fixed[2,],
                                          mod_vaxART_adults$summary.fixed[2:3,])
            fixed_effects_adults$outcome <- "adults"
        
        ### with no province replication but just children
        base_pcvart_childrenmodel_list <-list()
        threads=4
        for(i in 1:length(base_pcvart_children_formula_list)) {
          
          print(i)
          print("run nbinomial")
          base_pcvart_childrenmodel_list[[i]] <- inla.mod(base_pcvart_children_formula_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
        }      
        saveRDS(base_pcvart_childrenmodel_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvart_childrenmodel_list_weekly_",endyear,".rds"))
        mod_out_children <- evaluate_model_list(base_pcvart_childrenmodel_list, df, base_intercept, num_outcomes = 1)
        mod_out_children$outcome <- "children"
        ### extract fixed effects and save
              ### fixed effect models no province replicates children
              modART_children <- base_pcvart_childrenmodel_list[[1]]
              modPCV_children <- base_pcvart_childrenmodel_list[[2]]
              mod_vaxperiod_children <- base_pcvart_childrenmodel_list[[3]]
              mod_2009_children <- base_pcvart_childrenmodel_list[[4]]
              mod_ARTnational_children <- base_pcvart_childrenmodel_list[[5]]
              mod_vaxART_children <- base_pcvart_childrenmodel_list[[6]]
              ### fixed table
              fixed_effects_children <- rbind(modART_children$summary.fixed[2,],
                                              modPCV_children$summary.fixed[2,],
                                              mod_vaxperiod_children$summary.fixed[2,],
                                              mod_2009_children$summary.fixed[2,],
                                              mod_ARTnational_children$summary.fixed[2,],
                                              mod_vaxART_children$summary.fixed[2:3,])
              fixed_effects_children$outcome <- "children"
        
        ### GOFS join  together ###
        mod_out <- rbind(mod_out_nvt, mod_out_vts, mod_out_adults, mod_out_children)
        ## save table
        mod_out_tmp <- mod_out[,c("num","base","waic","mae","cpo","rsq","cov")]
        write.table(mod_out,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvart_variousoutcomes_summary_statisics_weekly_",endyear,".csv"),quote=FALSE,row.names=FALSE,col.names=TRUE)
        
        
        ########## LOAD ORIGINAL MODELS AGAIN ####
        # ##### load only the models we need
        base_pcvart_model_list <- readRDS(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvart_model_list_weekly_",endyear,".rds"))
        covs <- sapply(base_pcvart_model_list, function(x) x$cov)
        
        ### base models
        mod_base <- base_pcvart_model_list[[2]]
        mod_base_pr <- base_pcvart_model_list[[3]]
        
        ### fixed effect models no province replicates
        modART <- base_pcvart_model_list[[4]]
        modPCV <- base_pcvart_model_list[[5]]
        mod_vaxperiod <- base_pcvart_model_list[[6]]
        mod_2009 <- base_pcvart_model_list[[7]]
        mod_ARTnational <- base_pcvart_model_list[[8]]
        mod_vaxART <- base_pcvart_model_list[[9]]
        ################ FIXED EFFECTS LINEAR - NO PROVINCE REPLICATION ##############
        fixed_effects <- rbind(modART$summary.fixed[2,],
                               modPCV$summary.fixed[2,],
                               mod_vaxperiod$summary.fixed[2,],
                               mod_2009$summary.fixed[2,],
                               mod_ARTnational$summary.fixed[2,],
                               mod_vaxART$summary.fixed[2:3,])
        fixed_effects$outcome <- "disease"
        
        
        ### fixed effects join together ###
        fixed_effects_all <- rbind(fixed_effects, fixed_effects_nvts, fixed_effects_vts, fixed_effects_adults, fixed_effects_children)
        write.table(fixed_effects_all,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvart_variousoutcomes_fixedeffects_adm2_weekly_",endyear,".csv"),quote=FALSE,row.names=FALSE,col.names=TRUE)

##################### WITH PROVINCE REPLICATION - SUBSET MODELS ###################
            ### with no province replication but just non-vaccine types
            base_pcvartprovrep_nvtmodel_list <-list()
            threads=4
            for(i in 1:length(base_pcvartprovrep_nvt_formula_list)) {
              print(i)
              print("run nbinomial")
              base_pcvartprovrep_nvtmodel_list[[i]] <- inla.mod(base_pcvartprovrep_nvt_formula_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
            }      
            saveRDS(base_pcvartprovrep_nvtmodel_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/base_pcvartprovrep_nvtmodel_list_weekly_",endyear,".rds"))
            mod_out_nvt <- evaluate_model_list(base_pcvartprovrep_nvtmodel_list, df, base_intercept, num_outcomes = 1)
            mod_out_nvt$outcome <- "nvts"
            ### extract fixed effects and save
            ### fixed effect models no province replicates non-vaccine types alone
            modART_nvt_pr <- base_pcvartprovrep_nvtmodel_list[[1]]
            modPCV_nvt_pr <- base_pcvartprovrep_nvtmodel_list[[2]]
            mod_vaxperiod_nvt_pr <- base_pcvartprovrep_nvtmodel_list[[3]]
            mod_2009_nvt_pr <- base_pcvartprovrep_nvtmodel_list[[4]]
            mod_ARTnational_nvt_pr <- base_pcvartprovrep_nvtmodel_list[[5]]
            mod_vaxART_nvt_pr <- base_pcvartprovrep_nvtmodel_list[[6]]
            ### fixed table
            fixed_effects_nvts_pr <- rbind(modART_nvt_pr$summary.fixed[2,],
                                        modPCV_nvt_pr$summary.fixed[2,],
                                        mod_vaxperiod_nvt_pr$summary.fixed[2,],
                                        mod_2009_nvt_pr$summary.fixed[2,],
                                        mod_ARTnational_nvt_pr$summary.fixed[2,],
                                        mod_vaxART_nvt_pr$summary.fixed[2:3,])
            fixed_effects_nvts_pr$outcome <- "nvts"
            
            
            ### with no province replication but just vaccine types
            base_pcvartprovrep_vtsmodel_list <-list()
            threads=4
            for(i in 1:length(base_pcvartprovrep_vts_formula_list)) {
              
              print(i)
              print("run nbinomial")
              base_pcvartprovrep_vtsmodel_list[[i]] <- inla.mod(base_pcvartprovrep_vts_formula_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
            }      
            saveRDS(base_pcvartprovrep_vtsmodel_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/base_pcvartprovrep_vtsmodel_list_weekly_",endyear,".rds"))
            mod_out_vts <- evaluate_model_list(base_pcvartprovrep_vtsmodel_list, df, base_intercept, num_outcomes = 1)
            mod_out_vts$outcome <- "vts"
            ### extract fixed effects and save 
            ### fixed effect models no province replicates vaccine types alone
            modART_vts_pr <- base_pcvartprovrep_vtsmodel_list[[1]]
            modPCV_vts_pr <- base_pcvartprovrep_vtsmodel_list[[2]]
            mod_vaxperiod_vts_pr <- base_pcvartprovrep_vtsmodel_list[[3]]
            mod_2009_vts_pr <- base_pcvartprovrep_vtsmodel_list[[4]]
            mod_ARTnational_vts_pr <- base_pcvartprovrep_vtsmodel_list[[5]]
            mod_vaxART_vts_pr <- base_pcvartprovrep_vtsmodel_list[[6]]
            ## fixed
            fixed_effects_vts_pr <- rbind(modART_vts_pr$summary.fixed[2,],
                                       modPCV_vts_pr$summary.fixed[2,],
                                       mod_vaxperiod_vts_pr$summary.fixed[2,],
                                       mod_2009_vts_pr$summary.fixed[2,],
                                       mod_ARTnational_vts_pr$summary.fixed[2,],
                                       mod_vaxART_vts_pr$summary.fixed[2:3,])
            fixed_effects_vts_pr$outcome <- "vts"
            
            ### with no province replication but just adults
            base_pcvartprovrep_adultmodel_list <-list()
            threads=4
            for(i in 1:length(base_pcvartprovrep_adults_formula_list)) {
              
              print(i)
              print("run nbinomial")
              base_pcvartprovrep_adultmodel_list[[i]] <- inla.mod(base_pcvartprovrep_adults_formula_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
            }      
            saveRDS(base_pcvartprovrep_adultmodel_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/base_pcvartprovrep_adultmodel_list_weekly_",endyear,".rds"))
            mod_out_adults <- evaluate_model_list(base_pcvartprovrep_adultmodel_list, df, base_intercept, num_outcomes = 1)
            mod_out_adults$outcome <- "adults"
            ### extract fixed effects and save
            ### fixed effect models no province replicates adults
            modART_adults_pr <- base_pcvartprovrep_adultmodel_list[[1]]
            modPCV_adults_pr <- base_pcvartprovrep_adultmodel_list[[2]]
            mod_vaxperiod_adults_pr <- base_pcvartprovrep_adultmodel_list[[3]]
            mod_2009_adults_pr <- base_pcvartprovrep_adultmodel_list[[4]]
            mod_ARTnational_adults_pr <- base_pcvartprovrep_adultmodel_list[[5]]
            mod_vaxART_adults_pr <- base_pcvartprovrep_adultmodel_list[[6]]
            ### fixed table
            fixed_effects_adults_pr <- rbind(modART_adults_pr$summary.fixed[2,],
                                          modPCV_adults_pr$summary.fixed[2,],
                                          mod_vaxperiod_adults_pr$summary.fixed[2,],
                                          mod_2009_adults_pr$summary.fixed[2,],
                                          mod_ARTnational_adults_pr$summary.fixed[2,],
                                          mod_vaxART_adults_pr$summary.fixed[2:3,])
            fixed_effects_adults_pr$outcome <- "adults"
            
            ### with no province replication but just children
            base_pcvartprovrep_childrenmodel_list <-list()
            threads=4
            for(i in 1:length(base_pcvartprovrep_children_formula_list)) {
              
              print(i)
              print("run nbinomial")
              base_pcvartprovrep_childrenmodel_list[[i]] <- inla.mod(base_pcvartprovrep_children_formula_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
            }      
            saveRDS(base_pcvartprovrep_childrenmodel_list,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/base_pcvartprovrep_childrenmodel_list_weekly_",endyear,".rds"))
            mod_out_children <- evaluate_model_list(base_pcvartprovrep_childrenmodel_list, df, base_intercept, num_outcomes = 1)
            mod_out_children$outcome <- "children"
            ### extract fixed effects and save
            ### fixed effect models no province replicates children
            modART_children_pr <- base_pcvartprovrep_childrenmodel_list[[1]]
            modPCV_children_pr <- base_pcvartprovrep_childrenmodel_list[[2]]
            mod_vaxperiod_children_pr <- base_pcvartprovrep_childrenmodel_list[[3]]
            mod_2009_children_pr <- base_pcvartprovrep_childrenmodel_list[[4]]
            mod_ARTnational_children_pr <- base_pcvartprovrep_childrenmodel_list[[5]]
            mod_vaxART_children_pr <- base_pcvartprovrep_childrenmodel_list[[6]]
            ### fixed table
            fixed_effects_children_pr <- rbind(modART_children_pr$summary.fixed[2,],
                                            modPCV_children_pr$summary.fixed[2,],
                                            mod_vaxperiod_children_pr$summary.fixed[2,],
                                            mod_2009_children_pr$summary.fixed[2,],
                                            mod_ARTnational_children_pr$summary.fixed[2,],
                                            mod_vaxART_children_pr$summary.fixed[2:3,])
            fixed_effects_children_pr$outcome <- "children"
            
            ### GOFS join  together ###
            mod_out <- rbind(mod_out_nvt, mod_out_vts, mod_out_adults, mod_out_children)
            ## save table
            mod_out_tmp <- mod_out[,c("num","base","waic","mae","cpo","rsq","cov")]
            write.table(mod_out,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/base_pcvartprovrep_variousoutcomes_summary_statisics_weekly_",endyear,".csv"),quote=FALSE,row.names=FALSE,col.names=TRUE)
            
            
            
            ########## LOAD ORIGINAL MODELS AGAIN ####
            # ##### load only the models we need
            base_pcvart_model_list <- readRDS(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvart_model_list_weekly_",endyear,".rds"))
            base_pcvartprovrep_model_list <- readRDS(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/base_pcvartprovrep_model_list_weekly_",endyear,".rds"))
            covs <- sapply(base_pcvartprovrep_model_list, function(x) x$cov)
            
            ### base models
            mod_base <- base_pcvart_model_list[[2]]
            mod_base_pr <- base_pcvart_model_list[[3]]
            # saveRDS(mod_base_pr, file = paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/base_pcvartprovrep_noperturbations_weekly_",endyear,".rds"))
            
            ### fixed effect models no province replicates
            modART_pr <- base_pcvartprovrep_model_list[[1]]
            modPCV_pr <- base_pcvartprovrep_model_list[[2]]
            mod_vaxperiod_pr <- base_pcvartprovrep_model_list[[3]]
            mod_2009_pr <- base_pcvartprovrep_model_list[[4]]
            mod_ARTnational_pr <- base_pcvartprovrep_model_list[[5]]
            mod_vaxART_pr <- base_pcvartprovrep_model_list[[6]]
            ################ FIXED EFFECTS LINEAR - WITH PROVINCE REPLICATION ##############
            fixed_effects_pr <- rbind(modART_pr$summary.fixed[2,],
                                   modPCV_pr$summary.fixed[2,],
                                   mod_vaxperiod_pr$summary.fixed[2,],
                                   mod_2009_pr$summary.fixed[2,],
                                   mod_ARTnational_pr$summary.fixed[2,],
                                   mod_vaxART_pr$summary.fixed[2:3,])
            fixed_effects_pr$outcome <- "disease"
       
            ### fixed effects join together ###
            fixed_effects_all_pr <- rbind(fixed_effects_pr, fixed_effects_nvts_pr, fixed_effects_vts_pr, fixed_effects_adults_pr, fixed_effects_children_pr)
            write.table(fixed_effects_all_pr,file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvartprovrep_variousoutcomes_fixedeffects_adm2_weekly_",endyear,".csv"),quote=FALSE,row.names=FALSE,col.names=TRUE)

#################################################################################
# #### PLOT COMPARISONS OF DIFFERENT OUTCOMES AND PERTURBATIONS
#################################################################################
            endyear = 2023
            fixed_effects_all_pr <- fread(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/base_pcvartprovrep_variousoutcomes_fixedeffects_adm2_weekly_",endyear,".csv"))
            fixed_effects_all <- fread(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvart_models/base_pcvart_variousoutcomes_fixedeffects_adm2_weekly_",endyear,".csv"))

#################### NO PROVINCE REPLICATION FIXED EFFECTS #########
fixedeffs <- cbind(fixed_effects_all, model = rep(c("ART","PCV", "vaccine_periods", "2009vaccine","ARTnational","ART_Vax","ART_Vax"),5), 
                   effect = rep(c("ART","PCV","vaccine_periods","2009vaccine","ARTnational","vaccine_periods","ART"),5))
fixedeffs$model <- factor(fixedeffs$model, levels = c("ART_Vax","ART","ARTnational","2009vaccine","vaccine_periods","PCV"))
colnames(fixedeffs) <- c("mean","sd","lowerCI","median","upperCI","mode","kld","outcome","model","effect")
color_labels <- c("ART" = "#097969","PCV" = "#D91656", "ARTnational" = "grey", "ART_Vax" = "darkblue", "vaccine_periods"="orange","2009vaccine"="purple" )

fixed_national_all <- ggplot(fixedeffs)+
  geom_hline(yintercept=0,linetype="dashed",color="red")+
  geom_point(aes(x=effect,y=median, group=model, color= model),size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x=effect,ymin=lowerCI,ymax=upperCI, group=model, color= model), width=.3, position=position_dodge(width=0.5))+
  theme_bw()+
  xlab("Intervention")+
  ggtitle(paste0("2005-",endyear)) + 
  # ylim(-4,0.1)+
  ylab("Effect")+
  theme(axis.text.x = element_text(angle=45,hjust=1, size=13), axis.text = element_text(size=13), axis.title = element_text(size=13)) +
  scale_color_manual(values=color_labels)+
  labs(color="Model")+
  facet_wrap(outcome ~.)

unique(fixedeffs$outcome)
fixedeffs$outcome <- factor(fixedeffs$outcome,
                               levels = c("disease", "adults", "children", "vts", "nvts"), 
                               labels = c("disease", "adult (15-64)", "children (≤5)", "VT", "NVT"))

fixed_national_all_model <- ggplot(fixedeffs)+
  geom_hline(yintercept=0,linetype="dashed",color="red")+
  geom_point(aes(x=outcome,y=median, group=model, color= effect),size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x=outcome,ymin=lowerCI,ymax=upperCI, group=model, color= effect), width=.3, position=position_dodge(width=0.5))+
  theme_bw()+
  xlab("Outcome")+
  ggtitle(paste0("2005-",endyear)) + 
  # ggtitle("Province Replicated Fixed Effects") + 
  # ylim(-4,0.1)+
  ylab("Effect")+
  theme(axis.text.x = element_text(angle=45,hjust=1, size=13), axis.text = element_text(size=13), axis.title = element_text(size=13)) +
  scale_color_manual(values=color_labels)+
  labs(color="Model")+
  facet_wrap(model ~.)


pdf(paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/fixed_effect_perturbations_",endyear,".pdf"), width = 8, height =5)
print(fixed_national_all_model)
dev.off()

print(fixed_national_all_model)
ggsave(paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/fixed_effect_perturbations_",endyear,".png"), width = 8, height =5)
dev.off()


#################### WITH PROVINCE REPLICATION FIXED EFFECTS #########
fixedeffs_pr <- cbind(fixed_effects_all_pr, model = rep(c("ART","PCV", "vaccine_periods", "2009vaccine","ARTnational","ART_Vax","ART_Vax"),5), 
                   effect = rep(c("ART","PCV","vaccine_periods","2009vaccine","ARTnational","vaccine_periods","ART"),5))
fixedeffs_pr$model <- factor(fixedeffs_pr$model, levels = c("ART_Vax","ART","ARTnational","2009vaccine","vaccine_periods","PCV"))
colnames(fixedeffs_pr) <- c("mean","sd","lowerCI","median","upperCI","mode","kld","outcome","model","effect")
color_labels <- c("ART" = "#097969","PCV" = "#D91656", "ARTnational" = "grey", "ART_Vax" = "darkblue", "vaccine_periods"="orange","2009vaccine"="purple" )
fixedeffs_pr$effect <- factor(fixedeffs_pr$effect, levels = c("ART_Vax","ART","ARTnational","2009vaccine","vaccine_periods","PCV"))

fixedeffs_pr$outcome <- factor(fixedeffs_pr$outcome, levels = c("disease","adults","children","vts","nvts"))
fixed_national_all_pr <- ggplot(fixedeffs_pr)+
  geom_hline(yintercept=0,linetype="dashed",color="red")+
  geom_point(aes(x=effect,y=median, group=model, color= model ),size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x=effect,ymin=lowerCI,ymax=upperCI, group=model, color= model), width=.3, position=position_dodge(width=0.5))+
  theme_bw()+
  xlab("Intervention")+
  ggtitle(paste0("Province Replicated 2005-",endyear)) + 
  # ylim(-4,0.1)+
  ylab("Effect")+
  theme(axis.text.x = element_text(angle=45,hjust=1, size=13), axis.text = element_text(size=13), axis.title = element_text(size=13)) +
  scale_color_manual(values=color_labels)+
  labs(color="Model")+
  facet_wrap(outcome ~.)

unique(fixedeffs_pr$outcome)
fixedeffs_pr$outcome <- factor(fixedeffs_pr$outcome,
                               levels = c("disease", "adults", "children", "vts", "nvts"), 
                               labels = c("disease", "adult (15-64)", "children (≤5)", "VT", "NVT"))

fixed_national_all_pr_model <- ggplot(fixedeffs_pr)+
  geom_hline(yintercept=0,linetype="dashed",color="red")+
  geom_point(aes(x=outcome,y=median, group=model, color= effect),size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x=outcome,ymin=lowerCI,ymax=upperCI, group=model, color= effect), width=.3, position=position_dodge(width=0.5))+
  theme_bw()+
  xlab("Outcome")+
  ggtitle(paste0("Province Replicated 2005-",endyear)) + 
  # ggtitle("Province Replicated Fixed Effects") + 
  # ylim(-4,0.1)+
  ylab("Effect")+
  theme(axis.text.x = element_text(angle=45,hjust=1, size=13), axis.text = element_text(size=13), axis.title = element_text(size=13)) +
  scale_color_manual(values=color_labels)+
  labs(color="Model")+
  facet_wrap(model ~.)

### save
pdf(paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/fixed_effect_provrep_perturbations_",endyear,".pdf"), width = 8, height =5)
print(fixed_national_all_pr_model)
dev.off()
print(fixed_national_all_pr_model)
ggsave(paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/fixed_effect_provrep_perturbations_",endyear,".png"), width = 8, height =5)
dev.off()

# fixedeffs_pr[which(fixedeffs_pr$outcome=="NVT"),]
################################ VIF ###########################################

##### calculate Variance inflation factor between vaccine_periods and NVT coverage nationally
base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)',
                      'PCV_coverage',
                      'ART_coverage_national'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for vaccination_period and ART coverage nationally")
mod_vif <- inla.mod(base_form, fam = "nbinomial", df = df, nthreads=4, config=TRUE)



mod <- mod_vif  # your fitted INLA model
fixed_names <- rownames(mod$summary.fixed)[2:3]
n_fixed <- length(fixed_names)
# ️ Compute raw correlations and classical VIF
# extract predictor columns from your data
df2 <- data.frame(df)
X <- df2[, fixed_names]

# raw correlation
raw_cor <- cor(X)

# classical VIF
# VIF_j = 1 / (1 - R^2_j)
classic_vif <- numeric(n_fixed)
names(classic_vif) <- fixed_names

for(j in 1:n_fixed){
  y <- X[, j]
  X_others <- X[, -j, drop=FALSE]
  R2 <- summary(lm(y ~ ., data = as.data.frame(X_others)))$r.squared
  classic_vif[j] <- 1 / (1 - R2)
}


### plot for supplement######
### fixed effect models with province replicates
modART_pr <- base_pcvartprovrep_model_list[[1]]
modPCV_pr <- base_pcvartprovrep_model_list[[2]]
mod_vaxperiod_pr <- base_pcvartprovrep_model_list[[3]]
mod_2009_pr <- base_pcvartprovrep_model_list[[4]]
mod_ARTnational_pr <- base_pcvartprovrep_model_list[[5]]
mod_vaxART_pr <- base_pcvartprovrep_model_list[[6]]

baseidy <- mod_base_pr$summary.random$id_y
vaxperiodidy <- mod_vaxperiod_pr$summary.random$id_y
ARTidy <- modART_pr$summary.random$id_y
vaxARTidy <- mod_vaxART_pr$summary.random$id_y

colnames(vaxARTidy)<-colnames(ARTidy)<-colnames(vaxperiodidy)<-colnames(baseidy)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld")
## number of months and years
month_n<-length(unique(baseidy$ID))

## add region column
reg_month <- rep(c("Eastern Cape", "Free State", "Gauteng", 
                   "KwaZulu-Natal", "Limpopo", "Mpumalanga", 
                   "North West", "Northern Cape", "Western Cape"), each = month_n)
vaxARTidy$region<-ARTidy$region<-vaxperiodidy$region<-baseidy$region<-factor(reg_month)
baseidy$type <- "base"
vaxperiodidy$type <- "vaccine_periods"
ARTidy$type <- "ART"
vaxARTidy$type <- "ART_Vax"

## linetypes
lines<-c("base"="dashed", "ART" = "solid", "ART_Vax" = "solid", "vaccine_periods"="solid")
year_randomeffects <- rbind(vaxperiodidy,ARTidy,vaxARTidy,baseidy)
year_randomeffects$type <- factor(year_randomeffects$type, levels = c("base","vaccine_periods","ART","ART_Vax"))
year_randomeffects$year <- 2004 + year_randomeffects$ID
## 
color_labels <- c("base" = "black","ART" = "#097969","PCV" = "#D91656", "ART_Vax" = "darkblue", "vaccine_periods"="orange","2009vaccine"="purple" )

a<-ggplot()+
  # geom_line(data=baseidy,aes(x=ID,y=median, group=region),linetype="dashed")+
  geom_line(data=year_randomeffects,aes(x=year,y=median,group= interaction(region,type), color = type, linetype = type))+
  geom_hline(yintercept=0,color="darkred",alpha=0.8)+
  theme_classic()+
  theme( axis.text = element_text(size=15), axis.title = element_text(size=15),
         axis.text.x = element_text(angle=45, hjust =1))+
  labs(x = "Year",
       y = "Effect",
       linetype = "Model",
       color="Model") +
  # ylim(-1,1)+
  scale_color_manual(values= color_labels)+
  scale_linetype_manual(values = lines) +
  facet_wrap(region ~.)
a
pdf(paste0("/home/sbelman/Documents/env_sa_manuscript/models/pcvartprovrep_models/yearlyrandom_provRep_perturbations_",endyear,".pdf"), width = 8, height =5)
print(a)
dev.off()

################ FIXED EFFECTS LINEAR - PROVINCE REPLICATION ##############
# fixed_effects <- rbind(modART_pr$summary.fixed[2,],
#                        modPCV_pr$summary.fixed[2,],
#                        mod_ARTPCV_pr$summary.fixed[2:3,])
# fixedeffs <- cbind(fixed_effects, model = c("ART","PCV","ART_PCV","ART_PCV"), effect = c("ART","PCV","PCV","ART"))
# fixedeffs$model <- factor(fixedeffs$model, levels = c("ART_PCV","ART","PCV"))
# 
# colnames(fixedeffs) <- c("mean","sd","lowerCI","median","upperCI","mode","kld","model","effect")
# color_labels <- c("ART" = "#097969","PCV" = "#D91656", "ART_PCV" = "darkblue")
# 
# fixed_prov <- ggplot(fixedeffs)+
#   geom_point(aes(x=effect,y=median, group=model, color= model ),size=3, position=position_dodge(width=0.5)) +
#   geom_errorbar(aes(x=effect,ymin=lowerCI,ymax=upperCI, group=model, color= model), width=.3, position=position_dodge(width=0.5))+
#   theme_bw()+
#   xlab("Covariate")+
#   ylab("Fixed Effect")+
#   geom_hline(yintercept=0,linetype="dashed",color="red")+
#   scale_color_manual(values=color_labels)+
#   ggtitle("Fixed Effect\n(interannual effect province replicate)")+
#   theme(axis.text.x = element_text(angle=45,hjust=1), axis.text = element_text(size=13),axis.title=element_text(size=13)) +
#   labs(color="Model",
#        caption = "*PCV coverage data is only at the national level")
# fixed_national + fixed_prov
