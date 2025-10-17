

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
            endyear = 2019
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
  xlab("Intervention")+
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
  xlab("Intervention")+
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
# -----------------------------
# ️ Compute raw correlations and classical VIF
# -----------------------------
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

################# INTERANNUAL RANDOM EFFECTS ###################################
### interannual random effects from the same models
# plot_random_effects(mod_base, modART, title_a = "base", title_b = "ART", return = "year")
# plot_random_effects(mod_base, mod_ARTnational, title_a = "base", title_b = "ART national", return = "year")
# plot_random_effects(mod_base, modPCV, title_a = "base", title_b = "PCV", return = "year")
# plot_random_effects(mod_base, mod_vaxperiod, title_a = "base", title_b = "vaccine_periods", return = "year")
# plot_random_effects(mod_base, mod_2009, title_a = "base", title_b = "2009vaccine", return = "year")
# plot_random_effects(mod_base, mod_vaxART, title_a = "base", title_b = "ART and VaccinePeriods", return = "year")
# 
# ### interannual random effects from province replicate models
# plot_random_effect_provRep(mod_base_pr, modART_pr, title_a = "base", title_b = "ART", type = "annual")
# plot_random_effect_provRep(mod_base_pr, mod_ARTnational_pr, title_a = "base", title_b = "ART", type = "annual")
# plot_random_effect_provRep(mod_base_pr, modPCV_pr, title_a = "base", title_b = "PCV", type = "annual")
# plot_random_effect_provRep(mod_base_pr, mod_vaxperiod_pr, title_a = "base", title_b = "vaccine_periods", type = "annual")
# plot_random_effect_provRep(mod_base_pr, mod_2009_pr, title_a = "base", title_b = "2009vaccine", type = "annual")
# plot_random_effect_provRep(mod_base_pr, mod_vaxART_pr, title_a = "base", title_b = "ART and VaccinePeriods", type = "annual")

### plot for supplement
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
# 
# ################ NONLINEAR EFFECTS - NO PROVINCE REPLICATION ##############
# color_labels <- c("ART" = "#097969","PCV" = "#D91656", "ART_PCV" = "darkblue")
# nl_effects <- rbind(mod_ART_NL$summary.random[[4]],
#                        mod_PCV_NL$summary.random[[4]],
#                        mod_PCVART_NL$summary.random[[4]],
#                     mod_PCVART_NL$summary.random[[5]])
# colnames(nl_effects) <- c("coverage","mean","sd","lowerCI","median","upperCI","mode","kld")
# nl_effects <- nl_effects %>%
#   mutate(
#     RR_mean = exp(mean),
#     RR_lower = exp(lowerCI),
#     RR_upper = exp(upperCI)
#   )
# ## replicating by the number of cuts which was 5 but the PCV ones only ran 3 cuts due to the early zeroes so the model with both has 8 rows
# nl_effects <- cbind(nl_effects, model = c(rep("ART",5),rep("PCV",4),rep("ART_PCV",9)), effect = c(rep("ART",5),rep("PCV",4),rep("ART",5), rep("PCV",4)))
# nl_effects$model <- factor(nl_effects$model, levels = c("ART_PCV","ART","PCV"))
# nl_all <- ggplot(nl_effects) +
#   # geom_line(aes(x=coverage, y=median, group=model, group=effect, color=model))+
#   # geom_ribbon(aes(x=coverage,ymin=lowerCI,ymax=upperCI, fill=model,group=model, group=effect),alpha=0.2)+
#   geom_pointrange(aes(x=coverage, y=RR_mean, ymin=RR_lower,ymax=RR_upper, group=model, group=effect, color=model), position=position_dodge(width=0.1))+
#   theme_bw()+
#   xlab("Coverage")+
#   ylab("Effect")+
#   geom_hline(yintercept=1,linetype="dashed",color="red")+
#   theme(axis.text.x = element_text(angle=45,hjust=1), axis.text = element_text(size=13),axis.title=element_text(size=13)) +
#   scale_color_manual(values=color_labels)+
#   scale_fill_manual(values=color_labels)+
#   labs(fill="Model", color = "Model")+
#   facet_grid(effect~.)
#   
# 
# 
# ### just above and below 70% not all  for no province replication
# # Categorize coverage as <70% or >70%
# nl_effects <- nl_effects %>%
#   mutate(coverage_cat = ifelse(coverage < 0.7, "<70%", ">70%"))
# 
# # Summarize data: mean, lowerCI and upperCI per group and coverage category
# summary_df <- nl_effects %>%
#   group_by(model, effect, coverage_cat) %>%
#   summarise(
#     mean_effect = mean(RR_mean),
#     lower_CI = mean(RR_lower),
#     upper_CI = mean(RR_upper),
#     .groups = "drop"
#   )
# 
# # Create a new column for x-axis labels (coverage + model)
# summary_df <- summary_df %>%
#   mutate(x_label = paste(coverage_cat, effect, sep = "_"))
# write.csv(summary_df, file = "/home/sbelman/Documents/BRD/SouthAfrica/manuscript/figure2/coverage_cat_summarydf.csv")
# 
# # Define color palette
# color_labels <- c("ART" = "#097969", "PCV" = "#D91656", "ART_PCV" = "darkblue")
# 
# # Define dodge position object
# pd <- position_dodge(width = 0.6)
# 
# ggplot(summary_df, aes(x = effect, y = mean_effect, color = model)) +
#   geom_point(position = pd, size = 3) +
#   geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI),
#                 width = 0.2,
#                 position = pd) +
#   facet_wrap(~coverage_cat) +
#   scale_color_manual(values = color_labels) +
#   labs(
#     x = "Intervention",
#     y = "Relative Risk",
#     color = "Model"
#   ) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
#   theme_bw(base_size = 15) +
#   scale_y_continuous(trans = "log10" )+
#   # ylim(-0.7, 0.7) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "bottom"
#   )
# 
# 
# ################ NONLINEAR EFFECTS -  PROVINCE REPLICATION ##############
# color_labels <- c("ART" = "#097969","PCV" = "#D91656", "ART_PCV" = "darkblue")
# nl_effects <- rbind(mod_ART_prNL$summary.random[[4]],
#                     mod_PCV_prNL$summary.random[[4]],
#                     mod_PCVART_prNL$summary.random[[4]],
#                     mod_PCVART_prNL$summary.random[[5]])
# colnames(nl_effects) <- c("coverage","mean","sd","lowerCI","median","upperCI","mode","kld")
# ## replicating by the number of cuts which was 5 but the PCV ones only ran 35cuts due to the early zeroes so the model with both has 8 rows
# nl_effects <- cbind(nl_effects, model = c(rep("ART",5),rep("PCV",4),rep("ART_PCV",9)), effect = c(rep("ART",5),rep("PCV",4),rep("PCV",4), rep("ART",5)))
# nl_effects$model <- factor(nl_effects$model, levels = c("ART_PCV","ART","PCV"))
# 
# #### pull out the intercepts so that I can calculate the incidence
# nl_intercepts <- rbind(mod_ART_prNL$summary.fixed[1,c(1,3:5)],
#                     mod_PCV_prNL$summary.fixed[1,c(1,3:5)],
#                     mod_PCVART_prNL$summary.fixed[1,c(1,3:5)],
#                     mod_PCVART_prNL$summary.fixed[1,c(1,3:5)])
# colnames(nl_intercepts) <- c("intercept_mean", "intercept_lowerCI","intercept_median","intercept_upperCI")
# nl_intercepts$model <- c("ART","PCV","ART_PCV","ART_PCV")
# ### join the intercepts
# # nl_effects <- left_join(nl_effects,nl_intercepts,by="model")
# # ## calculate the cases at each coverage
# # nl_effects$exp_mod <- exp(nl_effects$mean)
# # nl_effects$est_weekly_cases <-  nl_effects$exp_mod *500
# # nl_effects$exp_int <- exp(nl_effects$intercept_mean) 
# # nl_effects$est_intercept_cases <-  nl_effects$exp_int *500
# # nl_effects$RRs <- nl_effects$est_weekly_cases/nl_effects$est_intercept_cases
# 
# nl_prov <- ggplot(nl_effects) +
#   # geom_line(aes(x=coverage, y=median, group=model, group=effect, color=model))+
#   # geom_ribbon(aes(x=coverage,ymin=lowerCI,ymax=upperCI, fill=model,group=model, group=effect),alpha=0.2)+
#   geom_pointrange(aes(x=coverage, y=median, ymin=lowerCI,ymax=upperCI, group=model, group=effect, color=model), position=position_dodge(width=0.1))+
#   
#   theme_bw()+
#   xlab("Coverage")+
#   ylab("Effect")+
#   geom_hline(yintercept=0,linetype="dashed",color="red")+
#   theme(axis.text.x = element_text(angle=45,hjust=1), axis.text = element_text(size=13),axis.title=element_text(size=13)) +
#   scale_color_manual(values=color_labels)+
#   scale_fill_manual(values=color_labels)+
#   labs(fill="Model", color = "Model")+
#   facet_grid(effect~.)
# nl_all + nl_prov
# 
# 
# ################################################################################
# #### EXTRACT INTERCEPT FOR EACH MODEL
# ################################################################################
# 
# ### intercept and estimate the actual impacts
# intercept_effs <- rbind(mod_base$summary.fixed[1,], modART$summary.fixed[1,],
#                         modPCV$summary.fixed[1,],
#                         mod_ARTPCV$summary.fixed[1,])
# intercept_effs <- cbind(intercept_effs, model = c("BASE","ART","PCV","ART_PCV"))
# intercept_effs$model <- factor(intercept_effs$model, levels = c("BASE","ART_PCV","ART","PCV"))
# intercept_effs$exp_mod <- exp(intercept_effs$mean)
# intercept_effs$est_weekly_cases <-  intercept_effs$exp_mod *500
# 
# ### intercept and estimate the actual impacts
# intercept_effs_pr <- rbind(mod_base_pr$summary.fixed[1,], modART_pr$summary.fixed[1,],
#                         modPCV_pr$summary.fixed[1,],
#                         mod_ARTPCV_pr$summary.fixed[1,])
# intercept_effs_pr <- cbind(intercept_effs_pr, model = c("BASE","ART","PCV","ART_PCV"))
# intercept_effs_pr$model <- factor(intercept_effs_pr$model, levels = c("BASE","ART_PCV","ART","PCV"))
# intercept_effs_pr$exp_mod <- exp(intercept_effs_pr$mean)
# intercept_effs_pr$est_weekly_cases <-  intercept_effs_pr$exp_mod *500
# 
# ### intercept from NL models
# intercept_effs_prNL <- rbind(mod_base_pr$summary.fixed[1,], mod_ART_NL$summary.fixed[1,],
#                              mod_PCV_NL$summary.fixed[1,],
#                              mod_PCVART_NL$summary.fixed[1,])
# intercept_effs_prNL <- cbind(intercept_effs_prNL, model = c("BASE","ART","PCV","ART_PCV"))
# intercept_effs_prNL$model <- factor(intercept_effs_prNL$model, levels = c("BASE","ART_PCV","ART","PCV"))
# intercept_effs_prNL$exp_mod <- exp(intercept_effs_prNL$mean)
# intercept_effs_prNL$est_weekly_cases <-  intercept_effs_prNL$exp_mod *500
# 
# ################################################################################
# #### EXTRACT INTERANNUAL RANDOM EFFECTS
# ################################################################################
# ####### NO PROVICNE REPLICATION ###########
# ####### plot all the yearly effects together to compare ###########
# first_year <- 2005
# # year_eff_base <- rbind(modART$summary.random$id_y, modPCV$summary.random$id_y, modPCV_timing$summary.random$id_y, mod_ARTPCV$summary.random$id_y, mod_base$summary.random$id_y)
# year_eff_base <- rbind(modART$summary.random$id_y, modPCV$summary.random$id_y,  mod_ARTPCV$summary.random$id_y, mod_base$summary.random$id_y)
# N_years <- nrow(modART$summary.random$id_y)
# colnames(year_eff_base)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld")
# year_eff_base <- year_eff_base %>%
#   mutate(
#     RR_mean = exp(mean),
#     RR_lower = exp(lowerCI),
#     RR_upper = exp(upperCI)
#   )
# 
# mod_names <- rep(c("ART_Coverage","PCV_Coverage","PCV_ART","Base"), each = N_years)
# year_eff_base$model <- mod_names
# year_eff_base$model <- factor(year_eff_base$model, levels = c("Base","ART_Coverage","PCV_Coverage", "PCV_ART"))
# year_eff_base$ID <- first_year + seq(0,(N_years-1),)
# 
# lines <- c("ART_Coverage" = "longdash","PCV_Coverage" = "longdash", "PCV_ART" = "dotted", "Base" = "solid")
# color_labels <- c("ART_Coverage" = "#097969","PCV_Coverage" = "#D91656", "PCV_ART" = "darkblue", "Base"= "black")
# random_national <- ggplot(year_eff_base)+
#   geom_hline(yintercept=1,linetype="dashed",color="grey",alpha=0.8)+
#   # geom_line(aes(x=ID,y=median,linetype=model, color = model, group=model))+
#   # geom_ribbon(aes(x=ID,ymin=lowerCI,ymax=upperCI, fill=model,group=model),alpha=0.2)+
#   geom_line(aes(x=ID,y=RR_mean,linetype=model, color = model, group=model))+
#   geom_ribbon(aes(x=ID,ymin=RR_lower,ymax=RR_upper, fill=model,group=model),alpha=0.2)+
#   scale_linetype_manual(values = lines)+
#   scale_color_manual(values=color_labels)+
#   scale_fill_manual(values = color_labels)+
#   theme_classic()+
#   theme(axis.text = element_text(size=13),axis.title = element_text(size = 13), axis.text.x =element_text(angle = 45, hjust =1)
#   )+
#   labs(x = "Year",
#        y = "Relative Risk") +
#   scale_y_continuous(trans = "log10")+
#   # ylim(-2,2)+
#   ggtitle("Interannual Effect")
# # 
# # png(file = "/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/base_models/effect_artpcv.png", width=400, height=700)
# # fixed_national/random_national
# # dev.off()
# 
# write.csv(year_eff_base, file = "/home/sbelman/Documents/BRD/SouthAfrica/manuscript/figure2/year_eff_base.csv")
# 
# ########PROVINCE REPLICATION ###########
# ####### plot all the yearly effects together to compare ##########
# first_year <- 2005
# year_eff_base <- rbind(modART_pr$summary.random$id_y, modPCV_pr$summary.random$id_y,  mod_ARTPCV_pr$summary.random$id_y, mod_base_pr$summary.random$id_y)
# N_years <- nrow(modART$summary.random$id_y)
# reg_year <- rep(c("Eastern Cape", "Free State", "Gauteng", 
#                   "KwaZulu-Natal", "Limpopo", "Mpumalanga", 
#                   "North West", "Northern Cape", "Western Cape"), each = N_years)
# year_eff_base$region<-rep(reg_year, 4)
# 
# colnames(year_eff_base)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld","region")
# mod_names <- rep(c("ART_Coverage","PCV_Coverage","PCV_ART","Base"), each = N_years*9)
# year_eff_base$model <- mod_names
# year_eff_base$model <- factor(year_eff_base$model, levels = c("Base","ART_Coverage","PCV_Coverage", "PCV_ART"))
# year_eff_base$ID <- first_year + seq(0,(N_years-1),)
# 
# lines <- c("ART_Coverage" = "dotted","PCV_Coverage" = "longdash", "PCV_ART" = "longdash", "Base" = "solid")
# color_labels <- c("ART_Coverage" = "#410445","PCV_Coverage" = "#D91656", "PCV_ART" = "darkblue", "Base"= "black")
# random_prov <- ggplot(year_eff_base)+
#   geom_hline(yintercept=0,linetype="dashed",color="grey",alpha=0.8)+
#   geom_line(aes(x=ID,y=median,linetype=model, color = model, group=model, group=region))+
#   # geom_ribbon(aes(x=ID,ymin=lowerCI,ymax=upperCI, fill=model,group=model),alpha=0.2)+
#   scale_linetype_manual(values = lines)+
#   scale_color_manual(values=color_labels)+
#   scale_fill_manual(values = color_labels)+
#   theme_classic()+
#   theme(axis.text = element_text(size=13), axis.text.x = element_text(angle=45,hjust=1),axis.title = element_text(size = 13))+
#   labs(x = "Year",
#        y = "Median Effect") +
#   # ylim(-2,2)+
#   ggtitle("Interannual Effect Replicated by Province") +
#   facet_wrap(.~region)
# 
# random_national + random_prov
# 
# 
# ################### PROVINCE WITH NONLINEAR MODEL ########################
# first_year <- 2005
# year_eff_base <- rbind(mod_ART_prNL$summary.random$id_y, mod_PCV_prNL$summary.random$id_y,  mod_PCVART_prNL$summary.random$id_y, mod_base_pr$summary.random$id_y)
# N_years <- nrow(modART$summary.random$id_y)
# reg_year <- rep(c("Eastern Cape", "Free State", "Gauteng", 
#                   "KwaZulu-Natal", "Limpopo", "Mpumalanga", 
#                   "North West", "Northern Cape", "Western Cape"), each = N_years)
# year_eff_base$region<-rep(reg_year, 4)
# 
# colnames(year_eff_base)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld","region")
# mod_names <- rep(c("ART_Coverage","PCV_Coverage","PCV_ART","Base"), each = N_years*9)
# year_eff_base$model <- mod_names
# year_eff_base$model <- factor(year_eff_base$model, levels = c("Base","ART_Coverage","PCV_Coverage", "PCV_ART"))
# year_eff_base$ID <- first_year + seq(0,(N_years-1),)
# 
# lines <- c("ART_Coverage" = "dotted","PCV_Coverage" = "longdash", "PCV_ART" = "longdash", "Base" = "solid")
# color_labels <- c("ART_Coverage" = "#410445","PCV_Coverage" = "#D91656", "PCV_ART" = "darkblue", "Base"= "black")
# random_prov_NL <- ggplot(year_eff_base)+
#   geom_hline(yintercept=0,linetype="dashed",color="grey",alpha=0.8)+
#   geom_line(aes(x=ID,y=median,linetype=model, color = model, group=model, group=region))+
#   # geom_ribbon(aes(x=ID,ymin=lowerCI,ymax=upperCI, fill=model,group=model),alpha=0.2)+
#   scale_linetype_manual(values = lines)+
#   scale_color_manual(values=color_labels)+
#   scale_fill_manual(values = color_labels)+
#   theme_classic()+
#   theme(axis.text = element_text(size=13), axis.text.x = element_text(angle=45,hjust=1),axis.title = element_text(size = 13))+
#   labs(x = "Year",
#        y = "Median Effect") +
#   # ylim(-2,2)+
#   ggtitle("Interannual Effect Replicated by Province (Nonlinear)") +
#   facet_wrap(.~region)
# 
# random_prov_NL + random_prov




######################### NON LINEAR MODELS FROM BEFORE ########################
# 
# ##### non-linear models #####
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'f(inla.group(PCV_coverage, n = 5 , method = "quantile"), model = "rw1")'),
#                     "disease")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("spatial, seasonal, and annual accounting for PCV non-linear")
# base_pcvart_formula_list[[9]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'f(inla.group(ART_coverage, n = 5 , method = "quantile"), model = "rw1")'),
#                     "disease")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("spatial, seasonal, and annual accounting for ART non-linear")
# base_pcvart_formula_list[[10]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'f(inla.group(ART_coverage, n = 5 , method = "quantile"), model = "rw1")',
#                       'f(inla.group(PCV_coverage, n = 5 , method = "quantile"), model = "rw1")'),
#                     "disease")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("spatial, seasonal, and annual accounting for PCV and ART non-linear")
# base_pcvart_formula_list[[11]]<-base_form
# 
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid",  replicate = id_prov)',
#                       'f(inla.group(PCV_coverage, n = 5 , method = "quantile"), model = "rw1")'),
#                     "disease")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for PCV non-linear")
# base_pcvart_formula_list[[12]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid",  replicate = id_prov)',
#                       'f(inla.group(ART_coverage, n = 5 , method = "quantile"), model = "rw1")'),
#                     "disease")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for ART non-linear")
# base_pcvart_formula_list[[13]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid", replicate = id_prov)',
#                       'f(inla.group(PCV_coverage, n = 5 , method = "quantile"), model = "rw1")',
#                       'f(inla.group(ART_coverage, n = 5 , method = "quantile"), model = "rw1")'),
#                     "disease")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for PCV and ART non-linear")
# base_pcvart_formula_list[[14]]<-base_form
# 
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid", replicate = id_prov)',
#                       'ART_coverage_national'),
#                     "disease")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for ART (national coverage)")
# base_pcvart_formula_list[[15]]<-base_form



# sapply(base_pcvart_formula_list, function(x) x$formula)









########### BEFORE MAKING THE FORMULAT THAT IS IN SOURCE FUNCTIONS
# 
# base_pcvart_nvt_formula_list <- list()
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'ART_coverage'),
#                     "nvt")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("NVT spatial, seasonal, and annual accounting for ART coverage")
# base_pcvart_nvt_formula_list[[1]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'PCV_coverage'),
#                     "nvt")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("NVT spatial, seasonal, and annual accounting for PCV coverage")
# base_pcvart_nvt_formula_list[[2]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'vaccination_period'),
#                     "nvt")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("NVT spatial, seasonal, and annual accounting for vaccination_period")
# base_pcvart_nvt_formula_list[[3]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'post_vaccination_2009'),
#                     "nvt")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("NVT spatial, seasonal, and annual accounting for 2009 vaccine")
# base_pcvart_nvt_formula_list[[4]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'ART_coverage_national'),
#                     "nvt")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("NVT spatial, seasonal, and annual accounting for ART Nationally")
# base_pcvart_nvt_formula_list[[5]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'vaccination_period',
#                       'ART_coverage'),
#                     "nvt")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("NVT spatial, seasonal, and annual accounting for vaccination_period and ART coverage")
# base_pcvart_nvt_formula_list[[6]]<-base_form
# 
# 
# ##### include VACCINE TYPES types only and no province replication ##### ##### 
# base_pcvart_vts_formula_list <- list()
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'ART_coverage'),
#                     "vts")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("vts spatial, seasonal, and annual accounting for ART coverage")
# base_pcvart_vts_formula_list[[1]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'PCV_coverage'),
#                     "vts")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("vts spatial, seasonal, and annual accounting for PCV coverage")
# base_pcvart_vts_formula_list[[2]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'vaccination_period'),
#                     "vts")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("vts spatial, seasonal, and annual accounting for vaccination_period")
# base_pcvart_vts_formula_list[[3]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'post_vaccination_2009'),
#                     "vts")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("vts spatial, seasonal, and annual accounting for 2009 vaccine")
# base_pcvart_vts_formula_list[[4]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'ART_coverage_national'),
#                     "vts")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("vts spatial, seasonal, and annual accounting for ART Nationally")
# base_pcvart_vts_formula_list[[5]]<-base_form
# 
# base_form<-list()
# form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
#                       'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                       'f(id_y, model = "iid")',
#                       'vaccination_period',
#                       'ART_coverage'),
#                     "vts")
# base_form$formula <- as.formula(form)
# base_form$covs<-paste0("vts spatial, seasonal, and annual accounting for vaccination_period and ART coverage")
# base_pcvart_vts_formula_list[[6]]<-base_form

