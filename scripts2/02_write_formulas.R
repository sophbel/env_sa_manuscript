################################################################################
#### LOAD DATA AND LIBRARIES
################################################################################
source("/home/sbelman/Documents/BRD/scripts/0_source_functions.R")
data<-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_weekly_lag_sc.csv")
# data<-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_monthly_lag.csv")

df<- data
df <- df %>%
mutate(date = as.Date(date))
## create interaction terms
df <- df %>%
  mutate(
    post_vaccination_2009 = ifelse(date >= as.Date("2009-01-01"), 1, 0),
    post_vaccination_2011 = ifelse(date >= as.Date("2011-01-01"), 1, 0),
    post_covid_2020 = ifelse(date >= as.Date("2020-01-01"), 1, 0),
    time_since_vaccination_2009 = as.numeric(date - as.Date("2009-01-01")) * post_vaccination_2009,
    time_since_vaccination_2011 = as.numeric(date - as.Date("2011-01-01")) * post_vaccination_2011,
    time_since_covid_2020 = as.numeric(date - as.Date("2020-01-01")) * post_covid_2020
  )

## vectors with variable names
all_names<-colnames(df)
### environmental
env_vars<-c(grep("sp|tas|hurs|absh|prlr|pm|o3|so|Wind",all_names,value=TRUE))
gene_vars<-grep("monthly_count",all_names,value=TRUE)

### include province as factors for replications
df$id_prov <- as.numeric(factor(df$NAME_1))

################################################################################
#### BASE MODEL TESTING
################################################################################
base_form_list<-list()
base_form<-list()
form <- as.formula(disease ~ 1)
base_form$formula <- as.formula(form)
base_form$covs<-paste0("Intercept")
base_form_list[[1]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial and seasonal")
base_form_list[[2]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("seasonal, and annual")
base_form_list[[3]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, and annual")
base_form_list[[4]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual")
base_form_list[[5]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (rep by prov)")
base_form_list[[6]]<-base_form


base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T,  adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'post_vaccination_2009'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for 2009 vaccine")
base_form_list[[7]]<-base_form


base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T,  adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'vaccination_period'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for both vaccines")
base_form_list[[8]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid",hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'vaccination_period',
                      'post_covid_2020' ),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for both vaccines and covid")
base_form_list[[9]]<-base_form




############ include replications ######
base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'vaccination_period'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal (rep by prov), and annual accounting for both vaccines")
base_form_list[[10]]<-base_form


base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'vaccination_period'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (rep by prov) accounting for both vaccines")
base_form_list[[11]]<-base_form


base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'post_vaccination_2009'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (rep by prov) accounting for 2009 vaccine")
base_form_list[[12]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'vaccination_period'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal (rep by prov), and annual (rep by prov) accounting for both vaccines")
base_form_list[[13]]<-base_form


base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'vaccination_period',
                      'post_covid_2020'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (rep by prov) accounting for both vaccines and covid")
base_form_list[[14]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'vaccination_period',
                      'population_density'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (rep by prov) accounting for both vaccines and population_density")
base_form_list[[15]]<-base_form


base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'vaccination_period',
                      'post_covid_2020'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial (rep by prov), seasonal (rep by prov), and annual (rep by prov) accounting for both vaccines and covid")
base_form_list[[16]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T,hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'PCV_coverage',
                      'population_density'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (rep by prov) accounting for PCV Coverage and population_density")
base_form_list[[17]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'ART_coverage',
                      'population_density'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (rep by prov) accounting for ART Coverage and population_density")
base_form_list[[18]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'ART_coverage',
                      'vaccination_period'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (rep by prov) accounting for ART Coverage and vaccination")
base_form_list[[19]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'ART_coverage',
                      'PCV_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (rep by prov) accounting for ART Coverage and PCV coverage")
base_form_list[[20]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'f(id_y, model = "iid", replicate = id_prov, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01))))',
                      'ART_coverage',
                      'vaccination_period',
                      'population_density'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (rep by prov) accounting for ART Coverage and vaccination period and population density")
base_form_list[[21]]<-base_form


saveRDS(base_form_list,file="/home/sbelman/Documents/env_sa_manuscript/formulas/base_form_list.rds")


################################################################################
#### SOCIODEMOGRAPHIC UNIVARIABLE MODELS
################################################################################
data<-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_weekly_lag.csv")
all_names <- colnames(data)
sociodem_vars <- grep("lag|count|GID_2|adm2_name",all_names, invert = TRUE, value = TRUE)
sociodem_vars <- sociodem_vars[c(12,19:35)]
sociodem_formulas <-sapply(sociodem_vars,function(x) construct_formula(x, "mod7"))
saveRDS(sociodem_formulas,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/formulas/sociodem_formulas.rds")


################################################################################
#### LINEARITY TESTING
################################################################################
cov_names<- c(grep("lag0",env_vars,value=TRUE,invert=FALSE))
lin_formulas <- list()

for(i in 1:length(cov_names)) {
  lin_list<-list()
  form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T,  adjust.for.con.comp=TRUE)',
                        'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                        'f(id_y, model = "iid")',
                        'vaccination_period',
                        'population_density',
                        paste(cov_names[i], collapse="+")),
                      "disease")
  lin_list$formula <-as.formula(form)
  lin_list$covs <- cov_names[i]
  lin_formulas[[i]]<-lin_list
}


# explore non-linear effects with 8 cuts----------------------------------------
nonlin_formulas<-list()

for(i in 1:length(cov_names)) {
  nonlin_list<-list()
  form <- reformulate(c(1,'f(id_u, model = "bym2", graph = g, scale.model = T,  adjust.for.con.comp=TRUE)',
                         'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                         'f(id_y, model = "iid", replicate = id_prov)',
                         'vaccination_period',
                        'population_density',
                        paste0("f(inla.group(", cov_names[i], ", n = 8", ",method='cut'), model='rw2')")), 
                      "disease")

  nonlin_list$formula <- as.formula(form)
  nonlin_list$covs <- cov_names[i]
  nonlin_formulas[[i]] <- nonlin_list
}

saveRDS(nonlin_formulas,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/formulas/nonlin_formulas.rds")
saveRDS(lin_formulas,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/formulas/lin_formulas.rds")

################################################################################
#### UNIVARIABLE MODELS
################################################################################
# #######################  MONTHLY
# data<-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_weekly_lag_sc.csv")
# 
# all_names<-colnames(data)
# cov_names<-grep("lag",all_names,value=TRUE)
# ## extract those which are nonlinear
# nls <- c("tas", "tasmin", "tasmax")
# # Create a pattern that matches any of the nonlinear elements
# pattern <- paste(paste0(nls,"_"), collapse = "|")
# # Filter cov_names based on the pattern
# cov_names[grepl(pattern,cov_names)] <-paste0(cov_names[grepl(pattern,cov_names)],"_nl")
# 
# ## write formulas with 0-12 lags
# univ_weekly_formulas<-sapply(cov_names,function(x) construct_formula(x, "mod7"))
# saveRDS(univ_weekly_formulas,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/formulas/univ_weekly_formulas_lags12.rds")
# 
# 
# ## write formulas with 6 lags
# cov_names_6 <- grep(paste0("lag",seq(7,12,1), collapse="|"),cov_names, value=TRUE, invert = TRUE)
# univ_weekly_formulas6<-sapply(cov_names_6,function(x) construct_formula(x, "mod7"))
# saveRDS(univ_weekly_formulas6,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/formulas/univ_weekly_formulas_lags6.rds")
# 
# ## write formulas with 7-12 lags
# cov_names_712 <- grep(paste0("lag",seq(7,12,1), collapse="|"),cov_names, value=TRUE, invert = FALSE)
# univ_weekly_formulas712<-sapply(cov_names_712,function(x) construct_formula(x, "mod7"))
# saveRDS(univ_weekly_formulas712,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/formulas/univ_weekly_formulas_lags712.rds")
# ################################################################################
# #### BIVARIABLE MODELS WITH 0 AND 2 MONTH LAGS ACROSS VARIABLES
# ################################################################################
# ## pull out all 0 and 2 month lags
# data<-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_monthly_lag_sc.csv")
# all_names<-colnames(data)
# cov_names<-grep(c("lag0|lag2"),all_names,value=TRUE)
# cov_names <- grep(c("o3|spi|spei|tasmin"),cov_names, invert=TRUE, value= TRUE)
# tas_names <- grep("tas", cov_names, value=TRUE)
# notas_names <- grep("tas", cov_names, value=TRUE, invert=TRUE)
# bivs <- expand.grid(tas_names,notas_names)
# 
# ## create formula 
# biv_form_list<-list()
# for(i in 1:nrow(bivs)){
#   base_form<-list()
#   form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T,  adjust.for.con.comp=TRUE)',
#                         'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                         'f(id_y, model = "iid", replicate = id_prov)',
#                         'post_vaccination_2009 * time_since_vaccination_2009',
#                         'post_vaccination_2011 * time_since_vaccination_2011',
#                         'post_covid_2020',
#                         paste(bivs[i,1],"+",bivs[i,2])),
#                     'disease')
#   base_form$formula <- as.formula(form)
#   base_form$covs<-paste0(bivs[i,1],"_",bivs[i,2])
#   biv_form_list[[i]]<-base_form
# }
# 
# saveRDS(biv_form_list,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/formulas/bivariable_formula_list_lags.rds")
# 
# 
# ################################################################################
# #### MULTIVARIATE FORMULATION FOR DIFFERENT DISEASE TYPES
# ################################################################################
# data<-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_monthly_lag_sc.csv")
# all_names<-colnames(data)
# all_names <- colnames(data)
# cov_names <- grep(c("lag0"), all_names, value = TRUE)
# cov_names <- grep(c("spi|spei|absh"), cov_names, invert = TRUE, value = TRUE)
# 
# ########### FOR DIFFERENT DISEASE TYPES ##############
# ### create multivariate formulas 
# multivar_form_list <- list()
# # cov_names <- cov_names[1:2]
# for (covariate in 1:length(cov_names)) {
#   base_form <- list()
#   
#   # form_str <- paste0("count ~ -1 +  outcome:",covariate," +",
#   form_str <- paste0("count ~ -1 +  outcome:cov +",
#     'f(id_u, model = "bym2", graph = g, replicate = outcome_idx) + ',  # Spatial effect per outcome
#     'f(id_m, model = "rw2", cyclic = TRUE, scale.model = TRUE, constr = TRUE, replicate = outcome_idx) + ',# seasonal effect per outcome
#     # 'f(id_y, model = "iid", replicate = outcome_idx)  +', #interannual effect per outcome
#     'f(id_y, model = "iid", replicate = outcome_idx_prov)  +', #interannual effect per outcome and province
#     'post_vaccination_2009 +', 
#     'post_vaccination_2011 +',
#     'population_density' )
#    
#   base_form$formula <- as.formula(form_str)
#   base_form$cov <- cov_names[covariate]
# 
#   multivar_form_list[[covariate]] <- base_form
# }
# 
# saveRDS(multivar_form_list,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/formulas/multivariate_univariable_formula_list.rds")
# 
# 
# 
# ########### FOR DIFFERENT EACH GPSC ##############
# ### create multivariate formulas 
# multivarGPSC_form_list <- list()
# # cov_names <- cov_names[1:2]
# for (covariate in 1:length(cov_names)) {
#   base_form <- list()
#   
#   # form_str <- paste0("count ~ -1 +  outcome:",covariate," +",
#   form_str <- paste0("count ~ -1 +  outcome:cov +",
#                      'f(id_u, model = "bym2", graph = g, replicate = outcome_idx, hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01)))) + ',  # Spatial effect per outcome
#                      'f(id_m, model = "rw2", cyclic = TRUE, scale.model = TRUE, constr = TRUE, replicate = outcome_idx, hyper = list(prec = list(prior = "pc.prec", param = c(4, 0.01)))) + ',# seasonal effect per outcome
#                      'f(id_y, model = "iid", replicate = outcome_idx) +', #interannual effect per outcome
#                      'post_vaccination_2009 +', 'post_vaccination_2011 +', 'population_density')
#   # 'post_covid_2020 +',
#   # 'offset(log(population_size/100000))')
#   
#   base_form$formula <- as.formula(form_str)
#   base_form$cov <- cov_names[covariate]
#   
#   multivarGPSC_form_list[[covariate]] <- base_form
# }
# 
# saveRDS(multivarGPSC_form_list,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/formulas/multivariateGPSC_univariable_formula_list.rds")
# 
# 
# ################################################################################
# #### GENOMIC OUTCOME MODELS
# ################################################################################
# data_monthly<-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_monthly_lag.csv")
# data_monthly$date <- as.Date(data_monthly$date)
# 
# df <- data_monthly %>%
#   mutate(
#     post_vaccination_2009 = ifelse(date >= as.Date("2009-01-01"), 1, 0),
#     post_vaccination_2011 = ifelse(date >= as.Date("2011-01-01"), 1, 0),
#     post_covid_2020 = ifelse(date >= as.Date("2020-01-01"), 1, 0),
#     time_since_vaccination_2009 = as.numeric(date - as.Date("2009-01-01")) * post_vaccination_2009,
#     time_since_vaccination_2011 = as.numeric(date - as.Date("2011-01-01")) * post_vaccination_2011,
#     time_since_covid_2020 = as.numeric(date - as.Date("2020-01-01")) * post_covid_2020
#   )
# 
# gpscs<-grep("count",grep("GPSC",names(df),value = TRUE),value=TRUE)
# gpsc_form_list<-list()
# for(i in 1:length(gpscs)){
#   base_form<-list()
#     form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T,  adjust.for.con.comp=TRUE)',
#                           'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                           'f(id_y, model = "iid")',
#                           'post_vaccination_2009 * time_since_vaccination_2009',
#                           'post_vaccination_2011 * time_since_vaccination_2011'),
#                         gpscs[i])
#     base_form$formula <- as.formula(form)
#     base_form$covs<-paste0(gpscs[i])
#     gpsc_form_list[[i]]<-base_form
# }
# saveRDS(gpsc_form_list,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/formulas/gpsc_outcome_formula_list.rds")
# 
# 
# #### NO SEASONAL
# gpscs<-grep("count",grep("GPSC",names(df),value = TRUE),value=TRUE)
# gpsc_form_list_novax<-list()
# for(i in 1:length(gpscs)){
#   base_form<-list()
#   form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T,  adjust.for.con.comp=TRUE)',
#                         'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
#                         'f(id_y, model = "iid")'),
#                         # 'post_vaccination_2009 * time_since_vaccination_2009',
#                         # 'post_vaccination_2011 * time_since_vaccination_2011'),
#                       gpscs[i])
#   base_form$formula <- as.formula(form)
#   base_form$covs<-paste0(gpscs[i])
#   gpsc_form_list_novax[[i]]<-base_form
# }
# saveRDS(gpsc_form_list_novax,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/formulas/gpsc_outcome_formulaNoVax_list.rds")

