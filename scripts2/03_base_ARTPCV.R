

################################################################################
#### LOAD DATA AND LIBRARIES
################################################################################
source("/home/sbelman/Documents/BRD/scripts/0_source_functions.R")
# data<-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_weekly_lag_sc.csv")
data<-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm2_weekly_lag.csv")
shp<-st_read("/home/sbelman/Documents/BRD/SouthAfrica/shps/gadm41_namematch_ZAF_2.shp")
## read adjacency matrix
g <- inla.read.graph(filename = "/home/sbelman/Documents/BRD/SouthAfrica/shps/sa_adjacency_map.adj")

################################################################################
#### PREPARE DATA
################################################################################
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

################################################################################
#### write a formula list ####
################################################################################
base_pcvart_formula_list <- list()
######### base models for comparison in this list #######
base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual")
base_pcvart_formula_list[[1]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication)")
base_pcvart_formula_list[[2]]<-base_form

######### include coverage of ART and PCV #########
base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")',
                      'ART_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for ART coverage")
base_pcvart_formula_list[[3]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")',
                      'PCV_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for PCV coverage")
base_pcvart_formula_list[[4]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")',
                      'PCV_coverage', 'ART_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for PCV and ART coverage")
base_pcvart_formula_list[[5]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)',
                      'ART_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for ART coverage")
base_pcvart_formula_list[[6]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)',
                      'PCV_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for PCV coverage")
base_pcvart_formula_list[[7]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)',
                      'PCV_coverage', 'ART_coverage'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for PCV and ART coverage")
base_pcvart_formula_list[[8]]<-base_form

##### non-linear models #####

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")',
                      'f(inla.group(PCV_coverage, n = 5 , method = "quantile"), model = "rw1")'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for PCV non-linear")
base_pcvart_formula_list[[9]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")',
                      'f(inla.group(ART_coverage, n = 5 , method = "quantile"), model = "rw1")'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for ART non-linear")
base_pcvart_formula_list[[10]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid")',
                      'f(inla.group(ART_coverage, n = 5 , method = "quantile"), model = "rw1")',
                      'f(inla.group(PCV_coverage, n = 5 , method = "quantile"), model = "rw1")'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual accounting for PCV and ART non-linear")
base_pcvart_formula_list[[11]]<-base_form


base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid",  replicate = id_prov)',
                      'f(inla.group(PCV_coverage, n = 5 , method = "quantile"), model = "rw1")'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for PCV non-linear")
base_pcvart_formula_list[[12]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid",  replicate = id_prov)',
                      'f(inla.group(ART_coverage, n = 5 , method = "quantile"), model = "rw1")'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for ART non-linear")
base_pcvart_formula_list[[13]]<-base_form

base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)',
                      'f(inla.group(PCV_coverage, n = 5 , method = "quantile"), model = "rw1")',
                      'f(inla.group(ART_coverage, n = 5 , method = "quantile"), model = "rw1")'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for PCV and ART non-linear")
base_pcvart_formula_list[[14]]<-base_form


base_form<-list()
form <- reformulate(c(1, 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE)',
                      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T)',
                      'f(id_y, model = "iid", replicate = id_prov)',
                      'ART_coverage_national'),
                    "disease")
base_form$formula <- as.formula(form)
base_form$covs<-paste0("spatial, seasonal, and annual (province replication) accounting for ART (national coverage)")
base_pcvart_formula_list[[15]]<-base_form

base_form<-list()
form <- as.formula(disease ~ 1)
base_form$formula <- as.formula(form)
base_form$covs<-paste0("Intercept")
base_pcvart_formula_list[[16]]<-base_form

sapply(base_pcvart_formula_list, function(x) x$formula)

################################################################################
#### RUN PCV ART MODELS
################################################################################

base_pcvart_model_list <-list()
threads=4
for(i in 1:length(base_pcvart_formula_list)) {

  print(i)
  print("run nbinomial")
  base_pcvart_model_list[[i]] <- inla.mod(base_pcvart_formula_list[[i]], fam = "nbinomial", df = df, nthreads=threads, config=FALSE)
}      
mod_out <- lapply(base_pcvart_model_list, function(mod) eval.mod(mod, df)) # output model evaluation metrics as a list
mod_out <- do.call(rbind, mod_out) # collapse into df
mod_out <- eval.ranks(mod_out) # evaluate ranks by summing across other metrics and rank scoring models

mod_out$model <- rownames(mod_out)
mod_out$cov <- gsub(" ","_",mod_out$cov)
mod_out$cov <- gsub(",","",mod_out$cov)

#### test r2 for all mdoels
mod_out$rsq <- apply(mod_out, 1, function(x) rsq(mod_out[mod_out$cov == x["cov"], ], base_pcvart_model_list[[16]], num_outcomes = 1))



mod_out$num <- gsub("mae","",mod_out$model)
mod_out$num[which(mod_out$num=="")] <- 0
mod_out$num <- as.numeric(mod_out$num)
mod_out$base <- paste0("base","_none")

mod_out_tmp <- mod_out[,c("num","base","waic","mae","cpo","rsq","cov")]

write.table(mod_out,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/pcvart_models/base_pcvart_summary_statisics_weekly.csv",quote=FALSE,row.names=FALSE,col.names=TRUE)
saveRDS(base_pcvart_model_list,file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/pcvart_models/base_pcvart_model_list_weekly.rds")
################################################################################
#### LOAD MODELS AND ASSIGN THEM TO NAMES
################################################################################
##### load only the models we need
base_pcvart_model_list <- readRDS(file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/pcvart_models/base_pcvart_model_list_weekly.rds")
covs <- sapply(base_pcvart_model_list, function(x) x$cov)

### base models 
mod_base <- base_pcvart_model_list[[1]]
mod_base_pr <- base_pcvart_model_list[[2]]

### fixed effect models
modART <- base_pcvart_model_list[[3]]
modPCV <- base_pcvart_model_list[[4]]
mod_ARTPCV <- base_pcvart_model_list[[5]]
modART_pr <- base_pcvart_model_list[[6]]
modPCV_pr <- base_pcvart_model_list[[7]]
mod_ARTPCV_pr <- base_pcvart_model_list[[8]]
modART_pr_national <- base_pcvart_model_list[[15]]

### non linear effect models
mod_PCV_NL <- base_pcvart_model_list[[9]]
mod_ART_NL <- base_pcvart_model_list[[10]]
mod_PCVART_NL <- base_pcvart_model_list[[11]]
mod_PCV_prNL <- base_pcvart_model_list[[12]]
mod_ART_prNL <- base_pcvart_model_list[[13]]
mod_PCVART_prNL <- base_pcvart_model_list[[14]]

################################################################################
#### EXTRACT VARIABLE EFFECTS 
################################################################################
################ FIXED EFFECTS LINEAR - NO PROVINCE REPLICATION ##############

fixed_effects <- rbind(modART$summary.fixed[2,],
                       modPCV$summary.fixed[2,],
                       mod_ARTPCV$summary.fixed[2:3,])
fixedeffs <- cbind(fixed_effects, model = c("ART","PCV","ART_PCV","ART_PCV"), effect = c("ART","PCV","PCV","ART"))
fixedeffs$model <- factor(fixedeffs$model, levels = c("ART_PCV","ART","PCV"))

colnames(fixedeffs) <- c("mean","sd","lowerCI","median","upperCI","mode","kld","model","effect")
fixedeffs_rr <-fixedeffs
# fixedeffs[, 1:5] <- exp(fixedeffs[,1:5])
color_labels <- c("ART" = "#097969","PCV" = "#D91656", "ART_PCV" = "darkblue")

fixed_national <- ggplot(fixedeffs)+
  geom_point(aes(x=effect,y=median, group=model, color= model ),size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x=effect,ymin=lowerCI,ymax=upperCI, group=model, color= model), width=.3, position=position_dodge(width=0.5))+
  theme_bw()+
  xlab("Intervention")+
  ylab("Effect")+
  # ggtitle("Fixed Effect (no province replicate)")+
  geom_hline(yintercept=0,linetype="dashed",color="red")+
  theme(axis.text.x = element_text(angle=45,hjust=1, size=13), axis.text = element_text(size=13), axis.title = element_text(size=13)) +
  scale_color_manual(values=color_labels)+
  labs(color="Model")


################ FIXED EFFECTS LINEAR - PROVINCE REPLICATION ##############
fixed_effects <- rbind(modART_pr$summary.fixed[2,],
                       modPCV_pr$summary.fixed[2,],
                       mod_ARTPCV_pr$summary.fixed[2:3,])
fixedeffs <- cbind(fixed_effects, model = c("ART","PCV","ART_PCV","ART_PCV"), effect = c("ART","PCV","PCV","ART"))
fixedeffs$model <- factor(fixedeffs$model, levels = c("ART_PCV","ART","PCV"))

colnames(fixedeffs) <- c("mean","sd","lowerCI","median","upperCI","mode","kld","model","effect")
color_labels <- c("ART" = "#097969","PCV" = "#D91656", "ART_PCV" = "darkblue")

fixed_prov <- ggplot(fixedeffs)+
  geom_point(aes(x=effect,y=median, group=model, color= model ),size=3, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x=effect,ymin=lowerCI,ymax=upperCI, group=model, color= model), width=.3, position=position_dodge(width=0.5))+
  theme_bw()+
  xlab("Covariate")+
  ylab("Fixed Effect")+
  geom_hline(yintercept=0,linetype="dashed",color="red")+
  scale_color_manual(values=color_labels)+
  ggtitle("Fixed Effect\n(interannual effect province replicate)")+
  theme(axis.text.x = element_text(angle=45,hjust=1), axis.text = element_text(size=13),axis.title=element_text(size=13)) +
  labs(color="Model", 
       caption = "*PCV coverage data is only at the national level")
fixed_national + fixed_prov

################ NONLINEAR EFFECTS - NO PROVINCE REPLICATION ##############
color_labels <- c("ART" = "#097969","PCV" = "#D91656", "ART_PCV" = "darkblue")
nl_effects <- rbind(mod_ART_NL$summary.random[[4]],
                       mod_PCV_NL$summary.random[[4]],
                       mod_PCVART_NL$summary.random[[4]],
                    mod_PCVART_NL$summary.random[[5]])
colnames(nl_effects) <- c("coverage","mean","sd","lowerCI","median","upperCI","mode","kld")
nl_effects <- nl_effects %>%
  mutate(
    RR_mean = exp(mean),
    RR_lower = exp(lowerCI),
    RR_upper = exp(upperCI)
  )
## replicating by the number of cuts which was 5 but the PCV ones only ran 3 cuts due to the early zeroes so the model with both has 8 rows
nl_effects <- cbind(nl_effects, model = c(rep("ART",5),rep("PCV",4),rep("ART_PCV",9)), effect = c(rep("ART",5),rep("PCV",4),rep("ART",5), rep("PCV",4)))
nl_effects$model <- factor(nl_effects$model, levels = c("ART_PCV","ART","PCV"))
nl_all <- ggplot(nl_effects) +
  # geom_line(aes(x=coverage, y=median, group=model, group=effect, color=model))+
  # geom_ribbon(aes(x=coverage,ymin=lowerCI,ymax=upperCI, fill=model,group=model, group=effect),alpha=0.2)+
  geom_pointrange(aes(x=coverage, y=RR_mean, ymin=RR_lower,ymax=RR_upper, group=model, group=effect, color=model), position=position_dodge(width=0.1))+
  theme_bw()+
  xlab("Coverage")+
  ylab("Effect")+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  theme(axis.text.x = element_text(angle=45,hjust=1), axis.text = element_text(size=13),axis.title=element_text(size=13)) +
  scale_color_manual(values=color_labels)+
  scale_fill_manual(values=color_labels)+
  labs(fill="Model", color = "Model")+
  facet_grid(effect~.)
  


### just above and below 70% not all  for no province replication
# Categorize coverage as <70% or >70%
nl_effects <- nl_effects %>%
  mutate(coverage_cat = ifelse(coverage < 0.7, "<70%", ">70%"))

# Summarize data: mean, lowerCI and upperCI per group and coverage category
summary_df <- nl_effects %>%
  group_by(model, effect, coverage_cat) %>%
  summarise(
    mean_effect = mean(RR_mean),
    lower_CI = mean(RR_lower),
    upper_CI = mean(RR_upper),
    .groups = "drop"
  )

# Create a new column for x-axis labels (coverage + model)
summary_df <- summary_df %>%
  mutate(x_label = paste(coverage_cat, effect, sep = "_"))
write.csv(summary_df, file = "/home/sbelman/Documents/BRD/SouthAfrica/manuscript/figure2/coverage_cat_summarydf.csv")

# Define color palette
color_labels <- c("ART" = "#097969", "PCV" = "#D91656", "ART_PCV" = "darkblue")

# Define dodge position object
pd <- position_dodge(width = 0.6)

ggplot(summary_df, aes(x = effect, y = mean_effect, color = model)) +
  geom_point(position = pd, size = 3) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI),
                width = 0.2,
                position = pd) +
  facet_wrap(~coverage_cat) +
  scale_color_manual(values = color_labels) +
  labs(
    x = "Intervention",
    y = "Relative Risk",
    color = "Model"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  theme_bw(base_size = 15) +
  scale_y_continuous(trans = "log10" )+
  # ylim(-0.7, 0.7) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )


################ NONLINEAR EFFECTS -  PROVINCE REPLICATION ##############
color_labels <- c("ART" = "#097969","PCV" = "#D91656", "ART_PCV" = "darkblue")
nl_effects <- rbind(mod_ART_prNL$summary.random[[4]],
                    mod_PCV_prNL$summary.random[[4]],
                    mod_PCVART_prNL$summary.random[[4]],
                    mod_PCVART_prNL$summary.random[[5]])
colnames(nl_effects) <- c("coverage","mean","sd","lowerCI","median","upperCI","mode","kld")
## replicating by the number of cuts which was 5 but the PCV ones only ran 35cuts due to the early zeroes so the model with both has 8 rows
nl_effects <- cbind(nl_effects, model = c(rep("ART",5),rep("PCV",4),rep("ART_PCV",9)), effect = c(rep("ART",5),rep("PCV",4),rep("PCV",4), rep("ART",5)))
nl_effects$model <- factor(nl_effects$model, levels = c("ART_PCV","ART","PCV"))

#### pull out the intercepts so that I can calculate the incidence
nl_intercepts <- rbind(mod_ART_prNL$summary.fixed[1,c(1,3:5)],
                    mod_PCV_prNL$summary.fixed[1,c(1,3:5)],
                    mod_PCVART_prNL$summary.fixed[1,c(1,3:5)],
                    mod_PCVART_prNL$summary.fixed[1,c(1,3:5)])
colnames(nl_intercepts) <- c("intercept_mean", "intercept_lowerCI","intercept_median","intercept_upperCI")
nl_intercepts$model <- c("ART","PCV","ART_PCV","ART_PCV")
### join the intercepts
# nl_effects <- left_join(nl_effects,nl_intercepts,by="model")
# ## calculate the cases at each coverage
# nl_effects$exp_mod <- exp(nl_effects$mean)
# nl_effects$est_weekly_cases <-  nl_effects$exp_mod *500
# nl_effects$exp_int <- exp(nl_effects$intercept_mean) 
# nl_effects$est_intercept_cases <-  nl_effects$exp_int *500
# nl_effects$RRs <- nl_effects$est_weekly_cases/nl_effects$est_intercept_cases

nl_prov <- ggplot(nl_effects) +
  # geom_line(aes(x=coverage, y=median, group=model, group=effect, color=model))+
  # geom_ribbon(aes(x=coverage,ymin=lowerCI,ymax=upperCI, fill=model,group=model, group=effect),alpha=0.2)+
  geom_pointrange(aes(x=coverage, y=median, ymin=lowerCI,ymax=upperCI, group=model, group=effect, color=model), position=position_dodge(width=0.1))+
  
  theme_bw()+
  xlab("Coverage")+
  ylab("Effect")+
  geom_hline(yintercept=0,linetype="dashed",color="red")+
  theme(axis.text.x = element_text(angle=45,hjust=1), axis.text = element_text(size=13),axis.title=element_text(size=13)) +
  scale_color_manual(values=color_labels)+
  scale_fill_manual(values=color_labels)+
  labs(fill="Model", color = "Model")+
  facet_grid(effect~.)
nl_all + nl_prov


################################################################################
#### EXTRACT INTERCEPT FOR EACH MODEL
################################################################################

### intercept and estimate the actual impacts
intercept_effs <- rbind(mod_base$summary.fixed[1,], modART$summary.fixed[1,],
                        modPCV$summary.fixed[1,],
                        mod_ARTPCV$summary.fixed[1,])
intercept_effs <- cbind(intercept_effs, model = c("BASE","ART","PCV","ART_PCV"))
intercept_effs$model <- factor(intercept_effs$model, levels = c("BASE","ART_PCV","ART","PCV"))
intercept_effs$exp_mod <- exp(intercept_effs$mean)
intercept_effs$est_weekly_cases <-  intercept_effs$exp_mod *500

### intercept and estimate the actual impacts
intercept_effs_pr <- rbind(mod_base_pr$summary.fixed[1,], modART_pr$summary.fixed[1,],
                        modPCV_pr$summary.fixed[1,],
                        mod_ARTPCV_pr$summary.fixed[1,])
intercept_effs_pr <- cbind(intercept_effs_pr, model = c("BASE","ART","PCV","ART_PCV"))
intercept_effs_pr$model <- factor(intercept_effs_pr$model, levels = c("BASE","ART_PCV","ART","PCV"))
intercept_effs_pr$exp_mod <- exp(intercept_effs_pr$mean)
intercept_effs_pr$est_weekly_cases <-  intercept_effs_pr$exp_mod *500

### intercept from NL models
intercept_effs_prNL <- rbind(mod_base_pr$summary.fixed[1,], mod_ART_NL$summary.fixed[1,],
                             mod_PCV_NL$summary.fixed[1,],
                             mod_PCVART_NL$summary.fixed[1,])
intercept_effs_prNL <- cbind(intercept_effs_prNL, model = c("BASE","ART","PCV","ART_PCV"))
intercept_effs_prNL$model <- factor(intercept_effs_prNL$model, levels = c("BASE","ART_PCV","ART","PCV"))
intercept_effs_prNL$exp_mod <- exp(intercept_effs_prNL$mean)
intercept_effs_prNL$est_weekly_cases <-  intercept_effs_prNL$exp_mod *500

################################################################################
#### EXTRACT INTERANNUAL RANDOM EFFECTS
################################################################################
####### NO PROVICNE REPLICATION ###########
####### plot all the yearly effects together to compare ###########
first_year <- 2005
# year_eff_base <- rbind(modART$summary.random$id_y, modPCV$summary.random$id_y, modPCV_timing$summary.random$id_y, mod_ARTPCV$summary.random$id_y, mod_base$summary.random$id_y)
year_eff_base <- rbind(modART$summary.random$id_y, modPCV$summary.random$id_y,  mod_ARTPCV$summary.random$id_y, mod_base$summary.random$id_y)
N_years <- nrow(modART$summary.random$id_y)
colnames(year_eff_base)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld")
year_eff_base <- year_eff_base %>%
  mutate(
    RR_mean = exp(mean),
    RR_lower = exp(lowerCI),
    RR_upper = exp(upperCI)
  )

mod_names <- rep(c("ART_Coverage","PCV_Coverage","PCV_ART","Base"), each = N_years)
year_eff_base$model <- mod_names
year_eff_base$model <- factor(year_eff_base$model, levels = c("Base","ART_Coverage","PCV_Coverage", "PCV_ART"))
year_eff_base$ID <- first_year + seq(0,(N_years-1),)

lines <- c("ART_Coverage" = "longdash","PCV_Coverage" = "longdash", "PCV_ART" = "dotted", "Base" = "solid")
color_labels <- c("ART_Coverage" = "#097969","PCV_Coverage" = "#D91656", "PCV_ART" = "darkblue", "Base"= "black")
random_national <- ggplot(year_eff_base)+
  geom_hline(yintercept=1,linetype="dashed",color="grey",alpha=0.8)+
  # geom_line(aes(x=ID,y=median,linetype=model, color = model, group=model))+
  # geom_ribbon(aes(x=ID,ymin=lowerCI,ymax=upperCI, fill=model,group=model),alpha=0.2)+
  geom_line(aes(x=ID,y=RR_mean,linetype=model, color = model, group=model))+
  geom_ribbon(aes(x=ID,ymin=RR_lower,ymax=RR_upper, fill=model,group=model),alpha=0.2)+
  scale_linetype_manual(values = lines)+
  scale_color_manual(values=color_labels)+
  scale_fill_manual(values = color_labels)+
  theme_classic()+
  theme(axis.text = element_text(size=13),axis.title = element_text(size = 13), axis.text.x =element_text(angle = 45, hjust =1)
  )+
  labs(x = "Year",
       y = "Relative Risk") +
  scale_y_continuous(trans = "log10")+
  # ylim(-2,2)+
  ggtitle("Interannual Effect")
# 
# png(file = "/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/base_models/effect_artpcv.png", width=400, height=700)
# fixed_national/random_national
# dev.off()

write.csv(year_eff_base, file = "/home/sbelman/Documents/BRD/SouthAfrica/manuscript/figure2/year_eff_base.csv")

########PROVINCE REPLICATION ###########
####### plot all the yearly effects together to compare ##########
first_year <- 2005
year_eff_base <- rbind(modART_pr$summary.random$id_y, modPCV_pr$summary.random$id_y,  mod_ARTPCV_pr$summary.random$id_y, mod_base_pr$summary.random$id_y)
N_years <- nrow(modART$summary.random$id_y)
reg_year <- rep(c("Eastern Cape", "Free State", "Gauteng", 
                  "KwaZulu-Natal", "Limpopo", "Mpumalanga", 
                  "North West", "Northern Cape", "Western Cape"), each = N_years)
year_eff_base$region<-rep(reg_year, 4)

colnames(year_eff_base)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld","region")
mod_names <- rep(c("ART_Coverage","PCV_Coverage","PCV_ART","Base"), each = N_years*9)
year_eff_base$model <- mod_names
year_eff_base$model <- factor(year_eff_base$model, levels = c("Base","ART_Coverage","PCV_Coverage", "PCV_ART"))
year_eff_base$ID <- first_year + seq(0,(N_years-1),)

lines <- c("ART_Coverage" = "dotted","PCV_Coverage" = "longdash", "PCV_ART" = "longdash", "Base" = "solid")
color_labels <- c("ART_Coverage" = "#410445","PCV_Coverage" = "#D91656", "PCV_ART" = "darkblue", "Base"= "black")
random_prov <- ggplot(year_eff_base)+
  geom_hline(yintercept=0,linetype="dashed",color="grey",alpha=0.8)+
  geom_line(aes(x=ID,y=median,linetype=model, color = model, group=model, group=region))+
  # geom_ribbon(aes(x=ID,ymin=lowerCI,ymax=upperCI, fill=model,group=model),alpha=0.2)+
  scale_linetype_manual(values = lines)+
  scale_color_manual(values=color_labels)+
  scale_fill_manual(values = color_labels)+
  theme_classic()+
  theme(axis.text = element_text(size=13), axis.text.x = element_text(angle=45,hjust=1),axis.title = element_text(size = 13))+
  labs(x = "Year",
       y = "Median Effect") +
  # ylim(-2,2)+
  ggtitle("Interannual Effect Replicated by Province") +
  facet_wrap(.~region)

random_national + random_prov


################### PROVINCE WITH NONLINEAR MODEL ########################
first_year <- 2005
year_eff_base <- rbind(mod_ART_prNL$summary.random$id_y, mod_PCV_prNL$summary.random$id_y,  mod_PCVART_prNL$summary.random$id_y, mod_base_pr$summary.random$id_y)
N_years <- nrow(modART$summary.random$id_y)
reg_year <- rep(c("Eastern Cape", "Free State", "Gauteng", 
                  "KwaZulu-Natal", "Limpopo", "Mpumalanga", 
                  "North West", "Northern Cape", "Western Cape"), each = N_years)
year_eff_base$region<-rep(reg_year, 4)

colnames(year_eff_base)<-c("ID","mean","sd","lowerCI","median","upperCI","mode","kld","region")
mod_names <- rep(c("ART_Coverage","PCV_Coverage","PCV_ART","Base"), each = N_years*9)
year_eff_base$model <- mod_names
year_eff_base$model <- factor(year_eff_base$model, levels = c("Base","ART_Coverage","PCV_Coverage", "PCV_ART"))
year_eff_base$ID <- first_year + seq(0,(N_years-1),)

lines <- c("ART_Coverage" = "dotted","PCV_Coverage" = "longdash", "PCV_ART" = "longdash", "Base" = "solid")
color_labels <- c("ART_Coverage" = "#410445","PCV_Coverage" = "#D91656", "PCV_ART" = "darkblue", "Base"= "black")
random_prov_NL <- ggplot(year_eff_base)+
  geom_hline(yintercept=0,linetype="dashed",color="grey",alpha=0.8)+
  geom_line(aes(x=ID,y=median,linetype=model, color = model, group=model, group=region))+
  # geom_ribbon(aes(x=ID,ymin=lowerCI,ymax=upperCI, fill=model,group=model),alpha=0.2)+
  scale_linetype_manual(values = lines)+
  scale_color_manual(values=color_labels)+
  scale_fill_manual(values = color_labels)+
  theme_classic()+
  theme(axis.text = element_text(size=13), axis.text.x = element_text(angle=45,hjust=1),axis.title = element_text(size = 13))+
  labs(x = "Year",
       y = "Median Effect") +
  # ylim(-2,2)+
  ggtitle("Interannual Effect Replicated by Province (Nonlinear)") +
  facet_wrap(.~region)

random_prov_NL + random_prov




