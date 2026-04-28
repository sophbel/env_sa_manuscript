library(dlnm)
library(INLA)

setwd("/home/sbelman/Documents/env_sa_manuscript/")
####LOAD DATA & LIBRARIES #####################################################
source("scripts2/0_source_functions.R")
### set if timing is precovid is true or not
precov = TRUE
### set resolution
time = "weekly"
space = "adm2"

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

##### LOAD INTERCEPT MODELS FOR R2 CALCULATION ##################################
int_mod <- readRDS(file=paste0("models/base_models/base_model_main_",time,"_intercept_",space,"_",endyear,".rds"))
re_mod <- readRDS(file=paste0("models/base_models/base_model_",time,"_20092011_popdens_",space,"_",endyear,".rds"))

dlnm_results <- NULL

### -----------------------------
### 1. DEFINE LAG STRUCTURE
### -----------------------------
max_lag <- ifelse(time == "weekly", 8, 3)

### -----------------------------
### 2. CROSSBASIS: PM2.5
### -----------------------------
pm_knots <- quantile(df$pm2p5_lag0, probs = c(0.33, 0.66), na.rm = TRUE)

cb_pm <- crossbasis(
  x      = df$pm2p5_lag0,
  lag    = max_lag,
  argvar = list(fun = "ns", knots = pm_knots),
  arglag = list(fun = "ns", df = 2),
  group  = df$id_u
)

### -----------------------------
### 3. CROSSBASIS: TEMPERATURE
### -----------------------------
temp_knots <- quantile(df$tas_lag0, probs = c(0.33, 0.66), na.rm = TRUE)

cb_temp <- crossbasis(
  x      = df$tas_lag0,
  lag    = max_lag,
  argvar = list(fun = "ns", knots = temp_knots),
  arglag = list(fun = "ns", df = 2),
  group  = df$id_u
)

### -----------------------------
### 4. COMBINE DATA
### -----------------------------
df_model <- cbind(df, cb_pm, cb_temp)

### -----------------------------
### 5. DEFINE MODEL FORMULA
### -----------------------------
formula_dlnm <- as.formula(
  paste(
    "disease ~",
    paste(colnames(cb_pm), collapse = " + "), "+",
    paste(colnames(cb_temp), collapse = " + "), "+",
    
    # Spatial effect
    "f(id_u, model = 'bym2', graph = g,
       scale.model = TRUE, adjust.for.con.comp = TRUE,
       hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))) +",
    
    # Seasonal monthly trend
    "f(id_m, model = 'rw2', cyclic = TRUE,
       scale.model = TRUE, constr = TRUE,
       hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))) +",
    
    # Year effect
    "f(id_y, model = 'iid', replicate = id_prov,
       hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))) +",
    
    # Other fixed covariates
    "vaccination_period +",
    "population_density +",
    
    # Optional humidity smooth (if desired)
    "f(absh_grp, model = 'rw2',
       scale.model = TRUE,
       hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))"
  )
)

### -----------------------------
### 6. RUN INLA MODEL
### -----------------------------
mod <- INLA::inla(
  formula = formula_dlnm,
  family  = "nbinomial",
  data    = df_model,
  offset  = log(population_size / 100000),
  control.inla = list(strategy = "adaptive"),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  control.predictor = list(compute = TRUE),
  control.fixed = list(correlation.matrix = TRUE),
  verbose = FALSE
)

summary(mod)

### -------------------------------------------------
### 7. EXTRACT PM COEFFICIENTS FOR CROSSPRED
### -------------------------------------------------

# Identify PM coefficient names
pm_coef_names <- colnames(cb_pm)

# Extract coefficients
coef_pm <- mod$summary.fixed[pm_coef_names,c("mean")]

# Extract variance-covariance matrix
vcov_full <- mod$misc$lincomb.derived.covariance.matrix
vcov_pm   <- vcov_full[pm_coef_names, pm_coef_names]

### -------------------------------------------------
### 8. RUN CROSSPRED FOR PM
### -------------------------------------------------

# Set centering value (WHO guideline for PM2.5 example)
cen_pm <- 40   # change if needed

cp_pm <- crosspred(
  basis = cb_pm,
  coef  = coef_pm,
  vcov  = vcov_pm,
  cen   = cen_pm
)

summary(cp_pm)

### extract fit
df_all <- extract_cp_gpsc_data(cp_pm,"none","no_interaction", "pm2p5")

### extract goodness of fit metrics
cov<-mod$cov
#gof
gof <- eval.mod(mod,df)
gof$rsq <- rsq(gof, int_mod, 1)
gof <- gof[,c("cov","waic","mae","cpo","rsq")]

write.table(gof, file=paste0("models/dlnms/modifiers/mod_gof_dlnm2X_tasabsh",time,"_",space,"_",max_lag,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")
write.table(df_all, file=paste0("models/dlnms/modifiers/nointeraction_dlnm2x_results_fits_",time,"_",space,"_tasabsh_",max_lag,"_",endyear,".csv"), quote = FALSE, col.names = TRUE, row.names = TRUE, sep = ",")



# ##### compare CBpm2p5 + CBtas + absh model TO CBpm2p5 + tas + absh model
# dfth <- fread("models/dlnms/bivariable/nointeraction_results_fits_weekly_adm2_bivariable_8_2019.csv")
# dfth <- subset(dfth , dfth$covariate == "pm2p5_lag0" & dfth$cov2 == "absh_lag3")
# dfthC <- subset(dfth , dfth$lag_num == 3)
# ## read in double cross basis
# dfCB2 <- fread("models/dlnms/modifiers/nointeraction_dlnm2x_results_fits_weekly_adm2_tasabsh_8_2019.csv")
# dfCB2 <- subset(dfCB2 , dfCB2$covariate == "pm2p5" )
# dfCB2C <- subset(dfCB2 , dfCB2$covariate == "pm2p5" & dfCB2$lag_num == 3)
# 
# ## interaction version
# dfint <- fread("models/dlnms/modifiers/envmod_results_fits_weekly_adm2_allmodifiers_propprov_8_2019.csv")
# dfint <- subset(dfint , dfint$covariate == "pm2p5_lag0" &  dfint$interaction_level =="med")
# dfintC <- subset(dfint , dfint$lag_num == 3 )
# 
# ggplot() +
#   geom_line(data = dfthC, aes(x=predvar, y = cumulative_fit), color = "darkgreen") +
#   geom_ribbon(data=dfthC, aes(x=predvar, ymin = cum_lowerCI, ymax = cum_upperCI), fill = "darkgreen", alpha = 0.2) +
#   geom_line(data = dfCB2C, aes(x=predvar, y = cumulative_fit), color = "darkblue") +
#   geom_ribbon(data=dfCB2C, aes(x=predvar, ymin = cum_lowerCI, ymax = cum_upperCI), fill = "darkblue", alpha = 0.2) +
#   # geom_line(data = dfintC, aes(x=predvar, y = cumulative_fit), color = "red") +
#   # geom_ribbon(data=dfintC, aes(x=predvar, ymin = cum_lowerCI, ymax = cum_upperCI), fill = "red", alpha = 0.2) +
#   theme_bw()
# 
# 
# #### assess colinearity
# temp_lag3 <- df$tas_lag3
# 
# # extract crossbasis matrix only
# cb_temp_mat <- as.matrix(cb_temp)
# 
# # correlation between lag3 and each CB column
# cors <- apply(cb_temp_mat, 2, function(x)
#   cor(x, temp_lag3, use = "complete.obs")
# )
# round(cors, 3)
# lm_test <- lm(temp_lag3 ~ cb_temp_mat)
# summary(lm_test)$r.squared

# ggplot() +
#   geom_ribbon(data = dfth,aes(x = predvar,ymin = lowerCI,ymax = upperCI, group = lag_num),fill = "darkgreen",alpha = 0.2) +
#   geom_line(data = dfth,aes(x = predvar,y = fit, group = lag_num),color = "darkgreen") +
#   geom_ribbon(data = dfCB2,aes(x = predvar,ymin = lowerCI,ymax = upperCI, group = lag_num),fill = "darkblue",alpha = 0.2) +
#   geom_line(data = dfCB2,aes(x = predvar,y = fit, group = lag_num),color = "darkblue") +
#   geom_line(data = dfint, aes(x=predvar, y = fit, group = lag_num), color = "red") +
#   geom_ribbon(data=dfint, aes(x=predvar, ymin = lowerCI, ymax = upperCI, group = lag_num), fill = "red", alpha = 0.2) +
#   labs(x = "PM2.5",
#        y = "Effect (log RR)",
#        title = "") +
#   theme_bw() +
#   facet_grid(.~lag_num)+
#   theme(plot.title = element_text(hjust = 0.5))

