################## test spatial replicate by year versus interannual replicate by province ##############
####LOAD DATA##########
################################################################################
source("/home/sbelman/Documents/env_sa_manuscript/scripts2/0_source_functions.R")
# weekly=FALSE
precov=TRUE ### set whether the run ends in 2020 or 2023
time = "weekly"
space="adm2"
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
hurs_grp <- INLA::inla.group(df$hurs_lag3, n = 5)  # You can adjust 'n' (number of groups) depending on granularity
df$hurs_grp <- hurs_grp
absh_grp <- inla.group(df$absh_lag3, n = 5)  # You can adjust 'n' (number of groups) depending on granularity
df$absh_grp <- absh_grp

df_all <- df


############ CROSSBASIS ############
cov_oi <- "pm2p5_lag0"
max_lag <- 8
lag_knots <- logknots(max_lag, 2)
exp_knots <- quantile(df[[cov_oi]], probs = c(0.33, 0.66), na.rm = TRUE)

cb <- crossbasis(
  x = df_all[[cov_oi]],
  lag = max_lag,
  argvar = list(fun = "ns", knots = exp_knots),
  arglag = list(fun = "ns", df = 2),
  # arglag = list(fun = "ns", knots = lag_knots),
  group = df$id_u
)
######################################################################
dftest <- cbind(df_all, cb)
dftest$id_y2 <- dftest$id_y
################# run base formula used later ####
#### test
model_specs <- list(
  
  ybyprov = list(
    spatial_term = 'f(id_u, model = "bym2", graph = g, scale.model = T, adjust.for.con.comp=TRUE,
                      hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    
    yearly_term  = 'f(id_y, model = "iid", replicate = id_prov,
                      hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    
    desc = "spatial, seasonal, and annual (province replication)"
  ),
  
  spacebyyear = list(
    spatial_term = 'f(id_u, model = "bym2", replicate = id_y2, graph = g, scale.model = T,
                      adjust.for.con.comp=TRUE,
                      hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    
    # yearly_term  = 'f(id_y, model = "iid",
    #                   hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',

    # desc = "spatial (yearly replication) and seasonal"
    desc = "spatial (yearly replication), annual, and seasonal"
    
  )
)

run_inla_model <- function(spec, cb, data) {
  
  if(spec$desc == "spatial, seasonal, and annual (province replication)"){
  form <- reformulate(c(
    1,
    paste(colnames(cb), collapse = " + "),
    spec$spatial_term,
    'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T,
       hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
    spec$yearly_term,
    'vaccination_period',
    'population_density',
    'tas_lag0',
    'f(absh_grp, model = "rw2", scale.model = T,
       hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))'
  ), response = "disease")
  }else{
    form <- reformulate(c(
      1,
      paste(colnames(cb), collapse = " + "),
      spec$spatial_term,
      'f(id_m, model = "rw2", cyclic = T, scale.model = T, constr = T,
       hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))',
      # spec$yearly_term,
      'vaccination_period',
      'population_density',
      'tas_lag0',
      'f(absh_grp, model = "rw2", scale.model = T,
       hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))'
    ), response = "disease")
  }
  
  fit <- INLA::inla(
    formula = as.formula(form),
    family = "nbinomial",
    offset = log(population_size/100000),
    control.inla = list(strategy = 'adaptive'),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
                           config = FALSE, return.marginals = TRUE),
    control.predictor = list(link = 1, compute = TRUE),
    control.fixed = list(correlation.matrix = TRUE),
    inla.setOption(num.threads = 4),
    verbose = FALSE,
    data = data
  )
  
  fit$cov <- paste0(spec$desc," accounting for vaccination and population density")
  
  return(fit)
}

models <- lapply(model_specs, run_inla_model, cb = cb, data = dftest)


############## EXTRACT COEFFICIENTS TO COMPARE DIFFERENCE TO PM2.5 #############
extract_coefs <- function(mod, cb) {
  cb_names <- colnames(cb)
  idx <- match(cb_names, rownames(mod$summary.fixed))
  if(any(is.na(idx))) {
    stop("Crossbasis terms missing in model")
  }
  coef <- mod$summary.fixed$mean[idx]
  # Try preferred first
  vcov <- tryCatch(
    mod$misc$lincomb.derived.covariance[idx, idx],
    error = function(e) mod$misc$cov.fixed[idx, idx]
  )
  return(list(coef = coef, vcov = vcov))
}

coefs <- lapply(models, extract_coefs, cb = cb)

cen <- 40

preds <- lapply(coefs, function(x) {
  crosspred(
    basis = cb,
    cen = cen,
    coef = x$coef,
    vcov = x$vcov
  )
})

beta_base <- coefs$ybyprov$coef
beta_flex <- coefs$spacebyyear$coef

pct_diff <- (beta_flex - beta_base) / beta_base * 100

extract_curve <- function(cp, model_name) {
  data.frame(
    x   = cp$predvar,
    fit = exp(cp$allfit),     # convert to RR
    low = exp(cp$alllow),
    high= exp(cp$allhigh),
    model = model_name
  )
}

curve_df <- bind_rows(
  extract_curve(preds$ybyprov, "interannual rep by province"),
  extract_curve(preds$spacebyyear, "spatial (district) rep by year")
)

cumpm25diff <- ggplot(curve_df, aes(x = x, y = fit, color = model, fill = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.4)+
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, color = NA) +
  geom_line() +
  labs(y = "Cumulative Relative Risk",
    title = "Exposure-Response Curves by Model Specification"
  ) +
  xlab(expression(paste('Concentration PM'[2.5],' (', mu, 'g/m'^3, ')')))+
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())


pdf("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/sensitivity/spatial_sensitivity/cumpm2p5_spatialsens_weekly_adm2_2019.pdf", width = 4, height =4)
print(cumpm25diff)
dev.off()

ggplot(curve_df[which(curve_df$x%in%c(50,100)),]) +
  geom_pointrange(aes(x=as.character(x), y =fit, ymin=low, ymax=high, group = model, color =model), position = position_dodge(width = 0.3))+
  theme_bw()

models$ybyprov$summary.hyperpar
models$spacebyyear$summary.hyperpar

pmsum <- df %>%
  group_by(GID_2, year) %>%
  summarise(meanpm = mean(pm2p5_lag0), .groups = "drop")

#### examine the spatial effects in detail and their deviation ##############################
###### plot #### spatial effects  
plotinterrep <- plot_spatial_effects(models$ybyprov,shp,N=52, structured = FALSE)
pdf("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/sensitivity/spatial_sensitivity/spatialstatic_yearbyprov_weekly_adm2_2019.pdf", width = 6, height =4)
print(plotinterrep)
dev.off()

library(sf)
# Define constants
N_areas <- 52
years <- 2005:2019
N_years <- length(years)

# Extract and label the random effects
res_all <- models$spacebyyear$summary.random$id_u %>%
  mutate(# Each year block is 104 rows (52 unstructured + 52 structured)
    year = rep(years, each = N_areas * 2),
    # Within each year, the first 52 are 'Unstructured', the next 52 are 'Structured'
    type = rep(rep(c("Unstructured", "Structured"), each = N_areas), N_years),
    # The actual area ID (1 to 52) repeated for both types across all years
    area_idx = rep(1:N_areas, N_years * 2))

# Join with your shapefile
# Ensure shp is ordered correctly to match the 1:52 ID
shp_indexed <- shp %>% 
  arrange(GID_2) %>% 
  mutate(area_idx = 1:n())

map_data <- shp_indexed %>%
  left_join(res_all, by = "area_idx")

spacebyyear <- ggplot(data = filter(map_data, type == "Structured")) +
  geom_sf(aes(fill = mean), color = "black", size = 0.05) +
  scale_fill_gradient2(low = "blue", mid = "white",high = "red", midpoint = 0,name = "Structured Effect", limits = c(-4,4) ) +
  facet_wrap(~year, ncol = 5) +
  theme_minimal() +
  labs(title = "Evolution of Spatial Structured Effects",
    subtitle = "District-level BYM2 Model (2005-2019)") +
  theme(axis.text = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"))
pdf("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/sensitivity/spatial_sensitivity/spatialdynamic_weekly_adm2_2019.pdf", width = 8, height =5)
print(spacebyyear)
dev.off()

### now plot the difference for each year as compared to the model without hte yearly replication for the spatial effect
# Extract the 104 rows from the province/baseline model
# Assuming the first 52 are unstructured and the next 52 are structured
N_areas <- 52
shp$province <- shp$NAME_1
area_province <- shp %>%
  arrange(GID_2) %>%
  mutate(area_idx = 1:n(),
         province_idx = as.numeric(factor(NAME_1))) %>%
  dplyr::select(area_idx, province_idx)


## pull res idu
res_idu <- models$ybyprov$summary.random$id_u %>%
  slice(1:N_areas) %>%   # N_areas = 52
  mutate(area_idx = 1:N_areas) %>%
  dplyr::select(area_idx, mean_idu = mean)
# pull out res yearly
res_idy <- models$ybyprov$summary.random$id_y %>%
  mutate(year = rep(2005:2019, each = 9),  # 9 provinces
         province_idx = rep(unique(area_province$province), times = 15)) %>%
  dplyr::select(year, province_idx, mean_idy = mean)
res_idy_expanded <- area_province %>%
  left_join(res_idy, by = "province_idx")
res_combined <- res_idy_expanded %>%
  left_join(res_idu, by = "area_idx") %>%
  mutate(mean_total = mean_idu + mean_idy)

## old way
# res_baseline <- models$ybyprov$summary.random$id_u %>%
#   slice((N_areas + 1):(N_areas * 2)) %>% # Get only Structured effects
#   mutate(area_idx = 1:N_areas) %>%
#   dplyr::select(area_idx, mean_baseline = mean)

# 1. Process the Yearly Replicated Model (the 1560 rows)
res_yearly <- models$spacebyyear$summary.random$id_u %>%
  mutate(year = rep(2005:2019, each = N_areas * 2),
    type = rep(rep(c("Unstructured", "Structured"), each = N_areas), 15),
    area_idx = rep(1:N_areas, 15 * 2)) %>%
  filter(type == "Structured") %>%
  dplyr::select(year, area_idx, mean_yearly = mean)

# 2. Join the Baseline and Calculate Difference
res_diff <- res_yearly %>%
  left_join(res_combined, by = c("area_idx","year")) %>%
  mutate(diff = mean_yearly - mean_total)

# 3. Merge with Shapefile
map_diff_data <- shp %>%
  arrange(GID_2) %>%
  mutate(area_idx = 1:n()) %>%
  left_join(res_diff, by = "area_idx")

spacbyyeardiff <- ggplot(data = map_diff_data) +
  geom_sf(aes(fill = diff), color = "black", size = 0.05) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,name = "Yearly Deviation", limits = c(-4,4)) +
  facet_wrap(~year, ncol = 5) +
  theme_minimal() +
  labs(title = "Yearly Spatial Deviations from Static Spatial",
    subtitle = "Positive (Red) indicates higher risk than the yearly average for that district",
    caption = "Difference = Yearly Replicated Spatial Effect - Static Spatial Effect") +
  theme(axis.text = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"))
pdf("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/sensitivity/spatial_sensitivity/diffspacebyyear_idmidy_weekly_adm2_2019.pdf", width = 8, height =5)
print(spacbyyeardiff)
dev.off()



provincedeviations <- ggplot(map_diff_data, aes(x = year, y = diff, group = district)) +
  # Add a light grey line for every district to show the bundle
  geom_line(alpha = 0.3, color = "grey40") + 
  # Add a zero line to represent the "Global Provincial Baseline"
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  # Facet by Province - let each province have its own Y-axis if needed
  facet_wrap(~NAME_1, scales = "free_y", ncol = 3) + 
  theme_bw(base_size = 11) +
  labs(title = "Dynamic Spatial Effect Deviations",
    x = "Year",
    y = "Difference in Random Effect (Dynamic - Static)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 11), axis.title = element_text(size = 11))


provinceeffects <- ggplot(map_diff_data, aes(x = year, y = mean_yearly, group = district)) +
  # Add a light grey line for every district to show the bundle
  geom_line(alpha = 0.3, color = "grey40") + 
  # Add a zero line to represent the "Global Provincial Baseline"
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  # Facet by Province - let each province have its own Y-axis if needed
  facet_wrap(~NAME_1, scales = "free_y", ncol = 3) + 
  theme_bw(base_size = 11) +
  labs(title = "Dynamic Spatial Effect",
       x = "Year",
       y = "Annual Spatial Effect") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 11), axis.title = element_text(size = 11))

############# plot interannual random effect from other model ############
yp <- models$ybyprov$summary.random$id_y
yp$year <- 2004 + yp$ID
year_n <- nrow(yp)/9 
colnames(yp)[grep("quant",colnames(yp))] <- c("lowerCI","median","upperCI")
reg_month_year <- rep(c("Eastern_Cape", "Free_State", "Gauteng",
                        "KwaZulu-Natal", "Limpopo", "Mpumalanga",
                        "North_West", "Northern_Cape", "Western_Cape"), each = year_n)
yp$region <- factor(reg_month_year)
y <- ggplot(yp)+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  geom_ribbon(aes(x = year, ymin = lowerCI, ymax = upperCI, group = region, ), alpha = 0.3)+
  geom_line(aes(x = year, y = median, group=region))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 11), axis.title = element_text(size = 11))+
  xlab("Year")+
  ylab("Annual Effect")+
  labs(fill="Province",color = "Province", title = "Interannual Effect")+
  facet_wrap(region ~.)

library(patchwork)
pdf("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/sensitivity/spatial_sensitivity/diffspatial_province_idmidy_weekly_adm2_2019.pdf", width = 10, height =5)
print(provincedeviations + y)
dev.off()

pdf("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/sensitivity/spatial_sensitivity/spatialinterann_province_weekly_adm2_2019.pdf", width = 10, height =5)
print(y + provinceeffects )
dev.off()

pdf("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/sensitivity/spatial_sensitivity/diffspatial_provincewrap_weekly_adm2_2019.pdf", width = 5, height =5)
print(provincedeviations)
dev.off()


#### calculate the correlation between these effects
# Summarize your data to get the provincial mean per year
map_diff_dataonly <- st_drop_geometry(map_diff_data)
prov_data <- map_diff_dataonly %>%
  group_by(NAME_1,NAME_2, year) %>%
  summarise(provest = mean(mean_yearly, na.rm = TRUE), .groups = "drop")

# Suppose yp$mean is aligned by year for each province (needs to be the same length)
# Join if necessary
colnames(yp)[which(colnames(yp)=="region")] <- "NAME_1"
yp$NAME_1 <- gsub("_"," ",yp$NAME_1)
combined <- prov_data %>%
  left_join(yp, by = c("year","NAME_1"))  # assuming yp has a 'year' column and 'mean' column

# Calculate correlation for each province
cor_results <- combined %>%
  group_by(NAME_1,NAME_2) %>%
  summarise(
    test = list(cor.test(provest, mean, use = "complete.obs")),
    correlation = map_dbl(test, ~ .x$estimate),
    conf.low = map_dbl(test, ~ .x$conf.int[1]),
    conf.high = map_dbl(test, ~ .x$conf.int[2])
  ) %>%
  dplyr::select(-test)

ggplot(cor_results, aes(x = correlation, y = reorder(NAME_2, correlation))) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ NAME_1, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Correlation (95% CI)",
    y = "District"
  )
shpcor <- left_join(shp,cor_results,by="NAME_2")
ggplot(shpcor) +
  geom_sf(aes(fill = correlation), color = "black", size = 0.2) +
  scale_fill_viridis_c(option = "plasma", direction = 1, 
                       limits = c(-1,1), na.value = "grey80") +
  theme_minimal() +
  labs(fill = "Correlation")+

# ############# compare base models at each specification #############
# dftest <- df_all; dftest$id_y2 <- dftest$id_y
# 
# model_specs <- list(
#   ybyprov = list(
#     sp='f(id_u,model="bym2",graph=g,scale.model=T,adjust.for.con.comp=T)',
#     yr='f(id_y,model="iid",replicate=id_prov)'),
#   spacebyyear = list(
#     sp='f(id_u,model="bym2",replicate=id_y2,graph=g,scale.model=T,adjust.for.con.comp=T)')
#   # yr='f(id_y,model="iid")')
# )
# 
# run <- function(s)
#   INLA::inla(
#     as.formula(paste("disease ~ 1 +", s$sp,
#                      '+ f(id_m,model="rw2",cyclic=T,scale.model=T,constr=T) +',
#                      s$yr, '+ vaccination_period + population_density')),
#     family="nbinomial",
#     offset=log(population_size/100000),
#     data=dftest,
#     control.compute=list(dic=T,waic=T, cpo = T))
# 
# mods <- lapply(model_specs, run)
# 
# 
# int_mod <- readRDS(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/base_models/base_model_main_",time,"_intercept_",space,"_",endyear,".rds"))
# gof <- eval.mod(mods$ybyprov,dftest)
# gof$rsq <- rsq(gof, int_mod, 1)
# gof$type <- "ybprov"
# 
# gof1 <- eval.mod(mods$spacebyyear,dftest)
# gof1$rsq <- rsq(gof1, int_mod, 1)
# gof1$type <- "spacebyyear"
# rbind(gof,gof1)



################################################################################
## goodness of fit
rbind(eval.mod(base_ybp,dftest),
      eval.mod(base_sby,dftest))

### plot
n_districts <- 52
n_years <- 15
n_components <- 2   # structured + unstructured

re <- base_sby$summary.random$id_u

re$component <- rep(rep(c("structured", "unstructured"),
                        each = n_districts),
                    times = n_years)

re$year <- rep(1:n_years, each = n_districts * n_components)

re$district_id <- rep(1:n_districts, times = n_years * n_components)

re_total <- re %>%
  group_by(year, district_id) %>%
  summarise(mean = sum(mean),
            sd = sqrt(sum(sd^2)), .groups = "drop")

shp$district_id <- 1:nrow(shp)

plot_data <- re[which(re$component=="structured"),] %>%
  left_join(shp, by = "district_id") %>%
  sf::st_as_sf()
ggplot(plot_data) +
  geom_sf(aes(fill = mean), color = NA) +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~ year, ncol = 5) +
  theme_minimal() +
  labs(fill = "Spatial effect",
       title = "BYM2 Total Spatial Effect by Year")
