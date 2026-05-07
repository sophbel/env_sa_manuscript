################################################################################
#### DLNM SENSITIVITY ANALYSIS: WEEKLY ADM2 ####################################
################################################################################
setwd("/home/sbelman/Documents/env_sa_manuscript/")
source("scripts2/0_source_functions.R")

# Fixed Settings
time <- "weekly"
space <- "adm2"
precov <- TRUE

# Load Data
if(time == "weekly" & space == "adm2"){
  data <- fread(file="dataframes/sa_adm2_weekly_lag_sc.csv")
  data_unscaled <- fread(file="dataframes/sa_adm2_weekly_lag.csv")
  shp<-st_read("input_datasets/shps/gadm41_namematch_ZAF_2.shp")
  ## read adjacency matrix
  g <- inla.read.graph(filename = "input_datasets/shps/sa_adjacency_map.adj")
}

# Pre-processing
df <- data %>% mutate(date = as.Date(week))
df$id_prov <- as.numeric(factor(df$NAME_1))
df$vaccination_period <- as.factor(ifelse(df$date < as.Date("2009-01-01"), 1, 
                                          ifelse(df$date >= as.Date("2011-01-01"), 3, 2)))

if(precov) df <- subset(df, date < as.Date("2020-01-01"))

# Define variables to test
cov_names <- grep("tasmax|tas|tasmin|hurs|absh|prlrsum|prlrmean|sfcWind|spei3|spei6", 
                  colnames(df), value = TRUE)
cn <- grep("tas|tasmin|tasmax|hurs|absh|prlrsum|sfcWind|spei3|pm2p5|pm10|so2", 
                  colnames(df), value = TRUE)
cov_names <- grep("lag0", 
                  cn, value = TRUE)

# Load Intercept Models for R2
int_mod <- readRDS(file=paste0("models/base_models/base_model_main_weekly_intercept_adm2_2019.rds"))

# Define Sensitivity Scenarios
sens_grid <- list(
  "base1" = list(kv = 2, kl = 2), # DF for lag for everything else
  "base2" = list(kv = 2, kl = 3), # DF for lag for hurs
  "varflex" = list(kv = 3, kl = 3), # Increased DF for variable
  "lagflex" = list(kv = 2, kl = 5) # Increased DF for lag
)

# Results Storage
gof_results <- data.frame()
plot_data <- data.frame()

################################################################################
#### LOOP THROUGH VARIABLES AND SPECIFICATIONS #################################
################################################################################

for(c in seq_along(cov_names)) {
  cov_oi <- cov_names[c]
  unscaled_var <- data_unscaled[[cov_oi]]
  
  # Determine Center (Median)
  cen <- median(df[[cov_oi]], na.rm = TRUE)
  # For unscaling the plot later
  var_mean <- mean(unscaled_var, na.rm = TRUE)
  var_sd <- sd(unscaled_var, na.rm = TRUE)
  
  for(s in seq_along(sens_grid)) {
    spec_name <- names(sens_grid)[s]
    kv <- sens_grid[[s]]$kv
    kl <- sens_grid[[s]]$kl
    
    print(paste("Testing:", cov_oi, "| Spec:", spec_name))
    
    # Define Var Knots
    if(grepl("prlr", cov_oi)) {
      exp_knots <- logknots(max(df[[cov_oi]], na.rm=T), kv)
    } else {
      probs <- seq(1/(kv+1), kv/(kv+1), length.out = kv)
      exp_knots <- quantile(df[[cov_oi]], probs = probs, na.rm = TRUE)
    }
    
    # Create Crossbasis
    cb <- crossbasis(
      x = df[[cov_oi]],
      lag = 8,
      argvar = list(fun = "ns", knots = exp_knots),
      arglag = list(fun = "ns", df = kl),
      group = df$id_u
    )
    
    # Prepare Model Data
    cb_df <- cbind(df, cb)
    
    # Formula
    cb_form <- as.formula(paste(
      "disease ~", paste(colnames(cb), collapse = " + "),
      "+ f(id_u, model = 'bym2', graph = g, scale.model = T)",
      "+ f(id_m, model = 'rw2', cyclic = T, scale.model = T)",
      "+ f(id_y, model = 'iid', replicate = id_prov)",
      "+ vaccination_period + population_density"
    ))
    
    # Run INLA
    # mod <- INLA::inla(
    #   formula = cb_form, family = "nbinomial",
    #   offset = log(population_size/100000),
    #   control.compute = list(dic = TRUE, waic = TRUE),
    #   data = cb_df
    # )
    
    print("running INLA")
    mod <- INLA::inla(
      formula = cb_form,
      family = "nbinomial",  # One per outcome
      offset=log(population_size/100000),
      control.inla = list(strategy = 'adaptive'),
      control.compute = list(dic = TRUE, waic=TRUE, cpo=TRUE, config = FALSE, return.marginals = TRUE),
      control.predictor = list(link = 1, compute = TRUE),
      control.fixed = list(correlation.matrix = TRUE),
      inla.setOption(num.threads = 4),
      verbose = FALSE,
      data = cb_df)
    
    # 1. Goodness of Fit
    gof_metrics <- eval.mod(mod, cb_df)
    current_rsq <- rsq(gof_metrics, int_mod, 1)
    gof_results <- rbind(gof_results, data.frame(
      variable = cov_oi,
      spec = spec_name,
      kv = kv,
      kl = kl,
      waic = mod$waic$waic,
      dic = mod$dic$dic,
      deviance_dic = gof_metrics$deviance_dic,
      mae = gof_metrics$mae,
      cpo = gof_metrics$cpo,
      rsq = current_rsq,
      n_obs = gof_metrics$N
    ))
    # gof_results <- rbind(gof_results, data.frame(
    #   variable = cov_oi,
    #   spec = spec_name,
    #   waic = mod$waic$waic,
    #   dic = mod$dic$dic
    # ))
    
    # 2. Extract Cumulative Curves for Plotting
    # pred <- crosspred(cb, mod, cen = cen)
    original_coefs <- mod$summary.fixed$mean[c(1:ncol(cb)+1)]  # Extract fixed-effect coefficients
    original_vcov <- mod$misc$lincomb.derived.covariance[1:ncol(cb)+1,1:ncol(cb)+1] # Extract variance-covariance
    pred <- crosspred(basis = cb, cen=cen, coef = original_coefs, vcov = original_vcov) ## have to center it somewhere because its not linear so its not just at time 0.
    
    
    # Rescale x-axis for plotting
    real_x <- (pred$predvar * var_sd) + var_mean
    
    temp_plot <- data.frame(
      variable = cov_oi,
      spec = spec_name,
      x = real_x,
      rr = exp(pred$allfit),
      lower = exp(pred$alllow),
      upper = exp(pred$allhigh)
    )
    plot_data <- rbind(plot_data, temp_plot)
    
    rm(cb, mod, pred); gc() # Clean memory
  }
}

################################################################################
#### PLOTTING & EXPORT #########################################################
################################################################################
# Save GOF table
write.csv(gof_results, "models/dlnms/univariable/sensitivity/args_sensitivity_gof_results.csv", row.names = FALSE, quote = FALSE)

# Generate Sensitivity Plots
# One plot per variable showing the three cumulative curves
library(ggplot2)

for(v in cov_names) {
  p_subset <- subset(plot_data, variable == v)
  
  p <- ggplot(p_subset, aes(x = x, y = rr, color = spec)) +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = spec), alpha = 0.1, color = NA) +
    facet_wrap(~variable, scales = "free_x") +
    scale_y_continuous(trans = "log10") +
    labs(title = paste("Sensitivity Analysis:", v),
         subtitle = "Comparison of Cumulative Relative Risk (Lag 0-8)",
         x = "Exposure Value", y = "Cumulative RR (log scale)",
         fill = "Specification", color = "Specification") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(filename = paste0("figures/supplement/arglag_sensitivity/args_sensitivity_", v, ".png"), plot = p, width = 8, height = 6)
}

print("Sensitivity Analysis Complete. Plots and GOF saved.")



library(ggplot2)
library(dplyr)
library(patchwork)

### plot again but in panel
met_vars <- c("tas_lag0","tasmax_lag0","tasmin_lag0","hurs_lag0",
              "absh_lag0","prlrsum_lag0","sfcWind_lag0","spei3_lag0")

air_vars <- c("pm2p5_lag0","pm10_lag0","so2_lag0")

plot_data$variable_fact <- factor(
  plot_data$variable,
  levels = c(met_vars, air_vars),
  labels = c(
    "Mean Temperature (°C)",
    "Maximum Temperature (°C)",
    "Minimum Temperature (°C)",
    "Relative Humidity (%)",
    "Absolute Humidity (g/m³)",
    "Precipitation (mm)",
    "Wind Speed (km/hr)",
    "SPEI-3 Index",
    "PM2.5 (µg/m³)",
    "PM10 (µg/m³)",
    "SO2 (µg/m³)"
  )
)

met_plot <- ggplot(subset(plot_data, variable %in% met_vars),
  aes(x = x, y = rr, color = spec, fill = spec)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_y_log10() +
  facet_wrap(~variable_fact, scales = "free", ncol = 2) +
  labs(y = "Cumulative Relative Risk", x = NULL, tag = "a") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

air_plot <- ggplot(subset(plot_data, variable %in% air_vars),
  aes(x = x, y = rr, color = spec, fill = spec)) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_y_log10() +
  facet_wrap(~variable_fact, scales = "free", ncol = 3) +
  labs(y = "Cumulative Relative Risk", x = NULL, tag = "b") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("figures/supplement/arglag_sensitivity/met_variables_panel.png",
       met_plot, width = 6, height = 8, dpi = 600, bg = "white")

ggsave("figures/supplement/arglag_sensitivity/met_variables_panel.pdf",
       met_plot, device = cairo_pdf,
       width = 6, height = 8, bg = "white")

ggsave("figures/supplement/arglag_sensitivity/air_pollution_panel.png",
       air_plot, width = 9, height = 3, dpi = 600, bg = "white")

ggsave("figures/supplement/arglag_sensitivity/air_pollution_panel.pdf",
       air_plot, device = cairo_pdf,
       width = 9, height = 3, bg = "white")


