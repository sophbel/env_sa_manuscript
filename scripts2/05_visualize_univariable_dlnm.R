
####LOAD DATA & LIBRARIES #####################################################
path_to_package <- "/home/sbelman/Documents/Extra_Projects/IDExtremes/GHRmodel/ghrmodel_0.0.0.9000.tar.gz"

install.packages(path_to_package , 
                 repos = NULL, type = "source", INSTALL_opts = c("--no-multiarch", "--no-test-load"))
library(ghrmodel)
source("/home/sbelman/Documents/BRD/scripts/0_source_functions.R")

### set resolution
time = "weekly"
space = "adm1"
precov = TRUE
## load spatial data
if(space=="adm1"){
  shp<-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_1.shp")
  ## read adjacency matrix
  g <- inla.read.graph(filename = "/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/sa_adjacency_map_adm1.adj")
}
if(space=="adm2"){
  shp<-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_2.shp")
  ## read adjacency matrix
  g <- inla.read.graph(filename = "/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/sa_adjacency_map.adj")
}

# load  data depending on aggregations
if(time=="weekly" & space == "adm1"){
  data <-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly_lag_sc.csv")
  data_unscaled <- fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly_lag.csv")
}
if(time=="weekly" & space == "adm2"){
  data <-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm2_weekly_lag_sc.csv")
  data_unscaled <- fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm2_weekly_lag.csv")
}
if(time=="monthly" & space == "adm1"){
  data <-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_monthly_lag_sc.csv")
  data_unscaled <- fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_monthly_lag.csv")
}
if(time=="monthly" & space == "adm2"){
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
    vaccination_period = ifelse(date < as.Date("2009-01-01"), 1, ifelse(date >= as.Date("2011-01-01"), 3, 2)),
    post_vaccination_2011 = ifelse(date >= as.Date("2011-01-01"), 1, 0),
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
df$vaccination_period <- as.factor(df$vaccination_period)
colnames(df)[grep("present",colnames(df))] <- gsub("_count","", colnames(df)[grep("present",colnames(df))])

### include province as factors for replications
df$id_prov <- as.numeric(factor(df$NAME_1, levels = c("Eastern_Cape", "Free_State", "Gauteng", 
                                                      "KwaZulu-Natal", "Limpopo", "Mpumalanga", 
                                                      "North_West", "Northern_Cape", "Western_Cape")))

## include the proportion of each GPSC per the number sequenced each month
df <- df %>%
  group_by(year, month, NAME_1) %>%
  summarise(month_seq_prov = sum(sequenced), .groups = 'drop') %>%
  left_join(df, by = c("year","month","NAME_1")) 

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
all <- grep("lag0",all, value=TRUE)
cov_names <- grep("tasmax|tasmin|absh|hurs|prlrsum|sfcWind|pm2p5|pm10|o3|so2", all, value=TRUE)
cov_names_labels <- gsub("_lag0", "", cov_names)
### select which GPSCs will be includes
gpsc_vec <- grep("GPSC", all_gpscs, value=TRUE)
gpsc_vec <- grep("count", gpsc_vec, value=TRUE) ## if including the proportions 

################################################################################
############################## READ IN ALL THE FIT RESULTS #####################
################################################################################
#### read in the fits for all the models -- NO INTERACTION 
# st_vec <- c("monthly_adm1","monthly_adm2","weekly_adm1","weekly_adm2")
st_vec <- c("weekly_adm1","weekly_adm2")

stlist <- list()
for(s in 1:length(st_vec)){
  m1 <- fread(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/nointeraction_results_fits_",st_vec[s],"_8_2019.csv"))
  m1$space_time <- st_vec[s]
  stlist[[s]] <- m1
}
fit_list_norm <- rbindlist(stlist)
fit_list_norm$cov <- gsub("_lag0","",fit_list_norm$covariate)

### with the interactions -- WITH INTERACTION
st_vec <- c("weekly_adm1","weekly_adm2")
stlist <- list()

for(s in 1:length(st_vec)){
  m1 <- fread(file=paste0("/home/sbelman/Documents/env_sa_manuscript/models/dlnms/gpsc_results_fits_",st_vec[s],"_allGPSCs_propprov_8_2019.csv"))
  m1$space_time <- st_vec[s]
  stlist[[s]] <- m1
}
interfit_list <- rbindlist(stlist)
interfit_list$cov <- gsub("_lag0","",interfit_list$covariate)

## bind this with the no interaction models
allfits <- rbind(interfit_list, fit_list_norm)


################################################################################
####### CUMULATIVE EFFECTS FOR EACH VARIABLE NO INTERACTION ############
################################################################################
fit_list <- fit_list_norm
fit_list_w1 <- subset(fit_list, fit_list$space_time=="weekly_adm2" & fit_list$lag_num==2)
ggplot(fit_list_w1)+
  geom_hline(yintercept=1, linetype="dashed",color="red")+
  geom_line(aes(x = var, y = exp(cumulative_fit), group=cov))+
  geom_ribbon(aes(x = var, ymin = exp(cum_lowerCI), ymax= exp(cum_upperCI), group=cov), alpha=0.5)+
  theme_bw()+
  ylab("Relative Risk")+
  # ylim(0.1,3)+
  xlab("Variable")+
  facet_wrap(cov~., scales="free")

### both together
# fit_list$cum_upperCI[which(fit_list$cum_upperCI>0.9162907)] <- 0.9162907
  # png(file = paste0("/home/sbelman/Documents/BRD/SouthAfrica/manuscript/figures/supplements/space_time_dlnmeffect.png"), width= 800, height = 500)
ggplot(fit_list)+
  geom_hline(yintercept=1, linetype="dashed",color="red")+
  geom_line(aes(x = var, y = exp(cumulative_fit), group=interaction(cov,space_time), color=space_time))+
  geom_ribbon(aes(x = var, ymin = exp(cum_lowerCI), ymax= exp(cum_upperCI), fill=space_time, group=interaction(cov,space_time)), alpha=0.5)+
  theme_bw()+
  ylab("Relative Risk")+
  theme(axis.text = element_text(size=13), axis.title= element_text(size=13), legend.text = element_text(size=13), legend.title = element_text(size=13), strip.text = element_text(size=13))+
  # ylim(0.5,3)+
  xlab("Variable")+
  facet_wrap(cov~., scales="free")
# dev.off()
################################################################################
####### NO INTERACTION EFFECT PER VARIABLE (spatiotemporal resolutions) ########
################################################################################
### plot all covs and 6 weeks of lags humidity and pm
cov_names <- unique(fit_list$cov)
# cov_names <- c("tasmax","absh","hurs","pm2p5","pm10","so2","o3")
# ax_labs <- c(rep("Degrees Celsius",1), "Absolute Humidity","Relative Humidity (%)", 
#              rep(expression(paste('Concentration (', mu, 'g/m'^3, ')')),4))
plotcov_list <- list()
max_lag = 8
for(c in 1:length(cov_names)){
  # tmp <- subset(fit_list, fit_list$cov%in%cov_names[c] & fit_list$lag_num %in% c(0:8) & fit_list$space_time=="weekly_adm1" )
  # tmp <- subset(fit_list, fit_list$cov%in%cov_names[c] & fit_list$lag_num %in% c(0:6) & fit_list$space_time=="weekly_adm2" )
  tmp <- subset(fit_list, fit_list$cov%in%cov_names[c] & fit_list$lag_num %in% c(0:max_lag) )
  
  tmp$rr <- exp(tmp$fit)
  tmp$lowerCI_rr <- exp(tmp$lowerCI)
  tmp$upperCI_rr <- exp(tmp$upperCI)
  # tmp$lowerCI_rr[which(tmp$lowerCI_rr<0.8)] <- 0.8
  # tmp$upperCI_rr[which(tmp$upperCI_rr>1.1)] <- 1.1
  tmp$lag_week <- paste0("Week",tmp$lag_num)
  tmp$lag_week <- factor(tmp$lag_week, levels = paste0("Week",seq(0,max_lag,1)), labels = paste0("Week",seq(0,max_lag,1)))
  plotcov <- ggplot(tmp)+
    geom_line(aes(x=var,y=rr,group=interaction(lag_week,space_time), color=space_time))+
    geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
    geom_ribbon(aes(x=var, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,space_time), fill=space_time),alpha=0.3)+
    theme_bw()+
    xlab("Var")+
    ylab("Relative Risk")+
    scale_y_continuous(trans="log10")+
    # scale_y_continuous(trans="log10", limits = c(0.9,1.1), breaks = c(0.9, 1, 1.1))+
    ggtitle(cov_names[c])+
    facet_wrap(~lag_week, nrow = 1)+
    # xlab(ax_labs[c])+
    theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13))
  plotcov_list[[c]] <- plotcov
}

(plotcov_list[[1]]/plotcov_list[[2]])
(plotcov_list[[3]]/plotcov_list[[4]])
(plotcov_list[[4]]/plotcov_list[[5]])
(plotcov_list[[6]]/plotcov_list[[7]])
(plotcov_list[[8]]/plotcov_list[[9]])
(plotcov_list[[10]]/plotcov_list[[11]])
(plotcov_list[[12]]/plotcov_list[[13]])


### plot subset at 2 week alg
tmp <- subset(fit_list, fit_list$cov%in%c("hurs","tasmax","prlrsum","pm2p5","pm10","o3","so2","tasmin","tas","sfcWind") & fit_list$lag_num %in% c(2))
# tmp$cov <- factor(tmp$cov, levels = c("hurs","tasmax","pm2p5","pm10"), labels = c("Relative_Humidity", "Max_Temp","PM2.5","PM10"))
ggplot(tmp)+
  geom_line(aes(x=var,y=cumulative_fit,group=interaction(cov,space_time), color=space_time))+
  geom_hline(yintercept = 0, linetype = "dashed", color="red", alpha=0.6)+
  geom_ribbon(aes(x=var, ymin=cum_lowerCI,ymax=cum_upperCI, group=interaction(cov,space_time,lag_num), fill=space_time),alpha=0.3)+
  theme_bw()+
  xlab("Var")+
  ylab("Effect")+
  facet_wrap(lag_num~cov,scales="free", ncol=8)+
  labs(fill="Resolution",color="Resolution") +
  theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13))

tmp <- subset(fit_list, fit_list$cov%in%c("hurs","prlrsum","tasmax","pm2p5","pm10") & fit_list$lag_num %in% c(0:6) & fit_list$space_time=="weekly_adm1" )

### low humidity
tmphigh <- tmp[which(tmp$lag_num%in%c(4)&tmp$lowerCI>0&tmp$cov=="hurs"),] ## what percent humidity is there an increased risk 
tmphigh[which(tmphigh$var<31&tmphigh$var>29),]
exp(tmphigh[which(tmphigh$var<31&tmphigh$var>29),"fit"])
exp(tmphigh[which(tmphigh$var<31&tmphigh$var>29),"lowerCI"])
exp(tmphigh[which(tmphigh$var<31&tmphigh$var>29),"upperCI"])

### high humidity
tmplow <- tmp[which(tmp$lag_num%in%c(6)&tmp$upperCI<0 &tmp$cov=="hurs"),] ## what percent humidity is there an increased risk 
tmplow[which(tmplow$var>85&tmplow$var<88),]
1-exp(tmplow[which(tmplow$var>75&tmplow$var<86),"fit"])
exp(tmplow[which(tmplow$var>75&tmplow$var<86),"lowerCI"])
exp(tmplow[which(tmplow$var>75&tmplow$var<86),"upperCI"])

################################################################################
##### VISUALIZE DATA 
###############################################################################
data2<- fread(file="/home/sbelman/Documents/BRD/SouthAfrica/disease/SA_disease_point_base.csv",quote=FALSE, header = TRUE)
gpsc_df<-data.table(table(data2$GPSC))
gpsc_df <- gpsc_df[which(gpsc_df$N>100),]
gpsc_vec <- as.character(gpsc_df[order(-gpsc_df$N), V1])
gpsc_vec <- paste0("GPSC",gpsc_vec,"_count")
# ## visualize the weekly admin 1 gpsc data
gpsc_spatial_plot <- gpsc_annual_plot <- gpsc_season_plot <- list()
shpnew <- shp
shpnew$NAME_1 <- gsub(" ","_",shpnew$NAME_1)

for(gp in 1:length(gpsc_vec)){
  gpsc_season_plot[[gp]] <- df %>%
    group_by(month)%>%
    summarise(count = sum(!!sym(gpsc_vec[gp])),
              month_tot = sum(week_seq_prov), .groups = 'drop') %>%
    ungroup() %>%
    ggplot()+
    geom_line(aes(x=month,y=count/month_tot))+
    ggtitle(gpsc_vec[gp])+
    theme(legend.position = "none")

  gpsc_annual_plot[[gp]] <- df %>%
    group_by(date)%>%
    summarise(count = sum(!!sym(gpsc_vec[gp])),
              month_tot = sum(week_seq_prov), .groups = 'drop') %>%
    ungroup() %>%
    ggplot()+
    geom_line(aes(x=date,y=count/month_tot))+
    ggtitle(gsub("_count","",gpsc_vec[gp]))+
    theme_bw()+
    ylab("Proportion")+
    xlab("Date")+
    theme(legend.position = "none")

  gpsc_plot_data <- df %>%
    group_by(NAME_1) %>%
    summarise(count = sum(!!sym(gpsc_vec[gp])),
              month_tot = sum(week_seq_prov),
              prop = count/month_tot, .groups = 'drop')

  gpsc_spatial_plot[[gp]] <- shpnew %>%
    left_join(gpsc_plot_data, by = "NAME_1") %>%
    ggplot() +
    geom_sf(aes(fill = prop)) +
    ggtitle(gpsc_vec[gp]) +
    theme(legend.position = "right")
}
# grid.arrange(grobs=gpsc_season_plot)
grid.arrange(grobs=gpsc_annual_plot[c(1:12)])
grid.arrange(grobs=gpsc_spatial_plot)

##### for specific gpsc #####
# df_spatial <- df %>%
#   group_by(province) %>%
#   summarise(pm2p5_avg = mean(pm2p5_lag2, na.rm = TRUE),
#             lineage_prop = mean(GPSC8_count, na.rm = TRUE))  # adjust for GPSC of interest
# 
# cor(df_spatial$pm2p5_avg, df_spatial$lineage_prop, use = "complete.obs")
# ggplot(df_spatial, aes(x = pm2p5_avg, y = lineage_prop)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   labs(x = "Avg PM2.5", y = "Avg GPSC Proportion")


################################################################################
####### CALCULATE THE RELATIVE RISK RATIO FOR EACH GPSC-PM AND HURS ############
################################################################################
data2<- fread(file="/home/sbelman/Documents/BRD/SouthAfrica/disease/SA_disease_point_base.csv",quote=FALSE, header = TRUE)
dtgpsc <- data.table(table(data2$GPSC))[which(data.table(table(data2$GPSC))$N>100 | data.table(table(data2$GPSC))$V1%in%c("8","41"))]
gpsc_vec <- paste0("GPSC",dtgpsc[order(-N)]$V1)
########## calculate the relative risk ratio at 180ug/m3 for pm2.5 and pm10 for each GPSC high vs low and medium vs. low
madm1 <- interfit_list[which(interfit_list$space_time=="weekly_adm2"),]
madm1 <- data.frame(madm1)
# Pivot wider to get high and low side-by-side
madm1_wide <- madm1 %>%
  dplyr::select(var, lag, GPSC, cov,covariate, interaction_level, cumulative_fit, cum_lowerCI, cum_upperCI) %>%
  # dplyr::select(var, lag, GPSC, covariate, interaction_level, fit, lowerCI, upperCI) %>%
  pivot_wider(
    names_from = interaction_level,
    # values_from = c(fit, lowerCI, upperCI),
    values_from = c(cumulative_fit, cum_lowerCI, cum_upperCI),
    names_glue = "{.value}_{interaction_level}"
  )

# Calculate RR, RR-ratio and approximate 95% CI for the ratio
madm1_rr <- madm1_wide %>%
  mutate(
    # # Convert to relative risks
    # RR_high = exp(fit_high),
    # RR_low = exp(fit_low),
    # RR_ratio = exp(fit_high - fit_low),
    # 
    # # Delta method for SE of difference in log(RR)
    # SE_high = (upperCI_high - lowerCI_high) / (2 * 1.96),
    # SE_low  = (upperCI_low - lowerCI_low) / (2 * 1.96),
    # SE_diff = sqrt(SE_high^2 + SE_low^2),
    # 
    # # CI for RR-ratio
    # RR_ratio_lower = exp((fit_high - fit_low) - 1.96 * SE_diff),
    # RR_ratio_upper = exp((fit_high - fit_low) + 1.96 * SE_diff),

    ## cumulative
    RR_high = exp(cumulative_fit_high),
    RR_low = exp(cumulative_fit_low),
    RR_ratio = exp(cumulative_fit_high - cumulative_fit_low),
    # 
    # Delta method for SE of difference in log(RR)
    SE_high = (cum_upperCI_high - cum_lowerCI_high) / (2 * 1.96),
    SE_low  = (cum_upperCI_low - cum_lowerCI_low) / (2 * 1.96),
    SE_diff = sqrt(SE_high^2 + SE_low^2),

    RR_ratio = exp(cumulative_fit_high - cumulative_fit_low),
    RR_ratio_lower = exp((cumulative_fit_high - cumulative_fit_low) - 1.96 * SE_diff),
    RR_ratio_upper = exp((cumulative_fit_high - cumulative_fit_low) + 1.96 * SE_diff)
    
  )
madm1_rr$cov <- gsub("_lag0","",madm1_rr$covariate)
madm1_rr$GPSC <- factor(madm1_rr$GPSC, levels = gpsc_vec)


### calculate the RR at 150ug/m3 of pm2p5 for each
pmtmp <- madm1_rr[which(madm1_rr$cov=="pm2p5"),]
pmplot <- madm1_rr[which(madm1_rr$cov=="pm2p5"&madm1_rr$var%in%c(50)&madm1_rr$lag%in%c("lag2")),]#&madm1_rr$GPSC%notin%c("GPSC56","GPSC8")),]
pmplot[which(pmplot$var==50&pmplot$GPSC=="GPSC14"),"RR_ratio"]
pmplot[which(pmplot$var==50&pmplot$GPSC=="GPSC14"),"RR_ratio_lower"]
pmplot[which(pmplot$var==50&pmplot$GPSC=="GPSC14"),"RR_ratio_upper"]

pmplot[which(pmplot$var==50&pmplot$GPSC=="GPSC21"),"RR_ratio"]
pmplot[which(pmplot$var==50&pmplot$GPSC=="GPSC21"),"RR_ratio_lower"]
pmplot[which(pmplot$var==50&pmplot$GPSC=="GPSC21"),"RR_ratio_upper"]
pm2p5<- ggplot(pmplot)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  geom_pointrange(aes(x=GPSC, y=RR_ratio, ymin=RR_ratio_lower,ymax=RR_ratio_upper,group=interaction(var,lag)))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ylab("Relative Risk Ratio")+
  xlab("GPSC")+
  # ylim(0.85,1.25)+
  ggtitle("PM2.5")+
  facet_grid(lag~var)

pm10plot <- madm1_rr[which(madm1_rr$cov=="pm10"&madm1_rr$var%in%c(100)&madm1_rr$lag%in%c("lag2")),]#&madm1_rr$GPSC%notin%c("GPSC56","GPSC8")),]
pm10plot[which(pm10plot$var==180&pm10plot$GPSC=="GPSC14"),"RR_ratio"]
pm10plot[which(pm10plot$var==180&pm10plot$GPSC=="GPSC14"),"RR_ratio_lower"]
pm10plot[which(pm10plot$var==180&pm10plot$GPSC=="GPSC14"),"RR_ratio_upper"]

pm10plot[which(pm10plot$var==180&pm10plot$GPSC=="GPSC21"),"RR_ratio"]
pm10plot[which(pm10plot$var==180&pm10plot$GPSC=="GPSC21"),"RR_ratio_lower"]
pm10plot[which(pm10plot$var==180&pm10plot$GPSC=="GPSC21"),"RR_ratio_upper"]

pm10 <- ggplot(pm10plot)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  geom_pointrange(aes(x=GPSC, y=RR_ratio, ymin=RR_ratio_lower,ymax=RR_ratio_upper,group=interaction(var,lag)))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ylab("Relative Risk Ratio")+
  ggtitle("PM10")+
  xlab("GPSC")+
  # ylim(0.85,1.25)+
  facet_grid(lag~var)

hursplot <- madm1_rr[which(madm1_rr$cov=="hurs"&madm1_rr$lag%in%c("lag1")),]#&madm1_rr$GPSC%notin%c("GPSC56","GPSC8")),]
hursplot$var1 <- round(hursplot$var,0)
hursplot[which(hursplot$var1==80&hursplot$GPSC=="GPSC17"),"RR_ratio"]
hursplot[which(hursplot$var1==80&hursplot$GPSC=="GPSC17"),"RR_ratio_lower"]
hursplot[which(hursplot$var1==80&hursplot$GPSC=="GPSC17"),"RR_ratio_upper"]


hursplot2 <- hursplot[which(hursplot$var1%in%c(27,80)),]
# hursplot$RR_ratio[which(hursplot$RR_ratio>1.25)] <- 1.25
# hursplot$RR_ratio_lower[which(hursplot$RR_ratio_lower>1.25)] <- 1.25
# hursplot$RR_ratio_upper[which(hursplot$RR_ratio_upper>1.25)] <- 1.25
# hursplot$RR_ratio[which(hursplot$RR_ratio<0.85)] <- 0.85
# hursplot$RR_ratio_lower[which(hursplot$RR_ratio_lower<0.85)] <- 0.85
# hursplot$RR_ratio_upper[which(hursplot$RR_ratio_upper<0.85)] <- 0.85
hursp <- ggplot(hursplot2)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  geom_pointrange(aes(x=GPSC, y=RR_ratio, ymin=RR_ratio_lower,ymax=RR_ratio_upper,group=interaction(var1,lag)))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  ylab("Relative Risk Ratio")+
  ggtitle("Relative Humidity")+
  xlab("GPSC")+
  # ylim(0.85,1.25)+
  facet_grid(lag~var1)

(pm2p5+pm10)/
  (hursp )


################################################################################
####### PLOT GPSC-INTERACTION WITH ENV EFFECTS ############
################################################################################

# tmp2$GPSC <- as.factor(tmp2$GPSC)
gpsc_colors <- c(
  "GPSC1" = "#E41A1C", "GPSC10" = "orange", "GPSC13" = "yellow", 
  "GPSC14" = "#984EA3", "GPSC17" = "darkblue", "GPSC2" = "#A65628", 
  "GPSC21" = "#F781BF", "GPSC26" = "#999999", "GPSC3" = "darkgreen",
  "GPSC5" = "#FC8D62", "GPSC56" = "#8DA0CB", "GPSC8" = "darkgrey", 
  "no_interaction" = "black"
)

##### plot effects for GPSCs together with high and low prevalence################################################
p10_list <- p14_list <- p26_list <- p21_list <- plist <- list()
cov_names_labels <- gsub("_lag0","",cov_names)
cov_names_labels <- unique(interfit_list$cov)
cov_names_labels <- cov_names_labels[2:3]
for(c in 1:length(cov_names_labels)){
  ## plot the effect with medium proportion of each GPSC for HURS, pm2p5, pm10, tasmax
  tmp2 <- subset(allfits, allfits$cov %in%cov_names_labels[c] & allfits$space_time%in%c("weekly_adm2") & allfits$lag_num%in%c(0:8) & allfits$interaction_level %in% c("high","low", "none") & allfits$GPSC%in%c("GPSC10","GPSC14","GPSC21","GPSC26","no_interaction") )

  #### split and pivot so I can look at each GPSC individually####################
  # Step 1: Split into main data and no_interaction
  no_int <- tmp2 %>% filter(GPSC == "no_interaction") %>%
    rename_with(~paste0(., "_noint"), c("fit", "lowerCI", "upperCI"))
  
  main_gpsc <- tmp2 %>% filter(GPSC != "no_interaction")
  
  # Step 2: Join them back by cov, var, and space_time
  allfits_pivoted <- main_gpsc %>%
    left_join(no_int %>% 
                dplyr::select(space_time, cov, var,lag,fit_noint, lowerCI_noint, upperCI_noint),
              by = c("var","cov","space_time","lag"))
  
  p <- ggplot(allfits_pivoted)+
    geom_line(aes(x=var,y=fit, group=interaction(cov,space_time, GPSC, interaction_level),  color=interaction_level))+
    geom_hline(yintercept = 0, linetype = "dashed", color="red", alpha=0.6)+
    geom_line(aes(x = var, y = fit_noint, group=GPSC), color = "black", linetype="dotted") +
    geom_ribbon(aes(x = var, ymax = upperCI_noint, ymin = lowerCI_noint, group=GPSC), fill = "black", alpha=0.2) +
    geom_ribbon(aes(x=var, ymin=lowerCI,ymax=upperCI, group=interaction(cov,space_time, GPSC, interaction_level), fill=interaction_level),alpha=0.2)+
    # geom_ribbon(aes(x=var, ymin=lowerCI_noint,ymax=upperCI_noint, group=interaction(cov,space_time, GPSC)),fill="black",alpha=0.03)+
    scale_color_manual(values = c("high"="red","low"="blue")) +
    scale_fill_manual(values = c("high"="red","low"="blue")) +
    # ylim(-0.2,0.2)+
    ggtitle(cov_names_labels[c])+
    theme_bw()+
    xlab("Var")+
    ylab("Effect")+
    # facet_wrap(GPSC~.,scales="free_x")+
    facet_grid(GPSC~lag,scales="free_x")+
    theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13), strip.text.y = element_text(angle = 0))
  plist[[c]] <- p
  
  ## just for GPSC10
  tmp10 <- subset(allfits, allfits$cov %in%cov_names_labels[c]& allfits$lag_num%in%c(0:8) & allfits$space_time%in%c("weekly_adm2") & allfits$GPSC%in%c("GPSC10","no_interaction") & allfits$interaction_level %in% c("high","low", "none") )
  no_int <- tmp10 %>% filter(GPSC == "no_interaction") %>%
    rename_with(~paste0(., "_noint"), c("fit", "lowerCI", "upperCI"))
  main_gpsc <- tmp10 %>% filter(GPSC != "no_interaction")
  allfits_pivoted10 <- main_gpsc %>%
    left_join(no_int %>% 
                dplyr::select(space_time,lag, cov, var,fit_noint, lowerCI_noint, upperCI_noint),
              by = c("var","lag","cov","space_time"))
  p10 <- ggplot(allfits_pivoted10)+
    geom_line(aes(x=var,y=fit, group=interaction(cov,space_time, GPSC, interaction_level,lag),  color=interaction_level))+
    geom_hline(yintercept = 0, linetype = "dashed", color="red", alpha=0.6)+
    geom_line(aes(x = var, y = fit_noint, group=interaction(lag, interaction_level)), color = "black", linetype="dotted") +
    geom_ribbon(aes(x=var, ymin=lowerCI,ymax=upperCI, group=interaction(cov,space_time, GPSC, interaction_level,lag), fill=interaction_level),alpha=0.2)+
    geom_ribbon(aes(x=var, ymin=lowerCI_noint,ymax=upperCI_noint, group=interaction(cov,space_time, GPSC)),fill="black",alpha=0.03)+
    scale_color_manual(values = c("high"="red","low"="blue")) +
    scale_fill_manual(values = c("high"="red","low"="blue")) +
    ylim(-0.3,0.3)+
    theme_bw()+
    xlab(cov_names_labels[c])+
    ylab("Effect")+
    facet_grid(.~lag,scales="free_x")+
    ggtitle("GPSC10")+
    labs(fill="GPSC Prevalence", color = "GPSC Prevalence")+
    theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13))
  p10_list[[c]] <- p10
  
  ### just for GPSC14 (serotype 23F)
  tmp3 <- subset(allfits, allfits$cov %in%cov_names_labels[c] &allfits$lag_num%in%c(0:8) & allfits$space_time%in%c("weekly_adm2") & allfits$GPSC%in%c("GPSC14","no_interaction") & allfits$interaction_level %in% c("high","low", "none") )
  no_int <- tmp3 %>% filter(GPSC == "no_interaction") %>%
    rename_with(~paste0(., "_noint"), c("fit", "lowerCI", "upperCI"))
  main_gpsc <- tmp3 %>% filter(GPSC != "no_interaction")
  allfits_pivoted3 <- main_gpsc %>%
    left_join(no_int %>% 
                dplyr::select(space_time, lag,cov, lag, var,fit_noint, lowerCI_noint, upperCI_noint),
              by = c("var","lag","cov","space_time"))
  p14 <- ggplot(allfits_pivoted3)+
    geom_line(aes(x=var,y=fit, group=interaction(cov,space_time, GPSC, interaction_level,lag),  color=interaction_level))+
    geom_hline(yintercept = 0, linetype = "dashed", color="red", alpha=0.6)+
    geom_line(aes(x = var, y = fit_noint, group=interaction(lag,interaction_level)), color = "black", linetype="dotted") +
    geom_ribbon(aes(x=var, ymin=lowerCI,ymax=upperCI, group=interaction(cov,space_time, GPSC, interaction_level,lag), fill=interaction_level),alpha=0.2)+
    geom_ribbon(aes(x=var, ymin=lowerCI_noint,ymax=upperCI_noint, group=interaction(cov,space_time, GPSC)),fill="black",alpha=0.03)+
    scale_color_manual(values = c("high"="red","low"="blue")) +
    scale_fill_manual(values = c("high"="red","low"="blue")) +
    # ylim(-0.2,0.2)+
    xlab(cov_names_labels[c])+
    theme_bw()+
    ylab("Effect")+
    facet_grid(.~lag,scales="free_x")+
    ggtitle("GPSC14")+
    labs(fill="GPSC Prevalence", color = "GPSC Prevalence")+
    theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13))
  p14_list[[c]] <- p14
  
  ### just for GPSC21  
  tmp5<- subset(allfits, allfits$cov %in%cov_names_labels[c]& allfits$lag_num%in%c(0:8) & allfits$space_time%in%c("weekly_adm2") & allfits$GPSC%in%c("GPSC21","no_interaction") & allfits$interaction_level %in% c("high","low", "none") )
  no_int <- tmp5 %>% filter(GPSC == "no_interaction") %>%
    rename_with(~paste0(., "_noint"), c("fit", "lowerCI", "upperCI"))
  main_gpsc <- tmp5 %>% filter(GPSC != "no_interaction")
  allfits_pivoted5 <- main_gpsc %>%
    left_join(no_int %>% 
                dplyr::select(space_time,lag, cov, var,fit_noint, lowerCI_noint, upperCI_noint),
              by = c("var","lag","cov","space_time"))
  p21 <- ggplot(allfits_pivoted5)+
    geom_line(aes(x=var,y=fit, group=interaction(cov,space_time, GPSC, interaction_level,lag),  color=interaction_level))+
    geom_hline(yintercept = 0, linetype = "dashed", color="red", alpha=0.6)+
    geom_line(aes(x = var, y = fit_noint, group=interaction(lag, interaction_level)), color = "black", linetype="dotted") +
    geom_ribbon(aes(x=var, ymin=lowerCI,ymax=upperCI, group=interaction(cov,space_time, GPSC, interaction_level,lag), fill=interaction_level),alpha=0.2)+
    geom_ribbon(aes(x=var, ymin=lowerCI_noint,ymax=upperCI_noint, group=interaction(cov,space_time, GPSC)),fill="black",alpha=0.03)+
    scale_color_manual(values = c("high"="red","low"="blue")) +
    scale_fill_manual(values = c("high"="red","low"="blue")) +
    # ylim(-0.2,0.2)+
    theme_bw()+
    xlab(cov_names_labels[c])+
    ylab("Effect")+
    facet_grid(.~lag,scales="free_x")+
    ggtitle("GPSC21")+
    labs(fill="GPSC Prevalence", color = "GPSC Prevalence")+
    theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13))
  p21_list[[c]] <- p21
  
  
  ### just for GPSC26 serotype 19A
  tmp4 <- subset(allfits, allfits$cov %in%cov_names_labels[c]& allfits$lag_num%in%c(0:8) & allfits$space_time%in%c("weekly_adm2") & allfits$GPSC%in%c("GPSC26","no_interaction") & allfits$interaction_level %in% c("high","low", "none") )
  no_int <- tmp4 %>% filter(GPSC == "no_interaction") %>%
    rename_with(~paste0(., "_noint"), c("fit", "lowerCI", "upperCI"))
  main_gpsc <- tmp4 %>% filter(GPSC != "no_interaction")
  allfits_pivoted4 <- main_gpsc %>%
    left_join(no_int %>% 
                dplyr::select(space_time,lag, cov, var,fit_noint, lowerCI_noint, upperCI_noint),
              by = c("var","lag","cov","space_time"))
  p26 <- ggplot(allfits_pivoted4)+
    geom_line(aes(x=var,y=fit, group=interaction(cov,space_time, GPSC, interaction_level,lag),  color=interaction_level))+
    geom_hline(yintercept = 0, linetype = "dashed", color="red", alpha=0.6)+
    geom_line(aes(x = var, y = fit_noint, group=interaction(lag, interaction_level)), color = "black", linetype="dotted") +
    geom_ribbon(aes(x=var, ymin=lowerCI,ymax=upperCI, group=interaction(cov,space_time, GPSC, interaction_level,lag), fill=interaction_level),alpha=0.2)+
    geom_ribbon(aes(x=var, ymin=lowerCI_noint,ymax=upperCI_noint, group=interaction(cov,space_time, GPSC)),fill="black",alpha=0.03)+
    scale_color_manual(values = c("high"="red","low"="blue")) +
    scale_fill_manual(values = c("high"="red","low"="blue")) +
    ylim(-0.3,0.3)+
    theme_bw()+
    xlab(cov_names_labels[c])+
    ylab("Effect")+
    facet_grid(.~lag,scales="free_x")+
    ggtitle("GPSC26")+
    labs(fill="GPSC Prevalence", color = "GPSC Prevalence")+
    theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13))
  p26_list[[c]] <- p26
}

## plot the effects across pm2p5, hurs and pm10
grid.arrange(grobs=plist, nrow=1)
grid.arrange(grobs=p10_list, nrow=1)
grid.arrange(grobs=p14_list, nrow=1)
grid.arrange(grobs=p21_list, nrow=1)
grid.arrange(grobs=p26_list, nrow=1)

p10_list[[2]] + p26_list[[2]] + p21_list[[2]]
plist[[3]]+plist[[4]]+plist[[5]]


##########################################################################################################
####### cumulative for all #####################################################
##########################################################################################################
tmps <- subset(allfits, allfits$cov== "pm10" & allfits$interaction_level %in% c("low","high","none") & allfits$lag_num==3 & allfits$space_time=="weekly_adm1")
# tmps <- subset(allfits, allfits$cov == "pm10" & allfits$interaction_level %in% c("low","none") & allfits$lag_num==3 & allfits$space_time=="weekly_adm2")

no_int <- tmps %>% filter(GPSC == "no_interaction") %>%
  rename_with(~paste0(., "_noint"), c("cumulative_fit", "cum_lowerCI", "cum_upperCI"))

main_gpsc <- tmps %>% filter(GPSC != "no_interaction")

# Step 2: Join them back by cov, var, and space_time
allfits_pivoted <- main_gpsc %>%
  left_join(no_int %>% 
              dplyr::select(space_time, cov, var,lag,cumulative_fit_noint, cum_lowerCI_noint, cum_upperCI_noint),
            by = c("var","cov","space_time","lag"))

p14 <- ggplot(allfits_pivoted)+
  geom_line(aes(x=var,y=exp(cumulative_fit), group=interaction(cov, GPSC, interaction_level,lag),  color=interaction_level))+
  geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
  geom_line(aes(x = var, y = exp(cumulative_fit_noint), group=interaction(lag,interaction_level)), color = "black", linetype="solid") +
  geom_ribbon(aes(x=var, ymin=exp(cum_lowerCI),ymax=exp(cum_upperCI), group=interaction(cov, GPSC, interaction_level,lag), fill=interaction_level),alpha=0.2)+
  geom_ribbon(aes(x=var, ymin=exp(cum_lowerCI_noint) ,ymax= exp(cum_upperCI_noint), group=interaction(cov, GPSC)),fill="black",alpha=0.09)+
  scale_color_manual(values = c("high"="red","low"="blue", "med" = "darkgreen")) +
  scale_fill_manual(values = c("high"="red","low"="blue", "med" = "darkgreen")) +
  ylim(0,3)+
  # xlab(cov_names_labels[c])+
  theme_bw()+
  ylab("Effect")+
  facet_wrap(GPSC~.,scales="free_x")+
  labs(fill="GPSC Prevalence", color = "GPSC Prevalence")+
  theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13),
        strip.text.y = element_text(angle=0))
p14


################################################################################
############ COMPARISON OF GPSC14 AND GPSC3 FOR PRESENTATION
################################################################################
space = "adm2"
### GPSC21
predpm_GPSC21 <- allfits[which(allfits$GPSC%in%c("GPSC21")&allfits$covariate%in%c("pm2p5_lag0")&allfits$space_time=="weekly_adm1"&allfits$interaction_level%in%c("low","high")),]
predpm_GPSC21$interaction_level <- factor(predpm_GPSC21$interaction_level, levels = c("low","high"))
pmgpsc21 <- ggplot(predpm_GPSC21)+
  geom_line( aes(x=predvar, y=exp(fit), group=interaction(interaction_level,cov), color = interaction_level))+
  geom_ribbon( aes(x=predvar,group=interaction(interaction_level,cov), ymin=exp(lowerCI), ymax=exp(upperCI), fill = interaction_level), alpha=0.2)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  scale_y_continuous(trans="log10"
  , breaks = c(0.85,1,1.2), labels = c(0.8,1,1.2), limits = c(0.7,1.5))+
  theme_bw()+
  facet_grid(.~lag_num)+
  xlab("Particulate Matter ug/m3")+
  ylab("Relative Risk")+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14), strip.text.y = element_text(size=14, angle=0),
        legend.position = "right")+
  ggtitle("GPSC21 (99% 19F)")+
  labs(fill="GPSC Prevalence", color = "GPSC Prevalence")

### GPSC1
predpm_GPSC1 <- allfits[which(allfits$GPSC%in%c("GPSC1")&allfits$covariate%in%c("pm2p5_lag0")&allfits$space_time=="weekly_adm1"&allfits$interaction_level%in%c("low","high")),]
predpm_GPSC1$interaction_level <- factor(predpm_GPSC1$interaction_level, levels = c("low","high"))
pmgpsc1 <- ggplot(predpm_GPSC1)+
  geom_line( aes(x=predvar, y=exp(fit), group=interaction(interaction_level,cov), color = interaction_level))+
  geom_ribbon( aes(x=predvar,group=interaction(interaction_level,cov), ymin=exp(lowerCI), ymax=exp(upperCI), fill = interaction_level), alpha=0.2)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  scale_y_continuous(trans="log10"
  , breaks = c(0.85,1,1.2), labels = c(0.8,1,1.2), limits = c(0.7,1.3))+
  theme_bw()+
  facet_grid(.~lag_num)+
  xlab("Particulate Matter ug/m3")+
  ylab("Relative Risk")+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14), strip.text.y = element_text(size=14, angle=0),
        legend.position = "right")+
  ggtitle("GPSC1 (96% 19F)")+
  labs(fill="GPSC Prevalence", color = "GPSC Prevalence")

# 19F types
library(patchwork)
pmgpsc21 /pmgpsc1


### GPSC14 Relative Risk plot for lag 1 or lag2 pm2p5
predpm_none <- allfits[which(allfits$GPSC%in%c("no_interaction")&allfits$cov%in%c("pm2p5","pm10")&allfits$space_time=="weekly_adm1"&allfits$lag_num==2),]
pmgpscnone <- ggplot(predpm_none)+
  geom_line( aes(x=predvar, y=exp(fit), group=(cov)))+
  geom_ribbon( aes(x=predvar,group=(cov), ymin=exp(lowerCI), ymax=exp(upperCI)), alpha=0.2)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  scale_y_continuous(trans="log10", breaks = c(0.8,1,1.2), labels = c(0.8,1,1.2), limits = c(0.7,1.3))+
  theme_bw()+
  facet_grid(cov~.)+
  xlab("Particulate Matter ug/m3")+
  ylab("Relative Risk")+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14), strip.text.y = element_text(size=14, angle=0),legend.text = element_text(size=14))+
  ggtitle("No Interaction")

predpm_GPSC3 <- allfits[which(allfits$GPSC%in%c("GPSC3")&allfits$cov%in%c("pm2p5","pm10")&allfits$space_time=="weekly_adm2"&allfits$lag_num==2&allfits$interaction_level%in%c("low","high")),]
predpm_GPSC3$interaction_level <- factor(predpm_GPSC3$interaction_level, levels = c("low","high"))
pmgpsc3 <- ggplot(predpm_GPSC3)+
  geom_line( aes(x=predvar, y=exp(fit), group=interaction(interaction_level,cov), color = interaction_level))+
  geom_ribbon( aes(x=predvar,group=interaction(interaction_level,cov), ymin=exp(lowerCI), ymax=exp(upperCI), fill = interaction_level), alpha=0.2)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  scale_y_continuous(trans="log10", breaks = c(0.8,1,1.2), labels = c(0.8,1,1.2), limits = c(0.7,1.3))+
  theme_bw()+
  facet_grid(cov~.)+
  xlab("Particulate Matter ug/m3")+
  ylab("Relative Risk")+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14), 
        strip.text.y = element_text(size=14, angle=0),legend.text = element_text(size=14),
        legend.position = "none")+
  ggtitle("GPSC3 - Serotype 8")+
  labs(fill="GPSC Prevalence", color = "GPSC Prevalence")


predpm_GPSC13 <- allfits[which(allfits$GPSC%in%c("GPSC13")&allfits$cov%in%c("pm2p5","pm10")&allfits$space_time=="weekly_adm1"&allfits$lag_num==2&allfits$interaction_level%in%c("low","high")),]
predpm_GPSC13$interaction_level <- factor(predpm_GPSC13$interaction_level, levels = c("low","high"))
predpm_GPSC13[which(predpm_GPSC13$upperCI>log(1.5))] <-log(1.5)
predpm_GPSC13[which(predpm_GPSC13$upperCI<log(0.5))] <-log(0.5)
pmgpsc13 <- ggplot(predpm_GPSC13)+
  geom_line( aes(x=predvar, y=exp(fit), group=interaction(interaction_level,cov), color = interaction_level))+
  geom_ribbon( aes(x=predvar,group=interaction(interaction_level,cov), ymin=exp(lowerCI), ymax=exp(upperCI), fill = interaction_level), alpha=0.2)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  scale_y_continuous(trans="log10", breaks = c(0.8,1,1.2), labels = c(0.8,1,1.2), limits = c(0.7,1.3))+
  theme_bw()+
  facet_grid(cov~.)+
  xlab("Particulate Matter ug/m3")+
  ylab("Relative Risk")+
  # ylim(0.8,1.15)+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14), strip.text.y = element_text(size=14, angle=0),legend.text = element_text(size=14))+
  ggtitle("GPSC13")+
  labs(fill="GPSC Prevalence", color = "GPSC Prevalence")

predpm_GPSC26 <- allfits[which(allfits$GPSC%in%c("GPSC26")&allfits$cov%in%c("pm2p5","pm10")&allfits$space_time=="weekly_adm1"&allfits$lag_num==2&allfits$interaction_level%in%c("low","high")),]
predpm_GPSC26$interaction_level <- factor(predpm_GPSC26$interaction_level, levels = c("low","high"))
predpm_GPSC26[which(predpm_GPSC26$upperCI>log(1.6))] <-log(1.6)
predpm_GPSC26[which(predpm_GPSC26$upperCI<log(0.4))] <-log(0.4)
pmgpsc26 <- ggplot(predpm_GPSC26)+
  geom_line( aes(x=predvar, y=exp(fit), group=interaction(interaction_level,cov), color = interaction_level))+
  geom_ribbon( aes(x=predvar,group=interaction(interaction_level,cov), ymin=exp(lowerCI), ymax=exp(upperCI), fill = interaction_level), alpha=0.2)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  scale_y_continuous(trans="log10", breaks = c(0.8,1,1.2), labels = c(0.8,1,1.2), limits = c(0.7,1.3))+
  theme_bw()+
  facet_grid(cov~.)+
  xlab("Particulate Matter ug/m3")+
  ylab("Relative Risk")+
  # ylim(0.8,1.15)+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14), strip.text.y = element_text(size=14, angle=0),legend.text = element_text(size=14))+
  ggtitle("GPSC26")+
  labs(fill="GPSC Prevalence", color = "GPSC Prevalence")

predpm_GPSC14 <- allfits[which(allfits$GPSC%in%c("GPSC14")&allfits$cov%in%c("pm2p5","pm10")&allfits$space_time=="weekly_adm1"&allfits$lag_num==2&allfits$interaction_level%in%c("low","high")),]
predpm_GPSC14$interaction_level <- factor(predpm_GPSC14$interaction_level, levels = c("low","high"))
predpm_GPSC14[which(predpm_GPSC14$upperCI>log(1.6))] <-log(1.6)
predpm_GPSC14[which(predpm_GPSC14$upperCI<log(0.4))] <-log(0.4)
pmgpsc14 <- ggplot(predpm_GPSC14)+
  geom_line( aes(x=predvar, y=exp(fit), group=interaction(interaction_level,cov), color = interaction_level))+
  geom_ribbon( aes(x=predvar,group=interaction(interaction_level,cov), ymin=exp(lowerCI), ymax=exp(upperCI), fill = interaction_level), alpha=0.2)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  scale_y_continuous(trans="log10", breaks = c(0.8,1,1.2), labels = c(0.8,1,1.2), limits = c(0.7,1.3))+
  theme_bw()+
  facet_grid(cov~.)+
  xlab("Particulate Matter ug/m3")+
  ylab("Relative Risk")+
  # ylim(0.8,1.15)+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14), 
        strip.text.y = element_text(size=14, angle=0),legend.position = "none")+
  ggtitle("GPSC14 - Serotype 23F")+
  labs(fill="GPSC Prevalence", color = "GPSC Prevalence")

predpm_gpsc2 <- allfits[which(allfits$GPSC%in%c("GPSC2")&allfits$cov%in%c("pm2p5","pm10")&allfits$space_time=="weekly_adm1"&allfits$lag_num==2&allfits$interaction_level%in%c("low","high")),]
predpm_gpsc2$interaction_level <- factor(predpm_gpsc2$interaction_level, levels = c("low","high"))
predpm_gpsc2[which(predpm_gpsc2$upperCI>log(1.6))] <-log(1.6)
predpm_gpsc2[which(predpm_gpsc2$upperCI<log(0.4))] <-log(0.4)
pmgpsc2 <- ggplot(predpm_gpsc2)+
  geom_line( aes(x=predvar, y=exp(fit), group=interaction(interaction_level,cov), color = interaction_level))+
  geom_ribbon( aes(x=predvar,group=interaction(interaction_level,cov), ymin=exp(lowerCI), ymax=exp(upperCI), fill = interaction_level), alpha=0.2)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  scale_y_continuous(trans="log10", breaks = c(0.8,1,1.2), labels = c(0.8,1,1.2), limits = c(0.7,1.3))+
  theme_bw()+
  facet_grid(cov~.)+
  xlab("Particulate Matter ug/m3")+
  ylab("Relative Risk")+
  # ylim(0.8,1.15)+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14), 
        strip.text.y = element_text(size=14, angle=0),legend.position = "none")+
  ggtitle("GPSC2")+
  labs(fill="GPSC Prevalence", color = "GPSC Prevalence")

pmgpsc2 + pmgpsc14 + pmgpsc3 + pmgpsc13 +pmgpsc21 +pmgpsc26 + pmgpscnone

### GPSC10 no impact Relative Risk plot for lag 1 or lag2 pm2p5
predpm_GPSC10 <- allfits[which(allfits$GPSC%in%c("GPSC10")&allfits$cov%in%c("pm2p5","pm10")&allfits$space_time=="weekly_adm1"&allfits$lag_num==2 &allfits$interaction_level%in%c("low","high")),]
predpm_GPSC10[which(predpm_GPSC10$upperCI>log(1.6))] <-log(1.6)
predpm_GPSC10[which(predpm_GPSC10$upperCI<log(0.4))] <-log(0.4)

predpm_GPSC10$interaction_level <- factor(predpm_GPSC10$interaction_level, levels = c("low","high"))
pmgpsc10 <- ggplot(predpm_GPSC10)+
  geom_line( aes(x=predvar, y=exp(fit), group=interaction(interaction_level,cov), color = interaction_level))+
  geom_ribbon( aes(x=predvar,group=interaction(interaction_level,cov), ymin=exp(lowerCI), ymax=exp(upperCI), fill = interaction_level), alpha=0.2)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  scale_y_continuous(trans="log10", breaks = c(0.8,1,1.2), labels = c(0.8,1,1.2), limits = c(0.7,1.3))+
  theme_bw()+
  facet_grid(cov~.)+
  xlab("Particulate Matter ug/m3")+
  ylab("Relative Risk")+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14), strip.text.y = element_text(size=14, angle=0),
        legend.position = "none")+
  ggtitle("GPSC10 - Serotype 14")+
  labs(fill="GPSC Prevalence", color = "GPSC Prevalence")
pmgpscnone + pmgpsc10 + pmgpsc14
 ################################################################################

####################### Relative Risk Plots GPSC 14 and GPSC26 with PM2.5 and PM10 #################
pred180 <- allfits[which(allfits$GPSC%in%c("GPSC14")&allfits$cov%in%c("pm2p5","pm10")&allfits$space_time=="weekly_adm2"&allfits$predvar==180&allfits$lag_num<5),]
ggplot(pred180)+
  geom_pointrange(aes(x=lag_num, y=exp(fit), group=cov, ymin=exp(lowerCI), ymax=exp(upperCI), color=interaction_level))+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  facet_grid(cov~interaction_level)+
  theme_bw()+
  xlab("Weekly Lag")+
  ylab("Relative Risk")

pred180_all <- fit_list[which(fit_list$cov%in%c("pm2p5","pm10")&fit_list$space_time=="weekly_adm2"&fit_list$predvar==180),]
ggplot(pred180_all)+
  geom_pointrange(aes(x=lag_num, y=exp(fit), group=cov, ymin=exp(lowerCI), ymax=exp(upperCI)))+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  facet_grid(cov~.)+
  theme_bw()+
  xlab("Weekly Lag")+
  ylab("Relative Risk")

###############################################################################
##### ###VISUALIZE THE DLNM FOR MULTIPLE SEROTYPES AND WITHOUT EFFECT FOR PM2.5 AND PM10
###############################################################################
    ## multiple outcomes
    dlnm_multiout_demog <- fread("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms_suboutcomes/dlnm_subsetoutcomes_results_fits_weekly_adm2_demographic_mixeddf.csv")
    dlnm_multiout_sero <- fread("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms_suboutcomes/dlnm_subsetoutcomes_results_fits_weekly_adm2_serotype_mixeddf.csv")
    dlnm_multiout <- rbind(dlnm_multiout_demog, dlnm_multiout_sero)
    dlnm_multiout$cov <- gsub("_lag0","",dlnm_multiout$covariate)

    max_lag = 12
    ## no interaction single outcome of disease
    dlnm_norm <- fread("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/nointeraction_results_fits_weekly_adm2_12week_mixeddf.csv")
    dlnm_norm$cov <- gsub("_lag0","",dlnm_norm$covariate)
    dlnm_norm_fit <- dlnm_norm[,c("fit","lowerCI","upperCI","cov","predvar","lag_num")]
    colnames(dlnm_norm_fit) <- c("fit_dis","lowerCI_dis","upperCI_dis", "cov","predvar","lag_num")
    
    ## pcv, nvt outcomes
    ## pcv7
    dlnm_pcv7 <- subset(dlnm_multiout, dlnm_multiout$GPSC%in%c("pcv7"))
    dlnm_pcv7 <- dlnm_pcv7[,c("fit","lowerCI","upperCI","cov","predvar","lag_num")]
    colnames(dlnm_pcv7) <- c("fit_pcv7","lowerCI_pcv7","upperCI_pcv7", "cov","predvar","lag_num")
    ##pcv13
    dlnm_pcv13 <- subset(dlnm_multiout, dlnm_multiout$GPSC%in%c("pcv13"))
    dlnm_pcv13 <- dlnm_pcv13[,c("fit","lowerCI","upperCI","cov","predvar","lag_num")]
    colnames(dlnm_pcv13) <- c("fit_pcv13","lowerCI_pcv13","upperCI_pcv13", "cov","predvar","lag_num")
    ## nvt
    dlnm_nvt <- subset(dlnm_multiout, dlnm_multiout$GPSC%in%c("nvt"))
    dlnm_nvt <- dlnm_nvt[,c("fit","lowerCI","upperCI","cov","predvar","lag_num")]
    colnames(dlnm_nvt) <- c("fit_nvt","lowerCI_nvt","upperCI_nvt", "cov","predvar","lag_num")
    dlnm_vt <- left_join(dlnm_pcv7,dlnm_pcv13, by = c("cov","predvar","lag_num"))
    dlnm_vt <- left_join(dlnm_vt,dlnm_nvt, by = c("cov","predvar","lag_num"))
    dlnm_multiout <- left_join(dlnm_multiout,dlnm_vt, by = c("cov","predvar","lag_num"))
    
    
    ## merge no interaction and single outcome here
    dlnm_multiout <- left_join(dlnm_multiout,dlnm_norm_fit, by = c("cov","predvar","lag_num"))
    fits <- subset(dlnm_multiout, dlnm_multiout$cov=="pm2p5" )
    ## clean up
    fits$lag_week <- paste0("Week",fits$lag_num)
    fits$lag_week <- factor(fits$lag_week, levels = paste0("Week",seq(0,max_lag,1)), labels = paste0("Week",seq(0,max_lag,1)))
    fits$rr <- exp(fits$fit)
    fits$lowerCI_rr <- exp(fits$lowerCI)
    fits$upperCI_rr <- exp(fits$upperCI)

    ## clean up for all covs
    fits_covs <- dlnm_multiout
    fits_covs$lag_week <- paste0("Week",fits_covs$lag_num)
    fits_covs$lag_week <- factor(fits_covs$lag_week, levels = paste0("Week",seq(0,max_lag,1)), labels = paste0("Week",seq(0,max_lag,1)))
    fits_covs$rr <- exp(fits_covs$fit)
    fits_covs$lowerCI_rr <- exp(fits_covs$lowerCI)
    fits_covs$upperCI_rr <- exp(fits_covs$upperCI)
    fits_covs$cum_rr <- exp(fits_covs$cumulative_fit)
    fits_covs$cum_lowerCI_rr <- exp(fits_covs$cum_lowerCI)
    fits_covs$cum_upperCI_rr <- exp(fits_covs$cum_upperCI)
    
    ### serotype
    groups <- unique(fits$GPSC)
    outcomes <- grep("Sero",groups, value = TRUE)
    fits_sero <- subset(fits, fits$GPSC%in%outcomes & fits$lag_num%in%seq(0,max_lag,1))
    fits_pcv_type <- subset(fits, fits$GPSC%in%c("pcv7","pcv13","nvt")& fits$lag_num%in%seq(0,max_lag,1))
    fits_distypes <- subset(fits, fits$GPSC%in%c("bact_count","mening_count","female_count")& fits$lag_num%in%seq(0,max_lag,1))
    fits_age <- subset(fits, fits$GPSC%in%c("age_gt65","age_gt80","age_lt5","age_18t64","female_count")& fits$lag_num%in%seq(0,max_lag,1))
    
    ## label group
    pcv13_vec <- paste0("Sero",c( "1", "3", "5", "6A", "7F", "19A"),"_count")
    pcv7_vec <- paste0("Sero",c("4", "6B", "9V", "14", "18C", "19F", "23F"),"_count")
    nvt_vec <- paste0("Sero",c("12F","15A","15BC","16F","8","9N"),"_count")
    fits_sero$pcv_type <- ifelse(fits_sero$GPSC%in%pcv7_vec,"PCV7",
                            ifelse(fits_sero$GPSC%in%pcv13_vec,"PCV13","NVT"))

    fits_sero$pcv_type <- factor(fits_sero$pcv_type, levels = c("PCV7","PCV13","NVT"))
    fits_sero$GPSC <- gsub("_count","",fits_sero$GPSC)
    
    
    ############### PLOT 1: VACCINE GROUP WITH SEROTYPES IN BACKGROUND##########
    ###### plot each vaccine group colored by serotype
    plot_vt7 <- ggplot(fits_sero[which(fits_sero$pcv_type=="PCV7"),])+
      geom_line(aes(x=predvar,y=rr,group=interaction(lag_week,GPSC,pcv_type), color=GPSC),alpha=0.7)+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,GPSC,pcv_type), fill = GPSC),alpha=0.2)+
      geom_line(aes(x=predvar,y=exp(fit_pcv7),group=interaction(lag_week,GPSC,pcv_type)), color="black", alpha=0.7)+
      geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_pcv7),ymax=exp(upperCI_pcv7), group=interaction(lag_week,GPSC,pcv_type)), fill = "black" ,alpha=0.2)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      scale_y_continuous(trans="log10"
                         # , limits = c(0.8,1.1), breaks = c(0.8, 1, 1.1)
                         )+
      ggtitle("PM2.5")+
      facet_grid(GPSC~lag_week)+
      labs(color="Serotype",fill="Serotype")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=8), axis.text.x = element_text(angle = 45,hjust=1, size=13))
    plot_vt13 <- ggplot(fits_sero[which(fits_sero$pcv_type=="PCV13"),])+
      geom_line(aes(x=predvar,y=rr,group=interaction(lag_week,GPSC,pcv_type), color=GPSC),alpha=0.7)+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,GPSC,pcv_type), fill = GPSC),alpha=0.2)+
      geom_line(aes(x=predvar,y=exp(fit_pcv13),group=interaction(lag_week,GPSC,pcv_type)), color="black", alpha=0.7)+
      geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_pcv13),ymax=exp(upperCI_pcv13), group=interaction(lag_week,GPSC,pcv_type)), fill = "black" ,alpha=0.2)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      scale_y_continuous(trans="log10"
                         # ,limits = c(0.8,1.08), breaks = c(0.8, 1, 1.08)
                         )+
      ggtitle("PM2.5")+
      facet_grid(GPSC~lag_week)+
      labs(color="Serotype",fill="Serotype")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=10), axis.text.x = element_text(angle = 45,hjust=1, size=13))
    plot_nvt <- ggplot(fits_sero[which(fits_sero$pcv_type=="NVT"),])+
      geom_line(aes(x=predvar,y=rr,group=interaction(lag_week,GPSC,pcv_type), color=GPSC),alpha=0.7)+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,GPSC,pcv_type), fill = GPSC),alpha=0.3)+
      geom_line(aes(x=predvar,y=exp(fit_nvt),group=interaction(lag_week,GPSC,pcv_type)), color="black", alpha=0.7)+
      geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_nvt),ymax=exp(upperCI_nvt), group=interaction(lag_week,GPSC,pcv_type)), fill = "black" ,alpha=0.3)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      scale_y_continuous(trans="log10"
                         # , limits = c(0.8,1.08), breaks = c(0.8, 1, 1.08)
                         )+
      ggtitle("PM2.5")+
      labs(color="Serotype",fill="Serotype")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=10), axis.text.x = element_text(angle = 45,hjust=1, size=13))
    # 
    plot_vt7
    plot_vt13
    plot_nvt
    
############### PLOT 2: DISEASE WITH VACCINE GROUP IN BACKGROUND##########
    plot_pcvtype <- ggplot(fits_pcv_type)+
      geom_line(aes(x=predvar,y=rr,group=interaction(lag_week,GPSC), color=GPSC))+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,GPSC), fill = GPSC),alpha=0.3)+
      # geom_line(aes(x=predvar,y=exp(fit_dis),group=interaction(lag_week,GPSC)), color="black", alpha=0.7)+
      # geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_dis),ymax=exp(upperCI_dis), group=interaction(lag_week,GPSC)), fill = "black" ,alpha=0.1)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      scale_y_continuous(trans="log10", limits = c(0.8,1.08), breaks = c(0.8, 1, 1.08))+
      ggtitle("PM2.5")+
      facet_wrap(.~lag_week, ncol=7)+
      labs(color="Vaccine Type",fill="Vaccine Type")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
            strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13),
            legend.position = "none")
    
    
    plot_pcvtype_nvt <- ggplot(fits_pcv_type[which(fits_pcv_type$GPSC=="nvt"),])+
      geom_line(aes(x=predvar,y=rr,group=interaction(lag_week,GPSC), color=GPSC), color="pink")+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,GPSC), fill = GPSC),fill="pink",alpha=0.5)+
      # geom_line(aes(x=predvar,y=exp(fit_dis),group=interaction(lag_week,GPSC)), color="black", alpha=0.7)+
      # geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_dis),ymax=exp(upperCI_dis), group=interaction(lag_week,GPSC)), fill = "black" ,alpha=0.1)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      scale_y_continuous(trans="log10", limits = c(0.8,1.2), breaks = c(0.8, 1, 1.2))+
      ggtitle("Non-Vaccine Type PM2.5 Effect")+
      facet_wrap(.~lag_week, ncol=13)+
      labs(color="Vaccine Type",fill="Vaccine Type")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
            strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13),
            legend.position = "none")
    plot_pcvtype_vt <- ggplot(fits_pcv_type[which(fits_pcv_type$GPSC=="pcv7"|fits_pcv_type$GPSC=="pcv13"),])+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,GPSC), fill = GPSC),alpha=0.6)+
      scale_color_manual(values=c("lightgreen","lightblue"), breaks =c("pcv7","pcv13") )+
      scale_fill_manual(values=c("lightgreen","lightblue"), breaks =c("pcv7","pcv13") )+
      geom_line(aes(x=predvar,y=rr,group=interaction(lag_week,GPSC), color=GPSC))+
      
      # geom_line(aes(x=predvar,y=exp(fit_dis),group=interaction(lag_week,GPSC)), color="black", alpha=0.7)+
      # geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_dis),ymax=exp(upperCI_dis), group=interaction(lag_week,GPSC)), fill = "black" ,alpha=0.1)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      scale_y_continuous(trans="log10", limits = c(0.8,1.2), breaks = c(0.8, 1, 1.2))+
      ggtitle("Vaccine Type PM2.5 Effect")+
      facet_wrap(.~lag_week, ncol =13)+
      labs(color="Vaccine Type",fill="Vaccine Type")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
            strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13),
            legend.position = "right")
    plot_pcvtype_nvt/plot_pcvtype_vt
    
    
    ### plot cumulative
    fits_pcv_type <- subset(fits_covs, fits_covs$GPSC%in%c("pcv7","pcv13","nvt")& fits_covs$lag_num%in%seq(0,max_lag,1))
    plot_pcvnvt_cum <- ggplot(fits_pcv_type)+
      geom_line(aes(x=predvar,y=exp(cumulative_fit),group=interaction(GPSC, cov), color=GPSC))+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=exp(cum_lowerCI),ymax=exp(cum_upperCI), group=interaction(GPSC, cov), fill = GPSC),alpha=0.5)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      # scale_y_continuous(trans="log10", limits = c(0.8,2.5), breaks = c(0.8, 1, 2.5))+
      # ggtitle("Non-Vaccine Type PM2.5 Effect")+
      facet_wrap(cov~GPSC, ncol = 3, scales="free_x")+
      labs(color="Vaccine Type",fill="Vaccine Type")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
            strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13),
            legend.position = "none")
    # fits_pcv_type[which(fits_pcv_type$var==180&fits_pcv_type$lag_num%in%c(2)),]
    # fits_pcv_type[which(fits_pcv_type$var==180&fits_pcv_type$lag_num%in%c(3)),]
    # 
    ############### PLOT 3: SPECIFIC SEROTYPES ##########
    plot_23F <- ggplot(fits_sero[which(fits_sero$GPSC=="Sero23F"),])+
      geom_line(aes(x=predvar,y=rr,group=interaction(lag_week,GPSC), color=pcv_type),alpha=0.7)+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,GPSC), fill = pcv_type),alpha=0.1)+
      geom_line(aes(x=predvar,y=exp(fit_dis),group=interaction(lag_week,GPSC)), color="black", alpha=0.7)+
      geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_dis),ymax=exp(upperCI_dis), group=interaction(lag_week,GPSC)), fill = "black" ,alpha=0.1)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      scale_y_continuous(trans="log10", limits = c(0.8,1.08), breaks = c(0.8, 1, 1.08))+
      ggtitle("PM2.5")+
      facet_wrap(GPSC~lag_week, ncol=13)+
      labs(color="Vaccine Type",fill="Vaccine Type")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13))
    
    
    plot_8 <- ggplot(fits_sero[which(fits_sero$GPSC=="Sero8"),])+
      geom_line(aes(x=predvar,y=rr,group=interaction(lag_week,GPSC), color=pcv_type),alpha=0.7)+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,GPSC), fill = pcv_type),alpha=0.1)+
      geom_line(aes(x=predvar,y=exp(fit_dis),group=interaction(lag_week,GPSC)), color="black", alpha=0.7)+
      geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_dis),ymax=exp(upperCI_dis), group=interaction(lag_week,GPSC)), fill = "black" ,alpha=0.1)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      # scale_y_continuous(trans="log10", limits = c(0.8,1.11), breaks = c(0.8, 1, 1.11))+
      scale_y_continuous(trans="log10")+
      ggtitle("Serotype 8 ~ PM2.5")+
      facet_wrap(.~lag_week, ncol=13)+
      labs(color="Vaccine Type",fill="Vaccine Type")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
            strip.text = element_text(size=9), axis.text.x = element_text(angle = 45,hjust=1, size=13),
            legend.position = "right")
    plot_19f <- ggplot(fits_sero[which(fits_sero$GPSC=="Sero19F"),])+
      geom_line(aes(x=predvar,y=rr,group=interaction(lag_week,GPSC), color=pcv_type),alpha=0.7)+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,GPSC), fill = pcv_type),alpha=0.1)+
      geom_line(aes(x=predvar,y=exp(fit_dis),group=interaction(lag_week,GPSC)), color="black", alpha=0.7)+
      geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_dis),ymax=exp(upperCI_dis), group=interaction(lag_week,GPSC)), fill = "black" ,alpha=0.1)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      # scale_y_continuous(trans="log10", limits = c(0.8,1.1), breaks = c(0.8, 1, 1.1))+
      ggtitle("Serotype 19F ~ PM2.5")+
      facet_wrap(.~lag_week, ncol=13)+
      labs(color="Vaccine Type",fill="Vaccine Type")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13))
   
    plot_1 <- ggplot(fits_sero[which(fits_sero$GPSC=="Sero1"),])+
      geom_line(aes(x=predvar,y=rr,group=interaction(lag_week,GPSC), color=pcv_type),alpha=0.7)+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,GPSC), fill = pcv_type),alpha=0.1)+
      geom_line(aes(x=predvar,y=exp(fit_dis),group=interaction(lag_week,GPSC)), color="black", alpha=0.7)+
      geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_dis),ymax=exp(upperCI_dis), group=interaction(lag_week,GPSC)), fill = "black" ,alpha=0.1)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      # scale_y_continuous(trans="log10", limits = c(0.8,1.1), breaks = c(0.8, 1, 1.1))+
      ggtitle("Serotype 1")+
      facet_wrap(.~lag_week, ncol=13)+
      labs(color="Vaccine Type",fill="Vaccine Type")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13))
    
############### PLOT 4: DISEASE disease type and female ##########
    plot_dis <- ggplot(fits_distypes)+
      geom_line(aes(x=predvar,y=rr,group=interaction(lag_week,GPSC), color=GPSC),alpha=0.7)+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,GPSC), fill = GPSC),alpha=0.2)+
      geom_line(aes(x=predvar,y=exp(fit_dis),group=interaction(lag_week,GPSC)), color="black", alpha=0.7)+
      geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_dis),ymax=exp(upperCI_dis), group=interaction(lag_week,GPSC)), fill = "black" ,alpha=0.1)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      # scale_y_continuous(trans="log10", limits = c(0.8,1.08), breaks = c(0.8, 1, 1.08))+
      scale_y_continuous(trans="log10")+
      
      ggtitle("PM2.5")+
      facet_wrap(GPSC~lag_week, ncol=13)+
      labs(color="Vaccine Type",fill="Vaccine Type")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
            strip.text = element_text(size=9), axis.text.x = element_text(angle = 45,hjust=1, size=13),
            legend.position = "none")
    
    
    fits_pcv_type[which(fits_pcv_type$var==180&fits_pcv_type$lag_num%in%c(2)),]
    fits_pcv_type[which(fits_pcv_type$var==180&fits_pcv_type$lag_num%in%c(3)),]
    
    plot_age <- ggplot(fits_age)+
      geom_line(aes(x=predvar,y=rr,group=interaction(lag_week,GPSC,cov), color=GPSC),alpha=0.7)+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,GPSC,cov), fill = GPSC),alpha=0.2)+
      # geom_line(aes(x=predvar,y=exp(fit_dis),group=interaction(lag_week,GPSC)), color="black", alpha=0.7)+
      # geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_dis),ymax=exp(upperCI_dis), group=interaction(lag_week,GPSC)), fill = "black" ,alpha=0.1)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      scale_y_continuous(trans="log10", limits = c(0.9,1.25), breaks = c(0.9, 1, 1.25))+
       # scale_y_continuous(trans="log10")+
      ggtitle("PM2.5")+
      facet_wrap(GPSC~lag_week , ncol=13)+
      labs(color="Vaccine Type",fill="Vaccine Type")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
            strip.text = element_text(size=9), axis.text.x = element_text(angle = 45,hjust=1, size=13), 
            legend.position = "none")
    plot_age
    plot_dis
    
    
    fits_age_cum <- fits_age[which(fits_age$lag_num==3),]
    fits_age_cum$GPSC <- factor(fits_age_cum$GPSC, levels = c("age_lt5","age_18t64","age_gt65", "age_gt80","female_count"))
    plot_age <- ggplot(fits_age_cum)+
      geom_line(aes(x=predvar,y=exp(cumulative_fit),group=interaction(lag_week,GPSC,cov), color=GPSC),alpha=0.7)+
      geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
      geom_ribbon(aes(x=predvar, ymin=exp(cum_lowerCI),ymax=exp(cum_upperCI), group=interaction(lag_week,GPSC,cov), fill = GPSC),alpha=0.2)+
      # geom_line(aes(x=predvar,y=exp(fit_dis),group=interaction(lag_week,GPSC)), color="black", alpha=0.7)+
      # geom_ribbon(aes(x=predvar, ymin=exp(lowerCI_dis),ymax=exp(upperCI_dis), group=interaction(lag_week,GPSC)), fill = "black" ,alpha=0.1)+
      theme_bw()+
      xlab(expression(paste('Concentration (', mu, 'g/m'^3, ')')))+
      ylab("Relative Risk")+
      # scale_y_continuous(trans="log10", limits = c(0.9,1.25), breaks = c(0.9, 1, 1.25))+
      # scale_y_continuous(trans="log10")+
      ggtitle("PM2.5")+
      # facet_wrap(GPSC~lag_week , ncol=13)+
      labs(color="Age Group",fill="Age Group")+
      theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
            strip.text = element_text(size=9), axis.text.x = element_text(angle = 45,hjust=1, size=13), 
            legend.position = "right")
  
    #   ################################################################################
    # ##### SENSITIVITY ANALYSIS ACROSS DEGREES OF FREEDOM
    # ###############################################################################
    # 
    # # ## read in different degrees of freedom
    # # st_vec <- c("weekly_adm1","weekly_adm2")
    # # df31 <- fread(file = paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/nointeraction_results_fits_", st_vec[1], "_12week.csv"))
    # # df31$space_time <- st_vec[1]
    # # df32 <- fread(file = paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/nointeraction_results_fits_", st_vec[2], "_12week.csv"))
    # # df32$space_time <- st_vec[2]
    # # df21 <- fread(file = paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/2degreesfreedom/nointeraction_results_fits_", st_vec[1], "_12week_df2.csv"))
    # # df21$space_time <- st_vec[1]
    # # df22 <- fread(file = paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/2degreesfreedom/nointeraction_results_fits_", st_vec[2], "_12week_df2.csv"))
    # # df22$space_time <- st_vec[2]
    # # dftest3 <- rbind(df31,df32)
    # # dftest3$df <- "3"
    # # dftest2 <- rbind(df21,df22)
    # # dftest2$df <- "2"
    # # dftest <- rbind(dftest2,dftest3)
    # 
    # ## subset
    # dftest <- subset(dftest,dftest$covariate %in% c("pm2p5_lag0","pm10_lag0","hurs_lag0"))
    # # add columns of interest
    # dftest$cov <- gsub("_lag0", "", dftest$covariate)
    # dftest$rr <- exp(dftest$fit)
    # dftest$lowerCI_rr <- exp(dftest$lowerCI)
    # dftest$upperCI_rr <- exp(dftest$upperCI)
    # dftest$lag_week <- paste0("Week",dftest$lag_num)
    # dftest$lag_week <- factor(dftest$lag_week, levels = paste0("Week",seq(0,12,1)), labels = paste0("Week",seq(0,12,1)))
    # 
    # ## plot pm, so2, and hurs across all lags to evaluate
    # ggplot(dftest[which(dftest$cov=="pm10"),])+
    #   geom_line(aes(x=var,y=rr,group=interaction(lag_week,space_time,df,cov), color=df, linetype = space_time))+
    #   geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
    #   geom_ribbon(aes(x=var, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,space_time,df,cov), fill=df, linetype = space_time),alpha=0.1)+
    #   theme_bw()+
    #   xlab(cov_oi)+
    #   ylab("Relative Risk")+
    #   scale_y_continuous(trans="log10")+
    #   # scale_y_continuous(trans="log10", limits = c(0.9,1.1), breaks = c(0.9, 1, 1.1))+
    #   # facet_wrap(.~lag_week, nrow = 1)+
    #   facet_grid(df~lag_week + cov)+
    #   
    #   theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
    #         strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13))
    # 
    # ################################################################################
    # ##### SENSITIVITY ANALYSIS ACROSS LAGS with all having 2 degrees of freedom
    # ###############################################################################
    # st_vec <- c("weekly_adm1","weekly_adm2")
    # lag_vec <- c("6week","8week","12week")
    # # Create all combinations
    # combos <- expand.grid(st = st_vec, lag = lag_vec, stringsAsFactors = FALSE)
    # 
    # # Read and annotate files
    # stlist <- Map(function(st, lag) {
    #   m1 <- fread(file = paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/nointeraction_results_fits_", st, "_", lag, "_mixeddf.csv"))
    #   m1$lag_time <- lag
    #   m1$space_time <- st
    #   return(m1)
    # }, combos$st, combos$lag)
    # 
    # # Combine and finalize
    # fit_list <- rbindlist(stlist)
    # 
    # # add columns of interest
    # fit_list$cov <- gsub("_lag0", "", fit_list$covariate)
    # fit_list$rr <- exp(fit_list$fit)
    # fit_list$lowerCI_rr <- exp(fit_list$lowerCI)
    # fit_list$upperCI_rr <- exp(fit_list$upperCI)
    # fit_list$lag_week <- paste0("Week",fit_list$lag_num)
    # fit_list$lag_week <- factor(fit_list$lag_week, levels = paste0("Week",seq(0,12,1)), labels = paste0("Week",seq(0,12,1)))
    # fit_list$lag_time <- factor(fit_list$lag_time, levels = c("6week","8week","12week"))
    # 
    # ### plot the combos together in a cumulative plot
    # ggplot(fit_list[which(fit_list$space_time=="weekly_adm1"),])+
    #   geom_hline(yintercept=1, linetype="dashed",color="red")+
    #   geom_line(aes(x = var, y = exp(cumulative_fit), group=interaction(cov,space_time, lag_time), color=lag_time))+
    #   geom_ribbon(aes(x = var, ymin = exp(cum_lowerCI), ymax= exp(cum_upperCI), fill=lag_time, group=interaction(cov,space_time,lag_time)), alpha=0.5)+
    #   theme_bw()+
    #   ylab("Relative Risk")+
    #   theme(axis.text = element_text(size=13), axis.title= element_text(size=13), legend.text = element_text(size=13), legend.title = element_text(size=13), strip.text = element_text(size=13))+
    #   # ylim(0.5,3)+
    #   xlab("Variable")+
    #   facet_wrap(cov~., scales="free")
    # 
    # ## plot pm, so2, and hurs across all lags to evaluate
    # cov_int <- c("pm2p5","pm10","hurs")
    # # cov_int <- unique(fit_list$cov)
    # cov_oi <- cov_int[2]
    # fit_var <- fit_list[which(fit_list$cov%in%cov_int),]
    # ggplot(fit_var)+
    #   geom_line(aes(x=var,y=rr,group=interaction(lag_week,space_time,lag_time), color=lag_time, linetype = space_time))+
    #   geom_hline(yintercept = 1, linetype = "dashed", color="red", alpha=0.6)+
    #   geom_ribbon(aes(x=var, ymin=lowerCI_rr,ymax=upperCI_rr, group=interaction(lag_week,space_time,lag_time), fill=lag_time, linetype = space_time),alpha=0.1)+
    #   theme_bw()+
    #   xlab(cov_oi)+
    #   ylab("Relative Risk")+
    #   scale_y_continuous(trans="log10")+
    #   # scale_y_continuous(trans="log10", limits = c(0.9,1.1), breaks = c(0.9, 1, 1.1))+
    #   # facet_wrap(.~lag_week, nrow = 1)+
    #   facet_grid(cov + lag_time~lag_week)+
    #   
    #   theme(axis.text = element_text(size=13),axis.title=element_text(size=13), 
    #         strip.text = element_text(size=13), axis.text.x = element_text(angle = 45,hjust=1, size=13))
    # 

###############################################################################
##### ### EVALUATE RANDOM EFFECTS WHEN THERE IS NO INTERACTION
# ###############################################################################
#         model_out <- readRDS(file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/model_out_summary_list_weekly_adm2_dlnm_none.rds"))
#         base_out <- readRDS(file="/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/base_models/base_model_weekly_20092011_popdens_adm2.rds")
#         shp<-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_2.shp")
#         
#         covs <- sapply(model_out, function(x) x$data$cov)
#         
#         year_n <- nrow(model_out[[1]]$summary.random$id_y)/9
#         ## for yearly replicate over province
#         reg_year <- rep(c("Eastern_Cape", "Free_State", "Gauteng", 
#                           "KwaZulu-Natal", "Limpopo", "Mpumalanga", 
#                           "North_West", "Northern_Cape", "Western_Cape"), each = year_n)
#         
#         ## create outputs and loop through each model
#         idm_all <- idy_all <- NULL
#         for(i in 1:length(covs)){
#           ## monthly
#           idm <- model_out[[i]]$summary.random$id_m
#           idm$cov <- covs[i]
#           idm$base <- base_out$summary.random$id_m$mean
#           
#           ## yearly
#           idy <- model_out[[i]]$summary.random$id_y
#           year_n <- nrow(model_out[[i]]$summary.random$id_y)/9
#        
#           idy$region<-factor(reg_year)
#           idy$cov <- covs[i]
#           idy$base <- base_out$summary.random$id_y$mean
#           
#           idm_all <- rbind(idm_all,idm)
#           idy_all <- rbind(idy_all,idy)
#         }
#         # basem <- base_out$summary.random$id_m
#         # basem$cov <- "base"
#         # idm_all <- rbind(idm_all[,c(1:8,10)],basem)
#         
#         # basey <- base_out$summary.random$id_y
#         # basey$region<-factor(reg_year)
#         # basey$cov <- "base"
#         # idy_all <- rbind(idy_all[,c(1:8,10:11)],basey)
#         
#         colnames(idm_all) <- c("month", "mean","sd","lowerCI","median","upperCI","mode","kld","cov","base")
#         colnames(idy_all) <- c("year", "mean","sd","lowerCI","median","upperCI","mode","kld","region","cov","base")
#         idm_all$cov <- gsub("_lag0","",idm_all$cov)
#         idm_all$month <- factor(idm_all$month, levels = seq(1,12,1), labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
#         #### plot seasonal with each variable
#         ggplot(idm_all[which(idm_all$cov%in%c("pm2p5","pm10","tasmax","hurs","prlrsum")),])+
#           geom_hline(yintercept=0,linetype="dashed",color="red",alpha=0.8)+
#           geom_line(aes(x=month,y = median, group=cov,color=cov))+
#           geom_ribbon(aes(x=month,ymin = lowerCI,ymax=upperCI, group=cov,fill=cov),alpha=0.2)+
#           geom_line(aes(x=month,y=base, group=cov),color="black")+
#           theme_bw()+
#           theme(axis.text.x = element_text(angle=45,hjust=1, size=13),axis.text = element_text(size=13), axis.title=element_text(size=13),legend.text = element_text(size=13))+
#           xlab("Month")+
#           ylab("Seasonal Effect")+
#           labs(fill="", color="")#+
#           facet_wrap(.~cov, ncol=1)
#         
#         #### plot seasonal with each variable
#         ggplot(idm_all[which(idm_all$cov%in%c("pm2p5","pm10","tasmax","hurs","prlrsum")),])+
#           geom_hline(yintercept=0,linetype="dashed",color="red",alpha=0.8)+
#           geom_line(aes(x=month,y = median, group=cov,color=cov))+
#           geom_ribbon(aes(x=month,ymin = lowerCI,ymax=upperCI, group=cov,fill=cov),alpha=0.2)+
#           geom_line(data=idm_all[which(idm_all$cov=="base"),],aes(x=month,y=median, group=cov),color="black")+
#           geom_ribbon(data=idm_all[which(idm_all$cov=="base"),],aes(x=month,ymin = lowerCI,ymax=upperCI, group=cov),fill="black",alpha=0.2)+
#           theme_bw()+
#           theme(axis.text.x = element_text(angle=45,hjust=1, size=13),axis.text = element_text(size=13), axis.title=element_text(size=13),legend.text = element_text(size=13))+
#           xlab("Month")+
#           ylab("Seasonal Effect")+
#           labs(fill="", color="")+
#           facet_wrap(cov~.)
#         
#         #### plot seasonal with each variable
#         ggplot(idy_all)+
#           geom_line(aes(x=year,y = median, group=cov,color=cov,group=region))+
#           geom_line(aes(x=year,y=base,group=region),color="black")+
#           theme_bw()+
#           facet_grid(region~cov)
# 
#         
#         ## plot spatial effect difference with each
#         ## blue is a reduction in spatial effect attributed to the covariate
#         ## excluding wind for talks 
#         plot_spatial_effects(model_out[[6]],shp, nrow(shp), structure=FALSE)
#         # nums <- c(2:9)
#         space_list <-list()
#         for(k in 1:length(covs)){
#           # for(k in 1:length(nums)){
#           space_list[[k]] <- plot_spatial_effects_diff(model_out[[k]],base_out,shp,nrow(shp),"Model","Base",TRUE)
#         }
#         grid.arrange(grobs = space_list, nrow = 2, bottom = "Blue areas indicate where the covariate explained spatial variation in pneumococcal disease.")
#         
#         
# ################# spatial effect difference with GPSC14 without ################
#         space<-"adm1"
#         model_outnone <- readRDS(file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/model_out_summary_list_",time,"_",space,"_dlnm_none.rds"))
#         model_outnonepm2p5 <- model_outnone[[6]]
#         sapply(model_outnone, function(x) x$data$cov)
#         
#         
#         model_out <- readRDS(file=paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/model_out_summary_list_",time,"_",space,"_dlnm_GPSC14_count.rds"))
#         sapply(model_out, function(x) x$data$cov)
#         model_outpm2p5 <- model_out[[4]]
# 
#         model_outpm2p5$summary.random$id_u
#         # plot_spatial_effects(model_outpm2p5,shp,9)
#         plot_spatial_effects_diff(model_outpm2p5,model_outnonepm2p5, shp, 9, title_a = "GPSC14", title_b = "none", structured = TRUE)

    
    ################################################################################
    ####### MAE FOR EACH ACROSS LOCATION ############
    ##### calculate the mae relative to the actual cases
    ################################################################################
    # shp<-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_1.shp")
    # cov_names <- c("hurs","pm2p5","pm10")
    # plist_covs <- list()
    # for(covs in 1:length(cov_names)){
    #   mod <- readRDS(paste0("/home/sbelman/Documents/BRD/SouthAfrica/models/outputs/dlnms/models/dlnm_model_univariable_",cov_names[covs],"_adm1_weekly.rds"))
    #   ### fitted values
    #   summary_fitted_values <- mod$summary.fitted.values
    #   
    #   ### true data
    #   df_mae <- cbind(df[,c("disease","province","year")], summary_fitted_values$mean)
    #   colnames(df_mae) <- c("disease","NAME_1","year","fit_mean")
    #   df_mae$fit_mean <- as.numeric(df_mae$fit_mean)
    #   df_mae$disease <- as.numeric(df_mae$disease)
    #   
    #   ### group by district and calculate the per district mae
    #   df_mae2 <- df_mae %>% group_by(NAME_1) %>%
    #     dplyr::summarise(mae = mean(abs(disease-fit_mean), na.rm = TRUE),
    #                      mape = mean(abs((disease - fit_mean) / disease), na.rm = TRUE),
    #                      rmse = sqrt(mean((disease - fit_mean)^2, na.rm = TRUE)),
    #                      mean_disease = mean(disease, na.rm = TRUE)) %>%
    #     mutate(mape = mape * 100) %>%
    #     mutate(norm_rmse = rmse / mean_disease)
    #   shp2 <- left_join(shp, df_mae2, by="NAME_1")
    #   plist_covs[[covs]] <- ggplot(shp2)+
    #     geom_sf(aes(fill=rmse)) + 
    #     scale_fill_viridis_c(option = "plasma", name = "norm_rmse") +
    #     ggtitle(cov_names[covs])+
    #     theme_bw()
    # }
    # grid.arrange(grobs = plist_covs, nrow = 1)
    # 