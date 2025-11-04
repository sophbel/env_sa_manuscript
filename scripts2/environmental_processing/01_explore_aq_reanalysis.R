library(data.table)
library(dplyr)
library(sf)  
library(tidyr)
library(lubridate)
library(ggplot2)


################################################################################
## COMPARE MERRA AND CAMS DATA SETS
################################################################################
aq_adm2_weekly_cams<-fread("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/reanalysis/aq_adm2_weekly_cams.csv")
aq_adm2_weekly_merra<-fread("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/reanalysis/aq_adm2_weekly_merra.csv")
aq_adm2_weekly_merraadj<-fread("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/reanalysis/aq_adm2_weekly_merra_adjusted.csv")
aq_adm2_weekly_ml1km<-fread("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/reanalysis/aq_adm2_weekly_ml1km.csv")

### compare across regions of South Africa using shp file
shp<-st_read("/home/sbelman/Documents/BRD/SouthAfrica/shps/gadm41_ZAF_2.shp")
shp1<-st_read("/home/sbelman/Documents/BRD/SouthAfrica/shps/gadm41_ZAF_1.shp")

## turn the shp file into an index for province
prov_idx <- shp[,c("GID_2","NAME_2","GID_1","NAME_1")]
prov_idx <- st_drop_geometry(prov_idx)

## merge datasets 
# colnames(aq_adm2_weekly_cams) <- c("epiweek", "adm2_id", "adm2_name", "pm2p5_cams", "pm10_cams", "o3_cams", "so2_cams")
# colnames(aq_adm2_weekly_merra) <- c("epiweek", "adm2_id", "adm2_name", "pm2p5_merra", "pm10_merra", "o3_merra", "so2_merra")
aq_adm2_weekly_cams$dataset <- "cams"
aq_adm2_weekly_merra$dataset <- "merra"

## include the same columns
aq_adm2_weekly_merraadj$pm10 <- NA
aq_adm2_weekly_merraadj$o3 <- NA
aq_adm2_weekly_merraadj$so2 <- NA
aq_adm2_weekly_merraadj$dataset <- "merra_adj"

aq_adm2_weekly_ml1km
aq_adm2_weekly_ml1km$pm10 <- NA
aq_adm2_weekly_ml1km$o3 <- NA
aq_adm2_weekly_ml1km$so2 <- NA
aq_adm2_weekly_ml1km$dataset <- "ml1km"

# ## convert to same unit set
## pm2p5 from kg m -3 to ug/m3
aq_adm2_weekly_merraadj$pm2p5 <- aq_adm2_weekly_merraadj$pm2p5 * 1e9
aq_adm2_weekly_merra$pm2p5 <- aq_adm2_weekly_merra$pm2p5 * 1e9
aq_adm2_weekly_cams$pm2p5 <- aq_adm2_weekly_cams$pm2p5 * 1e9
## NOTE: the machine learning method is already in ug/m3 can be checked here (/esarchive/recon/umd/global-ground-level/daily/pm2p5/)

## pm10 for cams from kg m-3  to ug/m3
aq_adm2_weekly_cams$pm10 <- aq_adm2_weekly_cams$pm10 * 1e9
## pm10 for cams from kg m-1  to ug/m3
air_density_reg <- 1.225
air_density_so2 <- 2.63  # SO2 density in kg/m³ (standard at sea level, 15°C)
air_density_o3 <- 2.14  # O3 density in kg/m³ (standard at sea level, 15°C)

# aq_adm2_weekly_merra$pm10 <- aq_adm2_weekly_merra$pm10 * air_density_reg * 1e9
aq_adm2_weekly_merra$pm10 <- aq_adm2_weekly_merra$pm10  * 1e9 ### wrong units

## o3 and so2 for all
aq_adm2_weekly_merra$o3 <- aq_adm2_weekly_merra$o3 * air_density_o3 * 1e9
aq_adm2_weekly_merra$so2 <- aq_adm2_weekly_merra$so2 * air_density_so2 * 1e9
aq_adm2_weekly_cams$o3 <- aq_adm2_weekly_cams$o3 * air_density_o3 * 1e9
aq_adm2_weekly_cams$so2 <- aq_adm2_weekly_cams$so2 * air_density_so2 * 1e9


## bind datasets
aq_adm2_weekly <- rbind(aq_adm2_weekly_cams,aq_adm2_weekly_merra, aq_adm2_weekly_merraadj, aq_adm2_weekly_ml1km)

colnames(prov_idx) <- c("adm2_id",colnames(prov_idx)[2:4])
aq_adm2_weekly <- left_join(aq_adm2_weekly,prov_idx, by="adm2_id")
# aq_adm2_weekly <- left_join(aq_adm2_weekly_cams,aq_adm2_weekly_merra, by= c("epiweek","adm2_id","adm2_name", "dataset"))

## split the epiweek column to create a date column with the week
aq_adm2_weekly[,c("week_num","year") := tstrsplit(epiweek, "-", fixed = TRUE)]
aq_adm2_weekly$week_num <- as.numeric(aq_adm2_weekly$week_num)
aq_adm2_weekly$year <- as.numeric(aq_adm2_weekly$year)
aq_adm2_weekly[, week := floor_date(make_date(year) + weeks(week_num - 1), "week")]
# # Check for duplicates
duplicates <- aq_adm2_weekly %>%
  group_by(adm2_id, week, dataset) %>%
  filter(n() > 1)
### take the means of the duplicate weeks
aq_adm2_weekly <- aq_adm2_weekly %>%
  group_by(week, adm2_id,dataset) %>%
  mutate(
    pm2p5 = mean(pm2p5, na.rm = TRUE),
    pm10 = mean(pm10, na.rm = TRUE),
    o3 = mean(o3, na.rm = TRUE),
    so2 = mean(so2, na.rm = TRUE),
  ) %>%
  ungroup()

aq_adm2_weekly <- subset(aq_adm2_weekly, aq_adm2_weekly$year<2021 )

#### plot comparison across dates
# tmp <- subset(aq_adm2_weekly, aq_adm2_weekly$adm2_name == "Mopani")
tmp <- aq_adm2_weekly
# Assuming your tibble is named `tmp`
tmp_long <- tmp %>%
  pivot_longer(
    cols = c(pm2p5, pm10, o3, so2), # Columns to pivot
    names_to = "variable",          # Name of the new column for variable names
    values_to = "value"             # Name of the new column for values
  )


##############################################################################
#### monthly data sets
##############################################################################
### merge to monthly
### aggregate them monthly
aq_adm2_monthly<-aq_adm2_weekly %>%
  mutate(year_month = floor_date(week, unit = "month")) %>%  # Convert date to start of the month 
  group_by(adm2_id, adm2_name, dataset, year_month, GID_1, NAME_1) %>%  # Group by GID_2, district, and weekly date
  summarise(
    pm2p5 = mean(pm2p5, na.rm = TRUE),  
    pm10 = mean(pm10, na.rm = TRUE), 
    o3 = mean(o3, na.rm = TRUE), 
    so2 = mean(so2, na.rm = TRUE)
  ) %>%
  ungroup()        

aq_adm1_monthly <- aq_adm2_monthly %>%
  group_by(year_month, GID_1, NAME_1, dataset) %>%
  summarize(
    pm2p5 = mean(pm2p5, na.rm = TRUE),
    pm10 = mean(pm10, na.rm = TRUE),
    o3 = mean(o3, na.rm = TRUE),
    so2 = mean(so2, na.rm = TRUE),
  ) %>%
  ungroup()
#### plot comparison across dates
aq_adm1_monthly$year <- year(aq_adm1_monthly$year_month)
tmp<- subset(aq_adm1_monthly, aq_adm1_monthly$year < 2020)
# Assuming your tibble is named `tmp`
tmp_long <- tmp %>%
  pivot_longer(
    cols = c(pm2p5, pm10, o3, so2), # Columns to pivot
    names_to = "variable",          # Name of the new column for variable names
    values_to = "value"             # Name of the new column for values
  )

##############################################################
### COMPARE MERRA AND CAMS WITH AN OBSERVATION STATION
##############################################################
files <- list.files("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/observations/")
files<-files[-grep("pull",files)]
obs_meta <- matrix(nrow=length(files),ncol=4)
# Split by commas and remove extra spaces
split_files <- lapply(files, function(x) strsplit(x, ",\\s*")[[1]])
# Bind into a table
result_table <- do.call(rbind, lapply(split_files, function(x) {
  # Ensure all rows have the same length
  length(x) <- max(lengths(split_files))
  return(x)
}))
# manually edit this table for better names
for(k in 1){ result_table[1,3] <- "Free State"
result_table[2,3] <- "Mpumalanga"
result_table[3,2] <- "cape town"
result_table[3,3] <- "Western Cape"
result_table[4,3] <- "KwaZulu Natal"
result_table[5,3] <- "Mpumalanga"
result_table[6,3] <- "Gauteng"
result_table[7,3] <- "Mpumalanga"
result_table[8,3] <- "Mpumalanga"
result_table[9,3] <- "KwaZulu Natal"
result_table[10,3] <- "KwaZulu Natal"
result_table[11,3] <- "Gauteng"
result_table[12,3] <- "Gauteng"
result_table[13,3] <- "Gauteng"
result_table[13,2] <- "centurion"
result_table[14,3] <- "KwaZulu Natal"
result_table[15,3] <- "KwaZulu Natal"
result_table[16,3] <- "Mpumalanga"
result_table[17,3] <- "Limpopo"
result_table[18,3] <- "Western Cape"
result_table[19,3] <- "KwaZulu Natal"
result_table[20,3] <- "Western Cape"
}
## call data table and rebind names
result_table <- data.table(result_table)
colnames(result_table) <- c("Municipality","District","Province")
result_table$file_names <- files

### save result table list of locations so that it can be read in later
# saveRDS(result_table, file = "/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/observation_file_list.rds")

############# GERT SIBANDE MPUMALANGA ##############
## assign names to provinces
aq_gertsibande_monthly <- aq_adm2_monthly[grep("gert", aq_adm2_monthly$adm2_name, ignore.case = TRUE),]
## read in gert sibande observational data
GS_files <- grep("gert", files, value=TRUE)
sites <- sub(",.*", "", GS_files)
obs_gs_list<-list()
for(s in 1:length(sites)){
  obs<- fread(paste0("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/observations/",GS_files[s]))
  obs$site <- sites[s]
  obs<-obs[,c("date","pm25","pm10","o3","so2","site")]
  obs_gs_list[[s]]<-obs
 }
## bind into one data table
gs_observations <- rbindlist(obs_gs_list)
gs_observations[, date := as.Date(date, format = "%Y/%m/%d")]

### take the daily mean 
gs_observations <- gs_observations %>%
  mutate(year_month = floor_date(date, "month")) %>% # Extract year-month
  group_by(year_month) %>%
  summarise(
    pm2p5 = mean(pm25, na.rm = TRUE),
    pm10 = mean(pm10, na.rm = TRUE),
    o3 = mean(o3, na.rm = TRUE),
    so2 = mean(so2, na.rm = TRUE)
  )
gs_observations$dataset <- "saaqis"
## add province IDs
ids<-unique(aq_gertsibande_monthly[,c("adm2_id","adm2_name","GID_1","NAME_1")])
gs_observations$adm2_id <- ids$adm2_id
gs_observations$adm2_name <- ids$adm2_name
gs_observations$GID_1 <- ids$GID_1
gs_observations$NAME_1 <- ids$NAME_1

## reorder to match
cols <- colnames(aq_gertsibande_monthly)
gs_observations<-gs_observations[,cols]
gs_observations$pm2p5[which(gs_observations$pm2p5>300)] <- 300


## bind observations with reanalysis
aq_gertsibande_monthly<-rbind(aq_gertsibande_monthly,gs_observations)
aq_gertsibande_monthly$location <- "Gert_Sibande"
aq_gertsibande_monthly$adm2_name <- gsub(" ","_", aq_gertsibande_monthly$adm2_name)

# write.table(aq_gertsibande_monthly, file= "/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/aq_gertsibande_monthly_obsreanalysis.csv", quote=FALSE, sep="," ,col.names = TRUE,row.names = FALSE)


# ## plot all together 
tmp_long <- aq_gertsibande_monthly %>%
  pivot_longer(
    cols = c(pm2p5), # Columns to pivot
    # cols = c(pm2p5, pm10, o3, so2), # Columns to pivot
    names_to = "variable",          # Name of the new column for variable names
    values_to = "value"             # Name of the new column for values
  )

# ## subset by the minimum from the observations
obs_min <- min(gs_observations$year_month)
# png("/home/sbelman/Documents/BRD/SouthAfrica/airquality/plots/gert_sibande_pm2p5_comparison.png", width=600,height=300)
tmp_long<- subset(tmp_long, tmp_long$year_month>=obs_min & year(tmp_long$year_month)<2021)
mpm_pm2 <- ggplot(tmp_long)+
  geom_point(aes(x=year_month, y=value, group=variable,   group= NAME_1, color=dataset, shape=dataset),alpha=0.7, size=3)+
  theme_bw()+
  theme(strip.text.y = element_text(angle=0))+
  scale_color_manual(
    name = "Dataset",  # Legend title
    values = c("cams" = "blue", "merra" = "red", "merra_adj" = "darkgreen", "ml1km" = "#AA336A","saaqis" = "black"),
    labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected", "Gapless_1km_ML","SAAQIS (Observations)"))+
  scale_shape_manual(
    name = "Dataset",  # Legend title
    values = c("cams" = 16, "merra" = 17, "merra_adj" = 15, "ml1km" = 3, "saaqis" = 4),
    labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected", "Gapless_1km_ML","SAAQIS (Observations)"))+
  xlab("Year (monthly data)") +
  ylab("Value (ug/m3)") +
  # ylim(0,500)+
  ggtitle("Gert Sibande, Mpumalanga")+
  facet_grid(variable ~ ., scales="free")+
  theme(legend.position = "none")
# dev.off()

############# GAUTENG ##############
## monthly
    ## assign names to provinces
    aq_gaut_monthly <- aq_adm2_monthly[grep("Gauteng", aq_adm2_monthly$NAME_1, ignore.case = TRUE),]
    aq_gaut_monthly<-aq_gaut_monthly %>%
    group_by(NAME_1, GID_1, dataset, year_month) %>%
      summarise(
        pm2p5 = mean(pm2p5, na.rm = TRUE),
        pm10 = mean(pm10, na.rm = TRUE),
        o3 = mean(o3, na.rm = TRUE),
        so2 = mean(so2, na.rm = TRUE))
        
    ## pull gauteng obs
    gaut_files <- result_table[which(result_table$Province=="Gauteng"),]$file_names
    gaut_mat<- result_table[which(result_table$Province=="Gauteng"),]
    obs_gaut_list<-list()
    for(s in 1:length(gaut_files)){
      obs<- fread(paste0("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/observations/",gaut_files[s]))
      obs$site <- gaut_mat[s,"Municipality"]
      if("pm25"%notin%colnames(obs)){
        obs$pm25 <- NA
      }
      obs<-obs[,c("date","pm25","pm10","o3","so2","site")]
      obs_gaut_list[[s]]<-obs
    }
    ## bind into one data table
    gaut_observations <- rbindlist(obs_gaut_list)
    gaut_observations[, date := as.Date(date, format = "%Y/%m/%d")]
    
    ### take the daily mean 
    gaut_observations <- gaut_observations %>%
      mutate(year_month = floor_date(date, "month")) %>% # Extract year-month
      group_by(year_month) %>%
      summarise(
        pm2p5 = mean(pm25, na.rm = TRUE),
        pm10 = mean(pm10, na.rm = TRUE),
        o3 = mean(o3, na.rm = TRUE),
        so2 = mean(so2, na.rm = TRUE)
      )
    gaut_observations$dataset <- "saaqis"
    ## add province IDs
    ids<-unique(aq_gaut_monthly[,c("GID_1","NAME_1")])
    gaut_observations$GID_1 <- ids$GID_1
    gaut_observations$NAME_1 <- ids$NAME_1
    
    ## reorder to match
    aq_gaut_monthly<-aq_gaut_monthly[,3:ncol(aq_gaut_monthly)]
    cols <- colnames(aq_gaut_monthly)
    gaut_observations<-gaut_observations[,cols]
    gaut_observations$pm2p5[which(gaut_observations$pm2p5>300)] <- 300
      
    ## bind observations with reanalysis
    aq_gaut_monthly<-rbind(aq_gaut_monthly,gaut_observations)
    aq_gaut_monthly$location <- "Gauteng"
    # write.table(aq_gaut_monthly , file= "/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/aq_gaut_monthly_obsreanalysis.csv", quote=FALSE, col.names = TRUE,row.names = FALSE)
    
    
    ########## plots #############
    ## plot all together 
    tmp_long <- aq_gaut_monthly %>%
      pivot_longer(
        cols = c(pm2p5, pm10, o3, so2), # Columns to pivot
        names_to = "variable",          # Name of the new column for variable names
        values_to = "value"             # Name of the new column for values
      )

    ## subset by the minimum from the observations
    obs_min <- min(gaut_observations$year_month)
    tmp_long<- subset(tmp_long, tmp_long$year_month>=obs_min & year(tmp_long$year_month)<2021)
    gaut_pm2 <- ggplot(tmp_long[which(tmp_long$variable=="pm2p5"),])+
      geom_point(aes(x=year_month, y=value, group=variable, color=dataset, shape = dataset),size=3,alpha=0.7)+
      theme_bw()+
      theme(strip.text.y = element_text(angle=0))+
      scale_color_manual(
        name = "Dataset",  # Legend title
        values = c("cams" = "blue", "merra" = "red", "merra_adj" = "darkgreen", "ml1km" = "#AA336A","saaqis" = "black"),
        labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected",  "Gapless_1km_ML","SAAQIS (Observations)"))+
      scale_shape_manual(
        name = "Dataset",  # Legend title
        values = c("cams" = 16, "merra" = 17, "merra_adj" = 15, "ml1km" = 3, "saaqis" = 4),
        labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected",  "Gapless_1km_ML","SAAQIS (Observations)"))+
      xlab("Year (monthly data)") +
      ylab("Value (ug/m3)") +
      ggtitle("Gauteng")+
      ylim(0,300)+
      facet_grid(variable ~ ., scales="free")+
      theme(legend.position = "none")
    # 
    # png("/home/sbelman/Documents/BRD/SouthAfrica/airquality/plots/gauteng_pm2p5_comparison.png", width=600,height=300)
    # ggplot(subset(tmp_long,tmp_long$variable=="pm2p5"))+
    #   geom_point(aes(x=year_month, y=value, group=variable, shape=dataset, color=dataset),alpha=0.7,size=4)+
    #   theme_bw()+
    #   theme(strip.text.y = element_text(angle=0), axis.title = element_text(size=12), axis.text = element_text(size=12))+
    #   scale_color_manual(
    #     name = "Dataset",  # Legend title
    #     values = c("cams" = "blue", "merra" = "red", "merra_adj" = "darkgreen", "ml1km" = "#AA336A","saaqis" = "black"),
    #     labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected", "Gapless_1km_ML","SAAQIS (Observations)"))+
    #   scale_shape_manual(
    #     name = "Dataset",  # Legend title
    #     values = c("cams" = 16, "merra" = 17, "merra_adj" = 15, "ml1km" = 3, "saaqis" = 4),
    #     labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected", "Gapless_1km_ML","SAAQIS (Observations)"))+
    #   xlab("Year (monthly data)") +
    #   ylab("Value (ug/m3)") +
    #   ggtitle("Gauteng")+
    #   ylim(0,300)+
    #   facet_grid(variable ~ ., scales="free")
    # dev.off()
    

    ########## plots #############
    ## plot all together 
    # tmp_long <- aq_gaut_weekly %>%
    #   pivot_longer(
    #     cols = c(pm2p5, pm10, o3, so2), # Columns to pivot
    #     names_to = "variable",          # Name of the new column for variable names
    #     values_to = "value"             # Name of the new column for values
    #   )
    
    ## subset by the minimum from the observations
    obs_min <- min(gaut_observations$year_month)
    tmp_long<- subset(tmp_long, tmp_long$year_month>=obs_min & year(tmp_long$year_month)<2021)
    ggplot(tmp_long[which(tmp_long$variable=="pm2p5"),])+
      geom_point(aes(x=year_month, y=value, group=variable, color=dataset, shape = dataset),alpha=0.7)+
      theme_bw()+
      theme(strip.text.y = element_text(angle=0))+
      scale_color_manual(
        name = "Dataset",  # Legend title
        values = c("cams" = "blue", "merra" = "red", "merra_adj" = "darkgreen", "ml1km" = "#AA336A","saaqis" = "black"),
        labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected", "Gapless_1km_ML","SAAQIS (Observations)"))+
      scale_shape_manual(
        name = "Dataset",  # Legend title
        values = c("cams" = 16, "merra" = 17, "merra_adj" = 15, "ml1km" = 3, "saaqis" = 4),
        labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected", "Gapless_1km_ML","SAAQIS (Observations)"))+
      xlab("Year (weekly data)") +
      ylab("Value (ug/m3)") +
      ggtitle("Gauteng")+
      ylim(0,200)+
      facet_grid(variable ~ ., scales="free")+
      theme(legend.position = "none")
    
############# KwaZulu Natal ##############
    ## assign names to provinces
    aq_kzn_monthly <- aq_adm2_monthly[grep("KwaZulu-Natal", aq_adm2_monthly$NAME_1, ignore.case = TRUE),]
    aq_kzn_monthly<-aq_kzn_monthly %>%
      group_by(NAME_1, GID_1, dataset, year_month) %>%
      summarise(
        pm2p5 = mean(pm2p5, na.rm = TRUE),
        pm10 = mean(pm10, na.rm = TRUE),
        o3 = mean(o3, na.rm = TRUE),
        so2 = mean(so2, na.rm = TRUE))
    
    ## pull KwaZulu Natal obs
    kzn_files <- result_table[which(result_table$Province=="KwaZulu Natal"),]$file_names
    kzn_mat<- result_table[which(result_table$Province=="KwaZulu Natal"),]
    obs_kzn_list<-list()
    for(s in 1:length(kzn_files)){
      obs<- fread(paste0("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/observations/",kzn_files[s]))
      obs$site <- kzn_mat[s,"Municipality"]
      if("pm25"%notin%colnames(obs)){
        obs$pm25 <- NA
      }
      if("o3"%notin%colnames(obs)){
        obs$o3 <- NA
      }
      obs<-obs[,c("date","pm25","pm10","o3","so2","site")]
      obs_kzn_list[[s]]<-obs
    }
    ## bind into one data table
    kzn_observations <- rbindlist(obs_kzn_list)
    kzn_observations[, date := as.Date(date, format = "%Y/%m/%d")]
    
    ### take the daily mean 
    kzn_observations <- kzn_observations %>%
      mutate(year_month = floor_date(date, "month")) %>% # Extract year-month
      group_by(year_month) %>%
      summarise(
        pm2p5 = mean(pm25, na.rm = TRUE),
        pm10 = mean(pm10, na.rm = TRUE),
        o3 = mean(o3, na.rm = TRUE),
        so2 = mean(so2, na.rm = TRUE)
      )
    kzn_observations$dataset <- "saaqis"
    ## add province IDs
    ids<-unique(aq_kzn_monthly[,c("GID_1","NAME_1")])
    kzn_observations$GID_1 <- ids$GID_1
    kzn_observations$NAME_1 <- ids$NAME_1
    
    ## reorder to match
    aq_kzn_monthly<-aq_kzn_monthly[,3:ncol(aq_kzn_monthly)]
    cols <- colnames(aq_kzn_monthly)
    kzn_observations<-kzn_observations[,cols]
    kzn_observations$pm2p5[which(kzn_observations$pm2p5>300)] <- 300
    ## bind observations with reanalysis
    aq_kzn_monthly<-rbind(aq_kzn_monthly,kzn_observations)
    aq_kzn_monthly$location <- "KwaZulu-Natal"
    # write.table(aq_kzn_monthly , file= "/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/aq_kzn_monthly_obsreanalysis.csv", quote=FALSE, col.names = TRUE,row.names = FALSE)
    
    
##### PLOT #####
## plot all together ##
    tmp_long <- aq_kzn_monthly %>%
      pivot_longer(
        cols = c(pm2p5, pm10, o3, so2), # Columns to pivot
        names_to = "variable",          # Name of the new column for variable names
        values_to = "value"             # Name of the new column for values
      )

    ## subset by the minimum from the observations
    obs_min <- min(kzn_observations$year_month)
    tmp_long<- subset(tmp_long, tmp_long$year_month>=obs_min & year(tmp_long$year_month)<2021 & tmp_long$value<400)
    # png("/home/sbelman/Documents/BRD/SouthAfrica/airquality/plots/kzn_pm2p5_comparison.png", width=600,height=300)
    kzn_pm2 <- ggplot(tmp_long[which(tmp_long$variable=="pm2p5"),])+
      geom_point(aes(x=year_month, y=value, group=variable,  color=dataset, shape=dataset),size=3,alpha=0.7)+
      theme_bw()+
      theme(strip.text.y = element_text(angle=0))+
      scale_color_manual(
        name = "Dataset",  # Legend title
        values = c("cams" = "blue", "merra" = "red", "merra_adj" = "darkgreen", "ml1km" = "#AA336A","saaqis" = "black"),
        labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected", "Gapless_1km_ML","SAAQIS (Observations)"))+
      scale_shape_manual(
        name = "Dataset",  # Legend title
        values = c("cams" = 16, "merra" = 17, "merra_adj" = 15, "ml1km" = 3, "saaqis" = 4),
        labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected", "Gapless_1km_ML","SAAQIS (Observations)"))+
      xlab("Year (monthly data)") +
      ylab("Value (ug/m3)") +
      # ylim(0,250)+
      ggtitle("KwaZulu Natal")+
      facet_grid(variable ~ ., scales="free")+
      theme(legend.position = "none")
# dev.off()
    
    
    ############# WESTERN CAPE ##############
    ## assign names to provinces
    aq_wesc_monthly <- aq_adm2_monthly[grep("Western Cape", aq_adm2_monthly$NAME_1, ignore.case = TRUE),]
    aq_wesc_monthly<-aq_wesc_monthly %>%
      group_by(NAME_1, GID_1, dataset, year_month) %>%
      summarise(
        pm2p5 = mean(pm2p5, na.rm = TRUE),
        pm10 = mean(pm10, na.rm = TRUE),
        o3 = mean(o3, na.rm = TRUE),
        so2 = mean(so2, na.rm = TRUE))
    
    ## pull KwaZulu Natal obs
    wesc_files <- result_table[which(result_table$Province=="Western Cape"),]$file_names
    wesc_mat<- result_table[which(result_table$Province=="Western Cape"),]
    obs_wesc_list<-list()
    for(s in 1:length(wesc_files)){
      obs<- fread(paste0("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/observations/",wesc_files[s]))
      obs$site <- wesc_mat[s,"Municipality"]
      if("pm25"%notin%colnames(obs)){
        obs$pm25 <- NA
      }
      if("o3"%notin%colnames(obs)){
        obs$o3 <- NA
      }
      obs<-obs[,c("date","pm25","pm10","o3","so2","site")]
      obs_wesc_list[[s]]<-obs
    }
    ## bind into one data table
    wesc_observations <- rbindlist(obs_wesc_list)
    wesc_observations[, date := as.Date(date, format = "%Y/%m/%d")]
    
    ### take the daily mean 
    wesc_observations <- wesc_observations %>%
      mutate(year_month = floor_date(date, "month")) %>% # Extract year-month
      group_by(year_month) %>%
      summarise(
        pm2p5 = mean(pm25, na.rm = TRUE),
        pm10 = mean(pm10, na.rm = TRUE),
        o3 = mean(o3, na.rm = TRUE),
        so2 = mean(so2, na.rm = TRUE)
      )
    wesc_observations$dataset <- "saaqis"
    ## add province IDs
    ids<-unique(aq_wesc_monthly[,c("GID_1","NAME_1")])
    wesc_observations$GID_1 <- ids$GID_1
    wesc_observations$NAME_1 <- ids$NAME_1
    
    ## reorder to match
    aq_wesc_monthly<-aq_wesc_monthly[,3:ncol(aq_wesc_monthly)]
    cols <- colnames(aq_wesc_monthly)
    wesc_observations<-wesc_observations[,cols]
    wesc_observations$pm2p5[which(wesc_observations$pm2p5>300)] <- 300
    ## bind observations with reanalysis
    aq_wesc_monthly<-rbind(aq_wesc_monthly,wesc_observations)
    aq_wesc_monthly$location <- "Western_Cape"
    # write.table(aq_wesc_monthly , file= "/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/aq_wesc_monthly_obsreanalysis.csv", quote=FALSE, col.names = TRUE,row.names = FALSE)
    
    
    ##### PLOT #####
    ## plot all together ##
    tmp_long <- aq_wesc_monthly %>%
      pivot_longer(
        cols = c(pm2p5, pm10, o3, so2), # Columns to pivot
        names_to = "variable",          # Name of the new column for variable names
        values_to = "value"             # Name of the new column for values
      )
    
    ## subset by the minimum from the observations
    obs_min <- min(wesc_observations$year_month)
    tmp_long<- subset(tmp_long, tmp_long$year_month>=obs_min & year(tmp_long$year_month)<2021 & tmp_long$value<400)
    # png("/home/sbelman/Documents/BRD/SouthAfrica/airquality/plots/wesc_pm2p5_comparison.png", width=600,height=300)
    wesc_pm2 <- ggplot(tmp_long[which(tmp_long$variable=="pm2p5"),])+
      geom_point(aes(x=year_month, y=value, group=variable,  color=dataset, shape=dataset),size=3,alpha=0.7)+
      theme_bw()+
      theme(strip.text.y = element_text(angle=0))+
      scale_color_manual(
        name = "Dataset",  # Legend title
        values = c("cams" = "blue", "merra" = "red", "merra_adj" = "darkgreen", "ml1km" = "#AA336A","saaqis" = "black"),
        labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected", "Gapless_1km_ML","SAAQIS (Observations)"))+
      scale_shape_manual(
        name = "Dataset",  # Legend title
        values = c("cams" = 16, "merra" = 17, "merra_adj" = 15, "ml1km" = 3, "saaqis" = 4),
        labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected", "Gapless_1km_ML","SAAQIS (Observations)"))+
      xlab("Year (monthly data)") +
      ylab("Value (ug/m3)") +
      ylim(0,75)+
      ggtitle("Western_Cape")+
      facet_grid(variable ~ ., scales="free")+
      theme(legend.position = "none")
    # dev.off()
    
    
    library(patchwork)
    
    pdf("/home/sbelman/Documents/BRD/SouthAfrica/airquality/plots/obs_versus_reanalysis_4sites.pdf", width = 7, height =5)
    print((gaut_pm2 + mpm_pm2)/(kzn_pm2 + wesc_pm2))
    dev.off()
    
    pdf("/home/sbelman/Documents/BRD/SouthAfrica/airquality/plots/obs_versus_reanalysis_wecsites.pdf", width = 4, height =4)
    print(wesc_pm2)
    dev.off()
    
################################################################################
    ## load reanalysis versus observation data
################################################################################
    aq_kzn_monthly <- fread(file="/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/aq_kzn_monthly_obsreanalysis.csv")
    aq_gertsibande_monthly <- fread(file="/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/aq_gertsibande_monthly_obsreanalysis.csv")
    aq_gaut_monthly <- fread(file="/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/aq_gaut_monthly_obsreanalysis.csv")
    aq_wesc_monthly <- fread(file="/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/aq_wesc_monthly_obsreanalysis.csv")
    

    aq_list <- list(aq_kzn_monthly,aq_gertsibande_monthly,aq_gaut_monthly,aq_wesc_monthly)
    ### measure the correlation between the observed and measured
    # Define error metrics
    mbe <- function(obs, pred) mean(pred - obs, na.rm = TRUE)
    mae <- function(obs, pred) mean(abs(pred - obs), na.rm = TRUE)
    rmse <- function(obs, pred) sqrt(mean((pred - obs)^2, na.rm = TRUE))

    dt_vec <- c("KwaZulu-Natal","Gert_Sibande,Mpumalanga","Gauteng","Western_Cape")
    obs_stats <- NULL
    for(dt in 1:length(dt_vec)){
      df_aq <- aq_list[[dt]]
      saaqis_data <- NULL
        # Filter SAAQIS data as the reference dataset
    saaqis_data <- df_aq %>% 
      filter(dataset == "saaqis") %>% 
      dplyr::select(year_month, pm2p5) 
    year_min <- min(year(saaqis_data$year_month))
    
    # Compute statistics for each dataset
    stat <- df_aq %>%
      filter(dataset != "saaqis", year(year_month)>=year_min) %>% 
      left_join(saaqis_data, by = "year_month", suffix = c("_pred", "_obs"))%>% 
      group_by(dataset) %>%
      summarise(
        Correlation = cor(pm2p5_pred, pm2p5_obs, use = "complete.obs", method= c("spearman")),
        MBE = mbe(pm2p5_obs, pm2p5_pred),
        MAE = mae(pm2p5_obs, pm2p5_pred),
        RMSE = rmse(pm2p5_obs, pm2p5_pred),
        Peak_Timing_Diff = which.max(pm2p5_pred) - which.max(pm2p5_obs)
      ) %>%
      ungroup()
    
    stat$location <- dt_vec[dt]
    obs_stats <- rbind(obs_stats,stat)
    }  


    ########################################################
    ### plot observational points overall
    ########################################################
    library(tidygeocoder)
    name_list <- paste0(result_table$Municipality,",",result_table$District,",",result_table$Province,",","South Africa" )
    # Create a dataframe with location names
    name_list[which(name_list=="club,-gert sibande,Mpumalanga,South Africa")] <- "Secunda, Gert Sibande, Mpumalanga, South Africa"
    name_list[which(name_list=="elandsfontein,-gert sibande,Mpumalanga,South Africa")] <- "Mpumalanga, South Africa"
    name_list[which(name_list=="esikhawini-- rbcaa,king cetshwayo,KwaZulu Natal,South Africa")] <- "KwaZulu-Natal, South Africa"
    name_list[which(name_list=="hammanskraal,-amajuba dm,KwaZulu Natal,South Africa")] <- "amajuba district municipality, KwaZulu-Natal, South Africa"
    name_list[which(name_list=="jabavu,-johannesburg metro,Gauteng,South Africa")] <- "Jabavu, Soweto, Gauteng, South Africa"
    name_list[which(name_list=="saltworks,-amajuba dm,KwaZulu Natal,South Africa")] <- "amajuba district municipality, KwaZulu-Natal, South Africa"
    name_list[which(name_list=="secunda,-gert sibande dm,Mpumalanga,South Africa")] <- "Secunda, Gert Sibande, Mpumalanga, South Africa"
    name_list[which(name_list=="table-view,cape town metro,Western Cape,South Africa")] <- "Table View, Cape Town, Western Cape, South Africa"
    name_list[which(name_list=="wentworth-reservior,ethekwini,KwaZulu Natal,South Africa")] <- "Durban, KwaZulu-Natal, South Africa"
    name_list <- data.frame(name = name_list)
    # Geocode using OpenStreetMap Nominatim
    geo_data <- name_list %>%
      geocode(name, method = "osm")

    color_labels = c("Eastern Cape" = "#008080", "Free State" = "#D75F00", "Gauteng" = "#967bb6" , "KwaZulu-Natal"= "#FFCCCB",
                     "Limpopo" = "#228B22", "Mpumalanga" = "#FFC000", "North West" = "#8B5A2B", "Northern Cape" = "#606060", "Western Cape" = "#1A2D6D" )
    # png("/home/sbelman/Documents/BRD/SouthAfrica/airquality/plots/air_quality_observations_map.png", width=500,height=300)
    ggplot(shp) +
      geom_sf(aes(fill=NAME_1))+
      geom_point(data=geo_data,aes(x=long,y=lat), size=2, fill = "#FF6961", color= "red") +
      scale_fill_manual(values=color_labels)+
      theme_bw()+
      labs(fill="Province")+
      xlab("Longitude")+
      ylab("Latitude") +
      theme(axis.text = element_text(size=12), axis.title = element_text(size=12))
    # dev.off()


################################################################################
### WEEKLY OBSERVATIONAL AIR QUALITY VERSUS DISEASE
################################################################################    
    library(ISOweek)
    library(stringr)
    ## read in disease data
    data <- fread(file="/home/sbelman/Documents/BRD/SouthAfrica/data/sa_adm1_weekly_lag.csv")
    prov_vec <- c("Gauteng", "Mpumalanga", "KwaZulu-Natal")
    result_table <- readRDS(file = "/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/observation_file_list.rds")
    result_table$Province[which(result_table$Province=="KwaZulu Natal")] <- "KwaZulu-Natal"
    obs_aq_list <- list()
    for(i in 1:length(prov_vec)){
          ## subset disease data
          data_prov <- data[which(data$NAME_1==prov_vec[i]),]
          ## assign names to provinces
          aq_gaut_weekly <- aq_adm2_weekly[grep(prov_vec[i], aq_adm2_weekly$NAME_1, ignore.case = TRUE),]
          aq_gaut_weekly <- aq_gaut_weekly %>%
            separate(epiweek, into = c("week", "year"), sep = "-", convert = TRUE) %>%
            mutate(
              week = str_pad(week, width = 2, pad = "0"),  # Pad week with zeroes
              iso_week_str = paste0(year, "-W", week, "-7"),  # Sunday of the week
              date = ISOweek2date(iso_week_str)  # Convert to Date
            )
          aq_gaut_weekly<-aq_gaut_weekly %>%
            group_by(NAME_1, GID_1, dataset, date) %>%
            summarise(
              pm2p5 = mean(pm2p5, na.rm = TRUE),
              pm10 = mean(pm10, na.rm = TRUE))#,
              # o3 = mean(o3, na.rm = TRUE),
              # so2 = mean(so2, na.rm = TRUE))
          
          ## pull gauteng obs
          gaut_files <- result_table[which(result_table$Province==prov_vec[i]),]$file_names
          gaut_mat<- result_table[which(result_table$Province==prov_vec[i]),]
          obs_gaut_list<-list()
          for(s in 1:length(gaut_files)){
            obs<- fread(paste0("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/observations/",gaut_files[s]))
            obs$site <- gaut_mat[s,"Municipality"]
            if("pm25"%notin%colnames(obs)){
              obs$pm25 <- NA
            }
            # obs<-obs[,c("date","pm25","pm10","o3","so2","site")]
            obs<-obs[,c("date","pm25","pm10","site")]
            obs_gaut_list[[s]]<-obs
          }
          ## bind into one data table
          gaut_observations <- rbindlist(obs_gaut_list)
          gaut_observations[, date := as.Date(date, format = "%Y/%m/%d")]
          
          ### take the daily mean 
          gaut_observations2 <- gaut_observations %>%
            mutate(week = floor_date(date, "week", week_start = 7)) %>% # Extract year-month
            group_by(week) %>%
            summarise(
              pm2p5 = mean(pm25, na.rm = TRUE),
              pm10 = mean(pm10, na.rm = TRUE)#,
              # o3 = mean(o3, na.rm = TRUE),
              # so2 = mean(so2, na.rm = TRUE)
            )
          gaut_observations2$dataset <- "saaqis"
          ## add province IDs
          ids<-unique(aq_gaut_weekly[,c("GID_1","NAME_1")])
          gaut_observations2$GID_1 <- ids$GID_1
          gaut_observations2$NAME_1 <- ids$NAME_1
          colnames(gaut_observations2)[which(colnames(gaut_observations2)=="week")] <- "date"
          ## reorder to match
          aq_gaut_weekly2<-aq_gaut_weekly[,3:ncol(aq_gaut_weekly)]
          cols <- colnames(aq_gaut_weekly2)
          gaut_observations2<-gaut_observations2[,cols]
          gaut_observations2$pm2p5[which(gaut_observations2$pm2p5>300)] <- 300
          
          ## bind observations with reanalysis
          aq_gaut_weekly3<-rbind(aq_gaut_weekly2,gaut_observations2)
          aq_gaut_weekly3$location <- prov_vec[i]
          
          # table(data$date%in%aq_gaut_weekly2$date)
          
          ### merge with disease
          data_prov$date <- as.Date(data_prov$date)  # Convert IDate to Date
          gauteng_adm2_weekly <- left_join(aq_gaut_weekly3[which(aq_gaut_weekly3$dataset=="saaqis"),], data_prov[,c("date","disease","GPSC21_count","GPSC1_count","Sero19F_count","population_density", "tasmax_lag0")], by = c("date"))
          obs_aq_list[[i]] <- gauteng_adm2_weekly
         print(ggplot(gauteng_adm2_weekly) +
            geom_line(aes(x=date, y= disease)))
    }
    
    ### all obs
    obs_aq_df <- rbindlist(obs_aq_list)
    # Reshape data to long format
    df_provs <- obs_aq_df %>%
      group_by(date, location) %>%
      summarise(
        pm2p5 = mean(pm2p5, na.rm = TRUE),
        pm10 = mean(pm10, na.rm = TRUE),
        disease = sum(disease, na.rm = TRUE),
        tasmax = mean(tasmax_lag0, na.rm = TRUE),
        population_density = mean(population_density, na.rm = TRUE)
      ) %>%
      arrange(location, date) %>%
      group_by(location) %>%
      mutate(disease_smooth = zoo::rollmean(disease, k = 7, fill = NA, align = "right"))      
  
    df_provs$pm2p5[which(df_provs$pm2p5 > 300)] <- 300
    df_provs$pm10[which(df_provs$pm10 > 300)] <- 300
    
    ### pm10
   ### Compute scale factor per location
    df_scaled <- df_provs %>%
      group_by(location) %>%
      mutate(disease_smooth = zoo::rollmean(disease, k = 7, fill = NA, align = "right")) %>%
      mutate(
        scale_factorpm10 = max(pm10, na.rm = TRUE) / max(disease_smooth, na.rm = TRUE),
        scale_factorpm2p5 = max(pm2p5, na.rm = TRUE) / max(disease_smooth, na.rm = TRUE),
        disease_scaled_pm2p5 = disease_smooth * scale_factorpm2p5,
        disease_scaled_pm10 = disease_smooth * scale_factorpm10
      )
    
    # ## Plot
    (ggplot(df_scaled) +
      geom_line(aes(x = date, y = pm10, color = "pm10", group = location)) +
      geom_line(aes(x = date, y = disease_scaled_pm10, color = "disease", group = location)) +
      facet_grid(location ~ ., scales = "free_y") + # Free y-axis per facet
      scale_color_manual(name = "", values = c("pm10" = "darkred", "disease" = "darkgreen")) +
      xlab("Date") +
      ylab("PM10 concentration") +
      theme_bw() ) +
    
      (ggplot(df_scaled) +
         geom_line(aes(x = date, y = pm2p5, color = "pm2p5", group = location)) +
         geom_line(aes(x = date, y = disease_scaled_pm2p5, color = "disease", group = location)) +
         facet_grid(location ~ ., scales = "free_y") + # Free y-axis per facet
         scale_color_manual(name = "", values = c("pm2p5" = "darkred", "disease" = "darkgreen")) +
         xlab("Date") +
         ylab("PM2.5 concentration") +
         theme_bw() )
    
    
    
    
    ######### decompose
    weekly_data <- gauteng_adm2_weekly %>%
      group_by(date) %>%
      summarise(
        pm2p5 = mean(pm2p5, na.rm = TRUE),
        pm10 = mean(pm10, na.rm = TRUE),
        o3 = mean(o3, na.rm = TRUE),
        so2 = mean(so2, na.rm = TRUE),
        disease = sum(disease, na.rm = TRUE)
      )
    library(zoo)
    
    weekly_data2 <- na.omit(weekly_data)
    weekly_data2 <- weekly_data2 %>%
      mutate(
        pm2p5_ma = rollmean(pm2p5, k = 3, fill = NA, align = "center"),
        disease_ma = rollmean(disease, k = 3, fill = NA, align = "center")
      )
    
    library(pracma)
    ## remove nas
    
    pm2p5_peaks <- findpeaks(weekly_data2$pm2p5_ma, threshold = 10)
    disease_peaks <- findpeaks(weekly_data2$disease_ma, threshold = 1)
    
    # Optional: add peak indicators back to your data
    weekly_data2$pm2p5_peak <- 0
    weekly_data2$disease_peak <- 0
    
    if (!is.null(pm2p5_peaks)) weekly_data2$pm2p5_peak[pm2p5_peaks[,2]] <- 1
    if (!is.null(disease_peaks)) weekly_data2$disease_peak[disease_peaks[,2]] <- 1
    
    scale_factor <- max(weekly_data2$pm2p5_ma, na.rm = TRUE) / max(weekly_data2$disease, na.rm = TRUE)
    
    ggplot(weekly_data2) +
      geom_line(aes(x = date, y = pm2p5_ma, color = "pm2.5"))+
      geom_point(aes(x = date, y = pm2p5_ma, color = "peak"), data = subset(weekly_data2, weekly_data2$pm2p5_peak == 1),
                  size = 2) +
      geom_line(aes(x = date, y = disease * scale_factor, color = "disease"))+
      xlab("Date")+
      theme_bw()+
      ylab("PM2.5 concentration") +
      scale_color_manual(
        name = "",
        values = c("pm2.5" = "darkred", "peak" = "blue", "disease" = "darkgreen")
      ) 
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  ####################################
    ### plot all the stations in Gauteng 
  #####################################
#     ## assign names to provinces
#     aq_gaut_monthly <- aq_adm2_monthly[grep("Gauteng", aq_adm2_monthly$NAME_1, ignore.case = TRUE),]
#     aq_gaut_monthly<-aq_gaut_monthly %>%
#       group_by(NAME_1, GID_1, dataset, year_month) %>%
#       summarise(
#         pm2p5 = mean(pm2p5, na.rm = TRUE),
#         pm10 = mean(pm10, na.rm = TRUE),
#         o3 = mean(o3, na.rm = TRUE),
#         so2 = mean(so2, na.rm = TRUE))
#     aq_gaut_monthly$site <- "reanalysis"
#     ## bind into one data table
#     gaut_observations <- rbindlist(obs_gaut_list)
#     gaut_observations[, date := as.Date(date, format = "%Y/%m/%d")]
#     
#     ### take the daily mean 
#     gaut_observations_sites <- gaut_observations %>%
#       mutate(year_month = floor_date(date, "month")) %>% # Extract year-month
#       group_by(year_month,site) %>%
#       summarise(
#         pm2p5 = mean(pm25, na.rm = TRUE),
#         pm10 = mean(pm10, na.rm = TRUE),
#         o3 = mean(o3, na.rm = TRUE),
#         so2 = mean(so2, na.rm = TRUE)
#       )
#     gaut_observations_sites$dataset <- "saaqis"
#     ## add province IDs
#     ids<-unique(aq_gaut_monthly[,c("GID_1","NAME_1")])
#     gaut_observations_sites$GID_1 <- ids$GID_1
#     gaut_observations_sites$NAME_1 <- ids$NAME_1
#     
#     ## reorder to match
#     aq_gaut_monthly<-aq_gaut_monthly[,3:(ncol(aq_gaut_monthly))]
#     cols <- colnames(aq_gaut_monthly)
#     gaut_observations_sites<-gaut_observations_sites[,cols]
#     
#     ## bind observations with reanalysis
#     aq_gaut_monthly_sites<-rbind(aq_gaut_monthly,gaut_observations_sites)
#     
#     ## plot all together 
#     tmp_long <- aq_gaut_monthly_sites %>%
#       pivot_longer(
#         cols = c(pm2p5, pm10, o3, so2), # Columns to pivot
#         names_to = "variable",          # Name of the new column for variable names
#         values_to = "value"             # Name of the new column for values
#       )
#     
#     ## subset by the minimum from the observations
#     obs_min <- min(gaut_observations_sites$year_month)
#     tmp_long<- subset(tmp_long, tmp_long$year_month>=obs_min & year(tmp_long$year_month)<2021)
#     ggplot(tmp_long[which(tmp_long$variable == "pm2p5"),])+
#       geom_point(aes(x=year_month, y=value, group=variable, color=dataset,group=site, shape=dataset),alpha=0.7)+
#       theme_bw()+
#       theme(strip.text.y = element_text(angle=0))+
#       # scale_color_manual(
#       #   name = "Dataset",  # Legend title
#       #   values = c("cams" = "blue", "merra" = "red", "merra_adj" = "darkgreen", "saaqis" = "black"),
#       #   labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA Adjusted","SAAQIS (Observations)"))+
#       scale_color_manual(
#         name = "Dataset",  # Legend title
#         values = c("cams" = "blue", "merra" = "red", "merra_adj" = "darkgreen", "ml1km" = "#AA336A","saaqis" = "black"),
#         labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected", "Gapless_1km_ML","SAAQIS (Observations)"))+
#       scale_shape_manual(
#         name = "Dataset",  # Legend title
#         values = c("cams" = 16, "merra" = 17, "merra_adj" = 15, "ml1km" = 3, "saaqis" = 4), 
#         labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA) Bias Corrected", "Gapless_1km_ML","SAAQIS (Observations)"))+
#       xlab("Year (monthly data)") +
#       ylab("Value (ug/m3)") +
#       ggtitle("Gauteng")+
#       ylim(0,300)
#     
#     
#     
#     ####################################
#     ### plot all the stations in KZN 
#   #####################################
#     ## assign names to provinces
#     aq_gaut_monthly <- aq_adm2_monthly[grep("KwaZulu-Natal", aq_adm2_monthly$NAME_1, ignore.case = TRUE),]
#     
#     aq_gaut_monthly<-aq_gaut_monthly %>%
#       group_by(NAME_1, GID_1, dataset, year_month) %>%
#       summarise(
#         pm2p5 = mean(pm2p5, na.rm = TRUE),
#         pm10 = mean(pm10, na.rm = TRUE),
#         o3 = mean(o3, na.rm = TRUE),
#         so2 = mean(so2, na.rm = TRUE))
#     aq_gaut_monthly$site <- "reanalysis"
#     ## bind into one data table
#     gaut_observations <- rbindlist(obs_kzn_list)
#     gaut_observations[, date := as.Date(date, format = "%Y/%m/%d")]
#     
#     ### take the daily mean 
#     gaut_observations_sites <- gaut_observations %>%
#       mutate(year_month = floor_date(date, "month")) %>% # Extract year-month
#       group_by(year_month,site) %>%
#       summarise(
#         pm2p5 = mean(pm25, na.rm = TRUE),
#         pm10 = mean(pm10, na.rm = TRUE),
#         o3 = mean(o3, na.rm = TRUE),
#         so2 = mean(so2, na.rm = TRUE)
#       )
#     gaut_observations_sites$dataset <- "saaqis"
#     ## add province IDs
#     ids<-unique(aq_gaut_monthly[,c("GID_1","NAME_1")])
#     gaut_observations_sites$GID_1 <- ids$GID_1
#     gaut_observations_sites$NAME_1 <- ids$NAME_1
#     
#     ## reorder to match
#     aq_gaut_monthly<-aq_gaut_monthly[,3:(ncol(aq_gaut_monthly))]
#     cols <- colnames(aq_gaut_monthly)
#     gaut_observations_sites<-gaut_observations_sites[,cols]
#     
#     ## bind observations with reanalysis
#     aq_gaut_monthly_sites<-rbind(aq_gaut_monthly,gaut_observations_sites)
#     
#     ## plot all together 
#     tmp_long <- aq_gaut_monthly_sites %>%
#       pivot_longer(
#         cols = c(pm2p5, pm10, o3, so2), # Columns to pivot
#         names_to = "variable",          # Name of the new column for variable names
#         values_to = "value"             # Name of the new column for values
#       )
#     
#     ## subset by the minimum from the observations
#     obs_min <- min(gaut_observations_sites$year_month)
#     tmp_long<- subset(tmp_long, tmp_long$year_month>=obs_min & year(tmp_long$year_month)<2021)
#     unique(tmp_long$site)
#     tmp_long$site <- factor(tmp_long$site, levels=rev(unique(tmp_long$site)))
#     p <- ggplot(tmp_long)+
#       geom_point(aes(x=year_month, y=value, group=variable, color=dataset,group=site),alpha=0.7)+
#       theme_bw()+
#       theme(strip.text.y = element_text(angle=0))+
#       scale_color_manual(
#         name = "Dataset",  # Legend title
#         values = c("cams" = "blue", "merra" = "red", "merra_adj" = "darkgreen", "saaqis" = "black"),
#         labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA Adjusted","SAAQIS (Observations)"))+
#       xlab("Year (monthly data)") +
#       ylab("Value (ug/m3)") +
#       ggtitle("KwaZulu-Natal")+
#       ylim(0,300)+
#       facet_grid(variable ~ site, scales="free_y",labeller = label_wrap_gen(width=10)) 
# p    
# 
# 
# 
# ####################################
# ### plot all the stations in Gert Sibande Mpumalanga
# #####################################
# ## assign names to provinces
# aq_gaut_monthly <- aq_adm2_monthly[grep("gert", aq_adm2_monthly$adm2_name, ignore.case = TRUE),]
# # aq_gertsibande_monthly <- aq_adm2_monthly[grep("gert", aq_adm2_monthly$adm2_name, ignore.case = TRUE),]
# 
# aq_gaut_monthly<-aq_gaut_monthly %>%
#   group_by(NAME_1, GID_1, dataset, year_month) %>%
#   summarise(
#     pm2p5 = mean(pm2p5, na.rm = TRUE),
#     pm10 = mean(pm10, na.rm = TRUE),
#     o3 = mean(o3, na.rm = TRUE),
#     so2 = mean(so2, na.rm = TRUE))
# aq_gaut_monthly$site <- "reanalysis"
# ## bind into one data table
# gaut_observations <- rbindlist(obs_gs_list)
# gaut_observations[, date := as.Date(date, format = "%Y/%m/%d")]
# 
# ### take the daily mean 
# gaut_observations_sites <- gaut_observations %>%
#   mutate(year_month = floor_date(date, "month")) %>% # Extract year-month
#   group_by(year_month,site) %>%
#   summarise(
#     pm2p5 = mean(pm25, na.rm = TRUE),
#     pm10 = mean(pm10, na.rm = TRUE),
#     o3 = mean(o3, na.rm = TRUE),
#     so2 = mean(so2, na.rm = TRUE)
#   )
# gaut_observations_sites$dataset <- "saaqis"
# ## add province IDs
# ids<-unique(aq_gaut_monthly[,c("GID_1","NAME_1")])
# gaut_observations_sites$GID_1 <- ids$GID_1
# gaut_observations_sites$NAME_1 <- ids$NAME_1
# 
# ## reorder to match
# aq_gaut_monthly<-aq_gaut_monthly[,3:(ncol(aq_gaut_monthly))]
# cols <- colnames(aq_gaut_monthly)
# gaut_observations_sites<-gaut_observations_sites[,cols]
# 
# ## bind observations with reanalysis
# aq_gaut_monthly_sites<-rbind(aq_gaut_monthly,gaut_observations_sites)
# 
# ## plot all together 
# tmp_long <- aq_gaut_monthly_sites %>%
#   pivot_longer(
#     cols = c(pm2p5, pm10, o3, so2), # Columns to pivot
#     names_to = "variable",          # Name of the new column for variable names
#     values_to = "value"             # Name of the new column for values
#   )
# 
# ## subset by the minimum from the observations
# obs_min <- min(gaut_observations_sites$year_month)
# tmp_long<- subset(tmp_long, tmp_long$year_month>=obs_min & year(tmp_long$year_month)<2021)
# unique(tmp_long$site)
# tmp_long$site <- factor(tmp_long$site, levels=rev(unique(tmp_long$site)))
# p <- ggplot(tmp_long)+
#   geom_point(aes(x=year_month, y=value, group=variable, color=dataset,group=site),alpha=0.7)+
#   theme_bw()+
#   theme(strip.text.y = element_text(angle=0))+
#   scale_color_manual(
#     name = "Dataset",  # Legend title
#     values = c("cams" = "blue", "merra" = "red", "merra_adj" = "darkgreen", "saaqis" = "black"),
#     labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA Adjusted","SAAQIS (Observations)"))+
#   xlab("Year (monthly data)") +
#   ylab("Value (ug/m3)") +
#   ggtitle("Gert Sibande, Mpumalanga")+
#   ylim(0,300)+
#   facet_grid(variable ~ site, scales="free_y",labeller = label_wrap_gen(width=10)) 
# p    
# 
# 
# 
# 
# 
# 
# ##############################################################
# ### plot the difference between the two datasets
# ## around a dashed line at 0
# ##############################################################
# aq_adm1_monthly <- aq_adm2_monthly %>%
#   group_by(year_month, GID_1, NAME_1, dataset) %>%
#   summarize(
#     pm2p5 = mean(pm2p5, na.rm = TRUE),
#     pm10 = mean(pm10, na.rm = TRUE),
#     o3 = mean(o3, na.rm = TRUE),
#     so2 = mean(so2, na.rm = TRUE),
#   ) %>%
#   ungroup()
# #### plot comparison across dates
# aq_adm1_monthly$year <- year(aq_adm1_monthly$year_month)
# tmp<- subset(aq_adm1_monthly, aq_adm1_monthly$year < 2020)
# ### summarize each month
# aq_adm1_monthly$month <- month(aq_adm1_monthly$year_month)
# # Assuming your tibble is named `tmp`
# # Add month column and convert to a factor
# month_labels <- c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")
# month_values <- c(5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4)
# aq_adm1_by_month <- aq_adm1_monthly %>%
#   group_by(month, GID_1, NAME_1, dataset) %>%
#   summarize(
#     pm2p5 = mean(pm2p5, na.rm = TRUE),
#     pm10 = mean(pm10, na.rm = TRUE),
#     o3 = mean(o3, na.rm = TRUE),
#     so2 = mean(so2, na.rm = TRUE),
#   ) %>%
#   pivot_longer(
#     cols = c(pm2p5, pm10, o3, so2), # Columns to pivot
#     names_to = "variable",          # Name of the new column for variable names
#     values_to = "value"             # Name of the new column for values
#   ) %>%
#   mutate(month_names = factor(month, levels = month_values, labels = month_labels))
# 
# 
# tmp_diff <- aq_adm1_by_month %>%
#   pivot_wider(
#     names_from = dataset,       # Use the 'dataset' column to create new columns
#     values_from = value         # Fill new columns with 'value'
#   ) 
# tmp_diff$diff_camsmerra <- tmp_diff$cams-tmp_diff$merra
# tmp_diff$diff_camsmerra_adj <- tmp_diff$cams-tmp_diff$merra_adj
# tmp_diff$diff_merras <- tmp_diff$merra-tmp_diff$merra_adj
# 
# 
# tmp_diff <- tmp_diff %>%
#   mutate(month_names = factor(month, levels = month_values, labels = month_labels))
# # diversity_results$month<-month(diversity_results$year_month)
# 
# 
# diff<- ggplot(tmp_diff)+
#   geom_hline(yintercept=0,linetype="dashed", color="red", linewidth=0.1) +
#   geom_point(aes(x=month_names, y=diff_camsmerra, group=variable, group=NAME_1), size=0.3)+
#   theme_classic()+
#   theme(strip.text.y= element_text(angle=0), axis.text.x = element_text(angle=90, size=9),
#         axis.text.y=element_text(size=12))+
#   xlab("Month") + ylab("Difference (ug/m3)\n(CAMS-MERRA)")+
#   labs(caption = ">0==CAMS overestimate; <0==MERRA overestimate")+
#   facet_grid(variable ~ NAME_1, scales="free")
# 
# diff_merras<- ggplot(subset(tmp_diff,tmp_diff$variable=="pm2p5"))+
#   geom_hline(yintercept=0,linetype="dashed", color="red", linewidth=0.1) +
#   geom_point(aes(x=month_names, y=diff_merras, group=variable, group=NAME_1), size=0.3)+
#   theme_classic()+
#   theme(strip.text.y= element_text(angle=0), axis.text.x = element_text(angle=90, size=9),
#         axis.text.y=element_text(size=12))+
#   xlab("Month") + ylab("Difference (ug/m3)\n(merra-merra_adj)")+
#   labs(caption = ">0==merra overestimate; <0==merra adjusted overestimate")+
#   facet_grid(variable ~ NAME_1, scales="free")
##############################################################
### plot the monthly means from both and plot against each other
##############################################################
# true <- ggplot(aq_adm1_by_month) +
#   geom_point(aes(x=month_names, y=value, color=dataset, group=variable, group=NAME_1), size=0.5)+
#   theme_classic()+
#   theme(strip.text.y= element_text(angle=0), axis.text.x = element_text(angle=90, size=9),
#         axis.text.y=element_text(size=12))+
#   scale_color_manual(
#     name = "Dataset",  # Legend title
#     values = c("cams" = "blue", "merra" = "red", "merra_adj" = "darkgreen"),  
#     labels = c("CAMS (Copernicus)", "MERRA (NASA)", "MERRA (NASA Adjusted)"))+
#   xlab("Month") + ylab("Concentration (ug/m3)")+  
#   facet_grid(variable ~ NAME_1, scales="free") 
# library(patchwork)
# true/diff
# 
    
    
    
    ################## CHECK POINT DATA FOR WEEKLY CYCLE
    # data <- fread("/home/sbelman/Documents/BRD/SouthAfrica/disease/SA_disease_point_base.csv")
    # data$month <- month(data$date)
    # data$year <- year(data$date)
    # dt  <- data %>% 
    #   group_by(month, year) %>%
    #   summarise(disease = sum(disease,na.rm = T))
    # 
    # ggplot(dt) +
    #   geom_line(aes(x=month, y = disease, group =year, color= year)) +
    #   facet_wrap(year~.)
    # 
    # data$agefct <- NA
    # data$agefct <- ifelse(data$ageyears <= 5, "<=5", ifelse(
    #   data$ageyears>5 & data$ageyears <=18, "5-18", ifelse(
    #     data$ageyears>18 & data$ageyears<65, "19-64", ifelse(
    #       data$ageyears>64, ">=65", NA
    #     )
    #   )
    # ))
    # data$agefct <- factor(data$agefct, levels = c("<=5", "5-18","19-64",">=65"))
    # 
    # tab <- table(data$agefct, data$spec_diagnosis)
    # tab/rowSums(tab)
    # 
    # c2 <- chisq.test(tab)
    # print(c2)
    # c2$statistic
    # c2$stdres
    
    
    
    