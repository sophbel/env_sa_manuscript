## ---------------------------
##
## Script name: start_airquality_cams
##
## Purpose of script: Getting the air quality variables variables from CAMS from Copernicus
## Author: Sophie Belman Adapted from Daniela Lührsen
## Date Created: 2024-09-18
## Email: sophie.belman@bsc.es
##
## ---------------------------

## load up the packages
packages <- c("startR", "exactextractr", "sf", "lubridate","raster","data.table")
# install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, require, character.only = TRUE)
print("load packages")
source("./scripts/start_function.R")




################################################################################
#### MERRA DATA FROM NASA
################################################################################
## ---------------------------
### Function to extract raster
# pm2p5_raster_merra<- start_merra(var="pm2p5",dataset="merra",local=TRUE,country="SA", start_year=2005, end_year=2006)
# pm10_raster_merra<- start_merra(var="pm10",dataset="merra",local=TRUE,country="SA", start_year=2005, end_year=2023)
# o3_raster_merra<- start_merra(var="o3",dataset="merra",local=TRUE,country="SA", start_year=2005, end_year=2023)
# so2_raster_merra<- start_merra(var="so2",dataset="merra",local=TRUE,country="SA", start_year=2005, end_year=2023)
# 
# # ## save rasters
# writeRaster(pm2p5_raster_merra, "./SouthAfrica/airquality/data/tifs/pm2p5_raster_merra_stack.tif", format = "raster", overwrite = TRUE)
# writeRaster(pm10_raster_merra, "./SouthAfrica/airquality/data/tifs/pm10_raster_merra_stack.tif", format = "raster", overwrite = TRUE)
# writeRaster(o3_raster_merra, "./SouthAfrica/airquality/data/tifs/o3_raster_merra_stack.tif", format = "raster", overwrite = TRUE)
# writeRaster(so2_raster_merra, "./SouthAfrica/airquality/data/tifs/so2_raster_merra_stack.tif", format = "raster", overwrite = TRUE)
# 
# # ## reload rasters if necessary
# pm2p5_raster_merra <- stack("./SouthAfrica/airquality/data/tifs/pm2p5_raster_merra_stack.grd")
# pm10_raster_merra <- stack("./SouthAfrica/airquality/data/tifs/pm10_raster_merra_stack.grd")
# # o3_raster_merra <- stack("./SouthAfrica/airquality/data/tifs/o3_raster_merra_stack.grd")
# # so2_raster_merra <- stack("./SouthAfrica/airquality/data/tifs/so2_raster_merra_stack.grd")
# 
# 
# 
# # obtain data for spatial unit: adm1 -------------------------------------------
# print("obtain data for spatial units")
# adm1 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_1.shp")
# adm1 <- st_transform(adm1, crs = st_crs(pm2p5_raster_merra))
# adm1_id <- adm1$GID_1
# adm1_name <- adm1$NAME_1
# adm1_rows <- nrow(adm1)
# 
# raster_df<-pm2p5_raster_merra
# daily_adm1 <- data.frame()
# for (l in 1:raster::nlayers(raster_df)) {
#   # for (l in 1:20) {
#   print(paste("Raster Layer",l,"admin1"))
#   name <- names(raster_df[[l]])
#   new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm1_rows),
#                        month = rep(substr(name, 7, 8), times = adm1_rows),
#                        year = rep(substr(name, 2, 5), times = adm1_rows),
#                        adm1_id = adm1_id,
#                        adm1_name = adm1_name,
#                        pm2p5 = exactextractr::exact_extract(pm2p5_raster_merra[[l]], adm1, "mean"),
#                        pm10 = exactextractr::exact_extract(pm10_raster_merra[[l]], adm1, "mean"),
#                        o3 = exactextractr::exact_extract(o3_raster_merra[[l]], adm1, "mean"),
#                        so2 = exactextractr::exact_extract(so2_raster_merra[[l]], adm1, "mean")
#                        )
#   daily_adm1 <- rbind(daily_adm1, new_df)
# }
# 
# # pm10_level70<-daily_adm1$pm10
# # pm10_level50<-daily_adm1$pm10
# # pm10_level1<-daily_adm1$pm10
# # test_levs<-cbind(pm10_level1,pm10_level50,pm10_level70)
# # colnames(test_levs)<-c("lev1","lev50","lev70")
# # test_levs<-data.table(test_levs)
# # library(patchwork)
# # sum(test_levs$lev1)
# # sum(test_levs$lev50)
# # sum(test_levs$lev70)
# 
# # obtain data for spatial unit: adm2 -------------------------------------------
# 
# adm2 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_2.shp")
# adm2_id <- adm2$GID_2
# adm2_name <- adm2$NAME_2
# adm2_rows <- nrow(adm2)
# 
# daily_adm2 <- data.frame()
# for (l in 1:raster::nlayers(raster_df)) {
#   print(paste("Raster Layer",l,"admin2"))
#   
#   name <- names(raster_df[[l]])
#   new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm2_rows),
#                        month = rep(substr(name, 7, 8), times = adm2_rows),
#                        year = rep(substr(name, 2, 5), times = adm2_rows),
#                        adm2_id = adm2_id,
#                        adm2_name = adm2_name,
#                        pm2p5 = exactextractr::exact_extract(pm2p5_raster_merra[[l]], adm2, "mean"),
#                        pm10 = exactextractr::exact_extract(pm10_raster_merra[[l]], adm2, "mean"),
#                        o3 = exactextractr::exact_extract(o3_raster_merra[[l]], adm2, "mean"),
#                        so2 = exactextractr::exact_extract(so2_raster_merra[[l]], adm2, "mean")
#                        )
#   daily_adm2 <- rbind(daily_adm2, new_df)
# }
# 
# 
# 
# # Weekly agregation -------------------------------------------------------
# print("aggregate")
# #adm1
# daily_adm1$date <- as.Date(paste0(daily_adm1$day, "-", daily_adm1$month, "-", daily_adm1$year), format = "%d-%m-%Y")
# daily_adm1$epiweek <- epiweek(daily_adm1$date)
# daily_adm1$weekyear <- paste0(daily_adm1$epiweek, "-", daily_adm1$year)
# 
# weekly_adm1 <- aggregate(list(pm2p5 = daily_adm1$pm2p5,
#                               pm10 = daily_adm1$pm10,
#                               o3 = daily_adm1$o3),
#                          by = list(epiweek = daily_adm1$weekyear,
#                                    adm1_id = daily_adm1$adm1_id,
#                                    adm1_name = daily_adm1$adm1_name),
#                          FUN = mean, na.rm = TRUE)
# 
# # adm2
# daily_adm2$date <- as.Date(paste0(daily_adm2$day, "-", daily_adm2$month, "-", daily_adm2$year), format = "%d-%m-%Y")
# daily_adm2$epiweek <- epiweek(daily_adm2$date)
# daily_adm2$weekyear <- paste0(daily_adm2$epiweek, "-", daily_adm2$year)
# 
# weekly_adm2 <- aggregate(list(pm2p5 = daily_adm2$pm2p5,
#                               pm10 = daily_adm2$pm10,
#                               o3 = daily_adm2$o3,
#                               so2 = daily_adm2$so2),
#                          by = list(epiweek = daily_adm2$weekyear,
#                                    adm2_id = daily_adm2$adm2_id,
#                                    adm2_name = daily_adm2$adm2_name),
#                          FUN = mean, na.rm = TRUE)
# 
# # Save files --------------------------------------------------------------
# print("save")
# write.csv(daily_adm1, paste0("./SouthAfrica/airquality/data/aq_adm1_daily_merra.csv"), row.names = FALSE)
# write.csv(daily_adm2,  paste0("./SouthAfrica/airquality/data/aq_adm2_daily_merra.csv"), row.names  =  FALSE)
# write.csv(weekly_adm1, paste0("./SouthAfrica/airquality/data/aq_adm1_weekly_merra.csv"), row.names = FALSE)
# write.csv(weekly_adm2,  paste0("./SouthAfrica/airquality/data/aq_adm2_weekly_merra.csv"), row.names  =  FALSE)
# 

################################################################################
#### CAMS DATA FROM COPERNICUS
################################################################################
## ---------------------------
## Function to extract raster
# pm2p5_raster_cams<- start_cams(var="pm2p5",dataset="camsra",local=FALSE,country="SA", start_year=2005, end_year=2023)
# saveRDS(pm2p5_raster_cams, "./SouthAfrica/airquality/data/tifs/pm2p5_raster_cams_stack.rds")
# # writeRaster(pm2p5_raster_cams, "./SouthAfrica/airquality/data/tifs/pm2p5_raster_cams_stack.tif", format = "raster", overwrite = TRUE)
# #
# pm10_raster_cams<- start_cams(var="pm10",dataset="camsra",local=FALSE,country="SA", start_year=2005, end_year=2023)
# saveRDS(pm10_raster_cams, "./SouthAfrica/airquality/data/tifs/pm10_raster_cams_stack.rds")
# # writeRaster(pm10_raster_cams, "./SouthAfrica/airquality/data/tifs/pm10_raster_cams_stack.tif", format = "raster", overwrite = TRUE)
# #
o3_raster_cams <- start_cams(var="o3",dataset="camsra",local=FALSE,country="SA", start_year=2005, end_year=2023)
saveRDS(o3_raster_cams, file = "./SouthAfrica/airquality/data/tifs/o3_raster_cams_stack.rds")
# writeRaster(o3_raster_cams, "./SouthAfrica/airquality/data/tifs/o3_raster_cams_stack.tif", format = "raster", overwrite = TRUE)
#
# so2_raster_cams<- start_cams(var="so2",dataset="camsra",local=FALSE,country="SA", start_year=2005, end_year=2023)
# saveRDS(so2_raster_cams, "./SouthAfrica/airquality/data/tifs/so2_raster_cams_stack.rds")
# # writeRaster(so2_raster_cams, "./SouthAfrica/airquality/data/tifs/so2_raster_cams_stack.tif", format = "raster", overwrite = TRUE)

# 
# # ## save rasters
# writeRaster(pm2p5_raster_cams, "./SouthAfrica/airquality/data/tifs/pm2p5_raster_cams_stack.tif", format = "raster", overwrite = TRUE)
# writeRaster(pm10_raster_cams, "./SouthAfrica/airquality/data/tifs/pm10_raster_cams_stack.tif", format = "raster", overwrite = TRUE)
# writeRaster(o3_raster_cams, "./SouthAfrica/airquality/data/tifs/o3_raster_cams_stack.tif", format = "raster", overwrite = TRUE)
# writeRaster(so2_raster_cams, "./SouthAfrica/airquality/data/tifs/so2_raster_cams_stack.tif", format = "raster", overwrite = TRUE)

## reload rasters if necessary
# pm2p5_raster_cams <- stack("./SouthAfrica/airquality/data/tifs/pm2p5_raster_cams_stack.grd")
# pm10_raster_cams <- stack("./SouthAfrica/airquality/data/tifs/pm10_raster_cams_stack.grd")
# o3_raster_cams <- stack("./SouthAfrica/airquality/data/tifs/o3_raster_cams_stack.grd")
# so2_raster_cams <- stack("./SouthAfrica/airquality/data/tifs/so2_raster_cams_stack.grd")

# load RDS
pm2p5_raster_cams <- readRDS("./SouthAfrica/airquality/data/tifs/pm2p5_raster_cams_stack.rds")
pm10_raster_cams <- readRDS("./SouthAfrica/airquality/data/tifs/pm10_raster_cams_stack.rds")
o3_raster_cams <- readRDS("./SouthAfrica/airquality/data/tifs/o3_raster_cams_stack.rds")
so2_raster_cams <- readRDS("./SouthAfrica/airquality/data/tifs/so2_raster_cams_stack.rds")

# check nlayers
nlayers(pm2p5_raster_cams)
nlayers(pm10_raster_cams)
nlayers(o3_raster_cams)
nlayers(so2_raster_cams)

# obtain data for spatial unit: adm1 -------------------------------------------
print("obtain data for spatial units")
adm1 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_1.shp")
adm1 <- st_transform(adm1, crs = st_crs(pm2p5_raster_cams))
adm1_id <- adm1$GID_1
adm1_name <- adm1$NAME_1
adm1_rows <- nrow(adm1)

raster_df<-pm2p5_raster_cams
daily_adm1 <- data.frame()
for (l in 1:raster::nlayers(raster_df)) {
  # for (l in 1:20) {
  print(paste("Raster Layer",l,"admin1"))
  name <- names(raster_df[[l]])
  new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm1_rows),
                       month = rep(substr(name, 7, 8), times = adm1_rows),
                       year = rep(substr(name, 2, 5), times = adm1_rows),
                       adm1_id = adm1_id,
                       adm1_name = adm1_name,
                       pm2p5 = exactextractr::exact_extract(pm2p5_raster_cams[[l]], adm1, "mean"),
                       pm10 = exactextractr::exact_extract(pm10_raster_cams[[l]], adm1, "mean"),
                       o3 = exactextractr::exact_extract(o3_raster_cams[[l]], adm1, "mean"),
                       so2 = exactextractr::exact_extract(so2_raster_cams[[l]], adm1, "mean")
  )
  daily_adm1 <- rbind(daily_adm1, new_df)
}

# pm10_level70<-daily_adm1$pm10
# pm10_level50<-daily_adm1$pm10
# pm10_level1<-daily_adm1$pm10
# test_levs<-cbind(pm10_level1,pm10_level50,pm10_level70)
# colnames(test_levs)<-c("lev1","lev50","lev70")
# test_levs<-data.table(test_levs)
# library(patchwork)
# sum(test_levs$lev1)
# sum(test_levs$lev50)
# sum(test_levs$lev70)

# obtain data for spatial unit: adm2 -------------------------------------------

adm2 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_2.shp")
adm2_id <- adm2$GID_2
adm2_name <- adm2$NAME_2
adm2_rows <- nrow(adm2)

daily_adm2 <- data.frame()
for (l in 1:raster::nlayers(raster_df)) {
  print(paste("Raster Layer",l,"admin2"))

  name <- names(raster_df[[l]])
  new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm2_rows),
                       month = rep(substr(name, 7, 8), times = adm2_rows),
                       year = rep(substr(name, 2, 5), times = adm2_rows),
                       adm2_id = adm2_id,
                       adm2_name = adm2_name,
                       pm2p5 = exactextractr::exact_extract(pm2p5_raster_cams[[l]], adm2, "mean"),
                       pm10 = exactextractr::exact_extract(pm10_raster_cams[[l]], adm2, "mean"),
                       o3 = exactextractr::exact_extract(o3_raster_cams[[l]], adm2, "mean"),
                       so2 = exactextractr::exact_extract(so2_raster_cams[[l]], adm2, "mean")
  )
  daily_adm2 <- rbind(daily_adm2, new_df)
}




# Weekly agregation -------------------------------------------------------
print("aggregate")
#adm1
daily_adm1$date <- as.Date(paste0(daily_adm1$day, "-", daily_adm1$month, "-", daily_adm1$year), format = "%d-%m-%Y")
daily_adm1$epiweek <- epiweek(daily_adm1$date)
daily_adm1$weekyear <- paste0(daily_adm1$epiweek, "-", daily_adm1$year)

weekly_adm1 <- aggregate(list(pm2p5 = daily_adm1$pm2p5,
                              pm10 = daily_adm1$pm10,
                              o3 = daily_adm1$o3,
                              so2 = daily_adm1$so2),
                         by = list(epiweek = daily_adm1$weekyear,
                                   adm1_id = daily_adm1$adm1_id,
                                   adm1_name = daily_adm1$adm1_name),
                         FUN = mean, na.rm = TRUE)

# adm2
daily_adm2$date <- as.Date(paste0(daily_adm2$day, "-", daily_adm2$month, "-", daily_adm2$year), format = "%d-%m-%Y")
daily_adm2$epiweek <- epiweek(daily_adm2$date)
daily_adm2$weekyear <- paste0(daily_adm2$epiweek, "-", daily_adm2$year)

weekly_adm2 <- aggregate(list(pm2p5 = daily_adm2$pm2p5,
                              pm10 = daily_adm2$pm10,
                              o3 = daily_adm2$o3,
                              so2 = daily_adm2$so2),
                         by = list(epiweek = daily_adm2$weekyear,
                                   adm2_id = daily_adm2$adm2_id,
                                   adm2_name = daily_adm2$adm2_name),
                         FUN = mean, na.rm = TRUE)

# Save files --------------------------------------------------------------
print("save")
write.csv(daily_adm2,  paste0("./SouthAfrica/airquality/data/aq_adm2_daily_cams.csv"), row.names  =  FALSE)
write.csv(weekly_adm2,  paste0("./SouthAfrica/airquality/data/aq_adm2_weekly_cams.csv"), row.names  =  FALSE)
write.csv(daily_adm1, paste0("./SouthAfrica/airquality/data/aq_adm1_daily_cams.csv"), row.names = FALSE)
write.csv(weekly_adm1, paste0("./SouthAfrica/airquality/data/aq_adm1_weekly_cams.csv"), row.names = FALSE)

################################################################################
#### ADJUSTED PM2P5 DATA FROM MERRA [2005/1-2023/9]
################################################################################
## ---------------------------
## Function to extract raster
# pm2p5_raster_merra_adj <- start_merra_adjusted(var="pm2p5",dataset="merra_v2-cnnhaqast",local=FALSE,country="SA", start_year=2005, end_year=2023)
# writeRaster(pm2p5_raster_merra_adj, "./SouthAfrica/airquality/data/tifs/pm2p5_raster_merra_adjusted_stack.tif", format = "raster", overwrite = TRUE)
# 
# 
# # pm2p5_raster_merra_adj <- stack("./SouthAfrica/airquality/data/tifs/pm2p5_raster_merra_adjusted_stack.grd")
# 
# # obtain data for spatial unit: adm1 -------------------------------------------
# print("obtain data for spatial units")
# adm1 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_1.shp")
# adm1 <- st_transform(adm1, crs = st_crs(pm2p5_raster_merra_adj))
# adm1_id <- adm1$GID_1
# adm1_name <- adm1$NAME_1
# adm1_rows <- nrow(adm1)
# 
# raster_df<-pm2p5_raster_merra_adj
# daily_adm1 <- data.frame()
# for (l in 1:raster::nlayers(raster_df)) {
#   # for (l in 1:20) {
#   print(paste("Raster Layer",l,"admin1"))
#   name <- names(raster_df[[l]])
#   new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm1_rows),
#                        month = rep(substr(name, 7, 8), times = adm1_rows),
#                        year = rep(substr(name, 2, 5), times = adm1_rows),
#                        adm1_id = adm1_id,
#                        adm1_name = adm1_name,
#                        pm2p5 = exactextractr::exact_extract(pm2p5_raster_merra_adj[[l]], adm1, "mean")
#   )
#   daily_adm1 <- rbind(daily_adm1, new_df)
# }
# 
# 
# # obtain data for spatial unit: adm2 -------------------------------------------
# 
# adm2 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_2.shp")
# adm2_id <- adm2$GID_2
# adm2_name <- adm2$NAME_2
# adm2_rows <- nrow(adm2)
# 
# daily_adm2 <- data.frame()
# for (l in 1:raster::nlayers(raster_df)) {
#   print(paste("Raster Layer",l,"admin2"))
#   
#   name <- names(raster_df[[l]])
#   new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm2_rows),
#                        month = rep(substr(name, 7, 8), times = adm2_rows),
#                        year = rep(substr(name, 2, 5), times = adm2_rows),
#                        adm2_id = adm2_id,
#                        adm2_name = adm2_name,
#                        pm2p5 = exactextractr::exact_extract(pm2p5_raster_merra_adj[[l]], adm2, "mean")
#   )
#   daily_adm2 <- rbind(daily_adm2, new_df)
# }
# 
# write.csv(daily_adm2,  paste0("./SouthAfrica/airquality/data/aq_adm2_daily_merra_adjusted.csv"), row.names  =  FALSE)
# write.csv(daily_adm1, paste0("./SouthAfrica/airquality/data/aq_adm1_daily_merra_adjusted.csv"), row.names = FALSE)
# 
# # Weekly agregation -------------------------------------------------------
# print("aggregate")
# #adm1
# daily_adm1$date <- as.Date(paste0(daily_adm1$day, "-", daily_adm1$month, "-", daily_adm1$year), format = "%d-%m-%Y")
# daily_adm1$epiweek <- epiweek(daily_adm1$date)
# daily_adm1$weekyear <- paste0(daily_adm1$epiweek, "-", daily_adm1$year)
# 
# weekly_adm1 <- aggregate(list(pm2p5 = daily_adm1$pm2p5),
#                          by = list(epiweek = daily_adm1$weekyear,
#                                    adm1_id = daily_adm1$adm1_id,
#                                    adm1_name = daily_adm1$adm1_name),
#                          FUN = mean, na.rm = TRUE)
# 
# # adm2
# daily_adm2$date <- as.Date(paste0(daily_adm2$day, "-", daily_adm2$month, "-", daily_adm2$year), format = "%d-%m-%Y")
# daily_adm2$epiweek <- epiweek(daily_adm2$date)
# daily_adm2$weekyear <- paste0(daily_adm2$epiweek, "-", daily_adm2$year)
# 
# weekly_adm2 <- aggregate(list(pm2p5 = daily_adm2$pm2p5),
#                          by = list(epiweek = daily_adm2$weekyear,
#                                    adm2_id = daily_adm2$adm2_id,
#                                    adm2_name = daily_adm2$adm2_name),
#                          FUN = mean, na.rm = TRUE)
# 
# # Save files --------------------------------------------------------------
# print("save")
# write.csv(weekly_adm1, paste0("./SouthAfrica/airquality/data/reanalysis/aq_adm1_weekly_merra_adjusted.csv"), row.names = FALSE)
# write.csv(weekly_adm2,  paste0("./SouthAfrica/airquality/data/reanalysis/aq_adm2_weekly_merra_adjusted.csv"), row.names  =  FALSE)

################################################################################
#### Gapless 1km2 pm2p5 ML model 2017-2022
################################################################################
## ---------------------------
## Function to extract raster
# pm2p5_raster_ml <- start_mlpm2p5_1km(var="pm2p5",dataset="ml_pm2p5_1km",local=FALSE,country="SA", start_year=2017, end_year=2022)
# writeRaster(pm2p5_raster_ml, "./SouthAfrica/airquality/data/tifs/pm2p5_raster_ml1km_stack.tif", format = "raster", overwrite = TRUE)

# 
# #### Looping through all years
# year_vec <- c(2017:2022)
# for(i in 1:length(year_vec)) {
#   pm2p5_raster_ml <- start_mlpm2p5_1km(var="pm2p5",dataset="ml_pm2p5_1km",local=FALSE,country="SA", start_year=year_vec[i], end_year=year_vec[i])
#   # writeRaster(pm2p5_raster_ml, paste0("./SouthAfrica/airquality/data/tifs/pm2p5_raster_ml1km_stack_",year_vec[i],".tif"), format = "raster", overwrite = TRUE)
#   saveRDS(pm2p5_raster_ml, file =  paste0("./SouthAfrica/airquality/data/tifs/pm2p5_raster_ml1km_stack_",year_vec[i],".rds"))
# 
#     # ## read and stack all files
#     # # Create file paths
#     # file_paths <- paste0("./SouthAfrica/airquality/data/tifs/pm2p5_raster_ml1km_stack_", year_vec, ".tif")
#     # # Read and stack the rasters
#     # pm2p5_raster_ml <- stack(file_paths)
#     # # Write the combined raster to a single file
#     # writeRaster(pm2p5_raster_ml, "./SouthAfrica/airquality/data/tifs/pm2p5_raster_ml1km_stack.grd", format = "raster", overwrite = TRUE)
#     # pm2p5_raster_ml <- stack("./SouthAfrica/airquality/data/tifs/pm2p5_raster_ml_stack.grd")
#     
#     # obtain data for spatial unit: adm1 -------------------------------------------
#     print("obtain data for spatial units")
#     adm1 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_1.shp")
#     adm1 <- st_transform(adm1, crs = st_crs(pm2p5_raster_ml))
#     adm1_id <- adm1$GID_1
#     adm1_name <- adm1$NAME_1
#     adm1_rows <- nrow(adm1)
#     
#     raster_df<-pm2p5_raster_ml
#     daily_adm1 <- data.frame()
#     for (l in 1:raster::nlayers(raster_df)) {
#       # for (l in 1:20) {
#       print(paste("Raster Layer",l,"admin1"))
#       name <- names(raster_df[[l]])
#       new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm1_rows),
#                            month = rep(substr(name, 7, 8), times = adm1_rows),
#                            year = rep(substr(name, 2, 5), times = adm1_rows),
#                            adm1_id = adm1_id,
#                            adm1_name = adm1_name,
#                            pm2p5 = exactextractr::exact_extract(pm2p5_raster_ml[[l]], adm1, "mean")
#       )
#       daily_adm1 <- rbind(daily_adm1, new_df)
#     }
#     
#     
#     # obtain data for spatial unit: adm2 -------------------------------------------
#     
#     adm2 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_2.shp")
#     adm2_id <- adm2$GID_2
#     adm2_name <- adm2$NAME_2
#     adm2_rows <- nrow(adm2)
#     
#     daily_adm2 <- data.frame()
#     for (l in 1:raster::nlayers(raster_df)) {
#       print(paste("Raster Layer",l,"admin2"))
#     
#       name <- names(raster_df[[l]])
#       new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm2_rows),
#                            month = rep(substr(name, 7, 8), times = adm2_rows),
#                            year = rep(substr(name, 2, 5), times = adm2_rows),
#                            adm2_id = adm2_id,
#                            adm2_name = adm2_name,
#                            pm2p5 = exactextractr::exact_extract(pm2p5_raster_ml[[l]], adm2, "mean")
#       )
#       daily_adm2 <- rbind(daily_adm2, new_df)
#     }
# 
# write.csv(daily_adm2,  paste0("./SouthAfrica/airquality/data/aq_adm2_daily_ml1km_",year_vec[i],".csv"), row.names  =  FALSE)
# write.csv(daily_adm1, paste0("./SouthAfrica/airquality/data/aq_adm1_daily_ml1km_",year_vec[i],".csv"), row.names = FALSE)
# }
# 
# 
# 
# # Weekly agregation -------------------------------------------------------
# print("aggregate")
# ### read in each year of the daily data and bind it
# local=TRUE
# daily_adm1_list_tmp <- daily_adm2_list_tmp <- list()
# for(i in 1:length(year_vec)){
#   if(local==TRUE){
#     daily_adm1_list_tmp[[i]] <- fread(paste0("./SouthAfrica/airquality/data/reanalysis/ml1km/aq_adm1_daily_ml1km_",year_vec[i],".csv"))
#     daily_adm2_list_tmp[[i]] <- fread(paste0("./SouthAfrica/airquality/data/reanalysis/ml1km/aq_adm2_daily_ml1km_",year_vec[i],".csv"))
#   }else{
#     daily_adm1_list_tmp[[i]] <- fread(paste0("./SouthAfrica/airquality/data/aq_adm1_daily_ml1km_",year_vec[i],".csv"))
#     daily_adm2_list_tmp[[i]] <- fread(paste0("./SouthAfrica/airquality/data/aq_adm2_daily_ml1km_",year_vec[i],".csv"))
#   }
# }
# daily_adm1_all <- rbindlist(daily_adm1_list_tmp)
# daily_adm2_all <- rbindlist(daily_adm2_list_tmp)
# 
# #adm1
# daily_adm1_all$date <- as.Date(paste0(daily_adm1_all$day, "-", daily_adm1_all$month, "-", daily_adm1_all$year), format = "%d-%m-%Y")
# daily_adm1_all$epiweek <- epiweek(daily_adm1_all$date)
# daily_adm1_all$weekyear <- paste0(daily_adm1_all$epiweek, "-", daily_adm1_all$year)
# 
# weekly_adm1 <- aggregate(list(pm2p5 = daily_adm1_all$pm2p5),
#                          by = list(epiweek = daily_adm1_all$weekyear,
#                                    adm1_id = daily_adm1_all$adm1_id,
#                                    adm1_name = daily_adm1_all$adm1_name),
#                          FUN = mean, na.rm = TRUE)
# 
# # adm2
# daily_adm2_all$date <- as.Date(paste0(daily_adm2_all$day, "-", daily_adm2_all$month, "-", daily_adm2_all$year), format = "%d-%m-%Y")
# daily_adm2_all$epiweek <- epiweek(daily_adm2_all$date)
# daily_adm2_all$weekyear <- paste0(daily_adm2_all$epiweek, "-", daily_adm2_all$year)
# 
# weekly_adm2 <- aggregate(list(pm2p5 = daily_adm2_all$pm2p5),
#                          by = list(epiweek = daily_adm2_all$weekyear,
#                                    adm2_id = daily_adm2_all$adm2_id,
#                                    adm2_name = daily_adm2_all$adm2_name),
#                          FUN = mean, na.rm = TRUE)
# 
# # Save files --------------------------------------------------------------
# print("save")
# write.csv(weekly_adm1, paste0("./SouthAfrica/airquality/data/aq_adm1_weekly_ml1km.csv"), row.names = FALSE)
# write.csv(weekly_adm2,  paste0("./SouthAfrica/airquality/data/aq_adm2_weekly_ml1km.csv"), row.names  =  FALSE)
# write.csv(daily_adm1_all, paste0("./SouthAfrica/airquality/data/aq_adm1_daily_ml1km.csv"), row.names = FALSE)
# write.csv(daily_adm2_all,  paste0("./SouthAfrica/airquality/data/aq_adm2_daily_ml1km.csv"), row.names  =  FALSE)
# #


















################################################################################
#### DEBUG
################################################################################
# 
# data<-readRDS("/home/sbelman/Documents/BRD/SouthAfrica/climate/data/data.rds")
# # Extract date
# b <- substr(data$attrs$Dates[1], 1, 10)


#### find missing dates
# # read rasters
# print("read rasters")
# prlr_raster<- readRDS(paste0("./",country,"/climate/data/prlr_raster.rds"))
# tasmin_raster<- readRDS(paste0("./",country,"/climate/data/tasmin_raster.rds"))
# tasmax_raster<- readRDS(paste0("./",country,"/climate/data/tasmax_raster.rds"))
# tas_raster<- readRDS(paste0("./",country,"/climate/data/tas_raster.rds"))

# ## debugging the naming 
# tas_names<-names(tas_raster)
# tasmin_names<-names(tasmin_raster)
# tasmax_names<-names(tasmax_raster)
# prlr_names<-names(prlr_raster)

# # missing from tasmin
# ## c("X2023.10.31.11.30.00","X2021.02.19","X2021.02.20.12.00.00") # the different format comes from prlr
# tas_names[which(tas_names%notin%tasmin_names)]
# prlr_names[which(prlr_names%notin%tasmin_names)]
# tasmax_names[which(tasmax_names%notin%tasmin_names)]

# # missing from tasmax
# ## c("X2023.10.31.11.30.00","X2021.02.19","X2021.02.20.12.00.00") # the different format comes from prlr
# tas_names[which(tas_names%notin%tasmax_names)]
# prlr_names[which(prlr_names%notin%tasmax_names)]
# tasmin_names[which(tasmin_names%notin%tasmax_names)]

# # missing from prlr
# ## c("X2021.02.19.11.30.00","X2021.02.20.11.30.00")
# tas_names[which(tas_names%notin%prlr_names)]
# tasmin_names[which(tasmin_names%notin%prlr_names)]
# tasmax_names[which(tasmax_names%notin%prlr_names)]

# #missing from tas
# ## c("X2021.02.19","X2021.02.20.12.00.00")# these both come from prlr
# prlr_names[which(prlr_names%notin%tas_names)]
# tasmin_names[which(tasmin_names%notin%tas_names)]
# tasmax_names[which(tasmax_names%notin%tas_names)]