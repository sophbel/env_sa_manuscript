## ---------------------------
##
## Script name: 01_climate
##
## Purpose of script: Getting the climate variables from era5land
## Author: Daniela Lührsen
## Date Created: 2023-08-01
## Email: daniela.luhrsen@bsc.es
##
## ---------------------------

## load up the packages
packages <- c("startR", "exactextractr", "sf", "lubridate","raster","dplyr","tidyr")
install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, require, character.only = TRUE)
print("load packages")
## ---------------------------


# define spatial and temporal extent --------------------------------------
country="SouthAfrica"
# local_path = "/home/sbelman/Documents/"
local_path = ""
lats_min <- -35
lats_max <- -20
lons_min <- 15
lons_max <- 33

start_month <- 1
end_month <- 12
start_year <- 2005
end_year <- 2023

dataset <- "era5land" #era5land, chrips

# obtain correct data format -------------
print("obtain correct data format")
for (mm in start_month:end_month){
  if (mm == start_month) {
    dates <- as.Date(paste0(start_year:end_year, "-", start_month, "-1"))
  } else {
    dates <- append(dates, as.Date(paste0(start_year:end_year, "-", mm, "-1")))
  }
}
dates <- dates[order(dates)]

dates_reanalysis <- array(substr(gsub("-", "", as.character(dates)), 1, 6),
                          dim = c(month = (end_month - start_month + 1),
                                  syear = (end_year - start_year + 1)))
dates_seasonalforecast <- array(paste0(substr(gsub("-", "", as.character(dates)), 1, 6), "01"),
                                dim = c(smonth = (end_month - start_month + 1),
                                        syear = (end_year - start_year + 1)))


# tas ---------------------------------------------------------------------
print("step: sfcWind run")
# var <- "sfcWind"
# # obtain path:
# if (dataset == "era5land") {
#   ff <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps" | var == "sfcWind") {"_f1h"} else {0}
#   gg <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps" | var == "sfcWind") {"_mean"} else {0}
#   path_dataset <- paste0(local_path,"/esarchive/recon/ecmwf/era5land/daily", gg, "/$var$", ff, "/$var$_$sdate$.nc")
#   selected_dates <- dates_reanalysis
#   forecast <- FALSE
# }  else {
#   print("Dataset not known")
# }
# 
# # load data:
# data <- Start(dat = path_dataset,
#               var = var,
#               sdate = selected_dates,
#               time = "all",
#               split_multiselected_dims = TRUE,
#               latitude = startR::values(list(lats_min, lats_max)),
#               latitude_reorder = Sort(decreasing = TRUE),
#               longitude = startR::values(list(lons_min, lons_max)),
#               longitude_reorder = CircularSort(-180, 180),
#               synonims = list(latitude = c("lat", "latitude"),
#                               longitude = c("lon", "longitude")),
#               return_vars = list(latitude = "dat",
#                                  longitude = "dat",
#                                  time = c("sdate")),
#               retrieve = TRUE)
# 
# # transform units
#   sfcWind <- data * 3.6
#   attr(sfcWind, "Variables")$common$sfcWind$units <- "km/hr"
#   sfcWind <- ClimProjDiags::Subset(sfcWind,
#                                 along = c("dat", "var"),
#                                 indices = list(1, 1),
#                                 drop = "selected")
# 
# #convert to raster
# print("convert sfcWind to raster")
# data <- sfcWind
# sfcWind_raster <- NULL
# # Assign latitude and longitude coordinates from data
# lon <- attr(data, "Variables")[["dat1"]][["longitude"]]
# lat <- attr(data, "Variables")[["dat1"]][["latitude"]]
# 
# # Set time variables
# nmonths <- dim(data)[["month"]]
# ndays <- dim(data)[["time"]]
# nyears <- dim(data)[["syear"]]
# ntimes <- nyears * nmonths * ndays
# 
# # Transform multidimensional array into list of rasters
# out <- NULL
# for (x in 1:ntimes){
# 
#   # Calculate i and j from x
#   d <- ((x - 1) %/% (nyears * nmonths)) + 1
#   m <- ((x - 1) %% nmonths) + 1
#   y <- (((x - 1) %/% nmonths)) %% nyears + 1
# 
#   # Extract date
#   b <- as.character(attr(data, "Variables")[["common"]][["time"]][[x]])
#   # Extract date
#   # b <- substr(data$attrs$Dates[x], 1, 10)
#   print(b)
# 
#   if (is.na(b)) {
#     print(paste(x, b))
#     next
#   }
#   # Extract array of values for year-month-day combination
#   a <-  data[m, y, d, , ] #month, year, day
# 
#   # Convert to a raster
#   out[[b]] <- raster::raster(a,
#                              xmn = min(lon), xmx = max(lon),
#                              ymn = min(lat), ymx = max(lat))
# 
#   # Assign CRS
#   crs(out[[b]]) <- "+init=epsg:4326"
# }
# 
# # Convert list of rasters into stack
# sfcWind_raster <- raster::stack(out)
# 
# 
# # # save rasters
# print("save rasters")
# saveRDS(sfcWind_raster, file=paste0("./",country,"/climate/data/sfcWind_raster.rds"))

# # read rasters
# print("read rasters")
sfcWind_raster<- readRDS(paste0("./",country,"/climate/data/sfcWind_raster.rds"))

nlayers(sfcWind_raster)

# obtain data for spatial unit: adm1 -------------------------------------------
print("obtain data for spatial units")
adm1 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_1.shp")
adm1_id <- adm1$GID_1
adm1_name <- adm1$NAME_1
adm1_rows <- nrow(adm1)

daily_adm1 <- data.frame()
for (l in 1:raster::nlayers(sfcWind_raster)) {
  print(paste("Raster Layer",l,"admin1"))
  name <- names(sfcWind_raster[[l]])
  new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm1_rows),
                       month = rep(substr(name, 7, 8), times = adm1_rows),
                       year = rep(substr(name, 2, 5), times = adm1_rows),
                       adm1_id = adm1_id,
                       adm1_name = adm1_name,
                       sfcWind = exactextractr::exact_extract(sfcWind_raster[[l]], adm1, "mean"))
  daily_adm1 <- rbind(daily_adm1, new_df)
}


# obtain data for spatial unit: adm2 -------------------------------------------

adm2 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_2.shp")
adm2_id <- adm2$GID_2
adm2_name <- adm2$NAME_2
adm2_rows <- nrow(adm2)

daily_adm2 <- data.frame()
for (l in 1:raster::nlayers(sfcWind_raster)) {
  print(paste("Raster Layer",l,"admin2"))
  
  name <- names(sfcWind_raster[[l]])
  new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm2_rows),
                       month = rep(substr(name, 7, 8), times = adm2_rows),
                       year = rep(substr(name, 2, 5), times = adm2_rows),
                       adm2_id = adm2_id,
                       adm2_name = adm2_name,
                       sfcWind = exactextractr::exact_extract(sfcWind_raster[[l]], adm2, "mean"))
  daily_adm2 <- rbind(daily_adm2, new_df)
}

# Weekly agregation ------------------------------------------------------------
daily_adm1 <- fread(paste0("./",country,"/climate/data/sfcWind_adm1_daily.csv"))
daily_adm2 <- fread(paste0("./",country,"/climate/data/sfcWind_adm2_daily.csv"))
print("aggregate")
#adm1
daily_adm1$date <- as.Date(paste0(daily_adm1$day, "-", daily_adm1$month, "-", daily_adm1$year), format = "%d-%m-%Y")
daily_adm1$epiweek <- epiweek(daily_adm1$date)
daily_adm1$weekyear <- paste0(daily_adm1$epiweek, "-", daily_adm1$year)

# Aggregate weekly data
weekly_adm1 <- daily_adm1 %>%
  group_by(weekyear, adm1_id, adm1_name) %>% 
  summarise(sfcWind = mean(sfcWind, na.rm =T)) 
  # summarise(
    # sfcWind = exactextractr::exact_extract(sfcWind_raster[[l]], adm1, "mean"))
  

# adm2
daily_adm2$date <- as.Date(paste0(daily_adm2$day, "-", daily_adm2$month, "-", daily_adm2$year), format = "%d-%m-%Y")
daily_adm2$epiweek <- epiweek(daily_adm2$date)
daily_adm2$weekyear <- paste0(daily_adm2$epiweek, "-", daily_adm2$year)

# Aggregate weekly data
weekly_adm2 <- daily_adm2 %>%
  group_by(weekyear, adm2_id, adm2_name) %>%
  summarise(sfcWind = mean(sfcWind, na.rm =T)) 
# 
#   summarise(
#     sfcWind = exactextractr::exact_extract(sfcWind_raster[[l]], adm2, "mean"))

# Save files -------------------------------------------------------------------
print("save")
write.csv(daily_adm1, paste0("./",country,"/climate/data/sfcWind_adm1_daily.csv"), row.names = FALSE)
write.csv(daily_adm2,  paste0("./",country,"/climate/data/sfcWind_adm2_daily.csv"), row.names  =  FALSE)
write.csv(weekly_adm1, paste0("./",country,"/climate/data/sfcWind_adm1_weekly.csv"), row.names = FALSE)
write.csv(weekly_adm2,  paste0("./",country,"/climate/data/sfcWind_adm2_weekly.csv"), row.names  =  FALSE)

