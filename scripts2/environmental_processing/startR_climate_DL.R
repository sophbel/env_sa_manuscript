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
print("step: tas run")
var <- "tas"
# obtain path:
if (dataset == "era5land") {
  ff <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps") {"_f1h"} else {0}
  gg <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps") {"_mean"} else {0}
  path_dataset <- paste0(local_path,"/esarchive/recon/ecmwf/era5land/daily", gg, "/$var$", ff, "/$var$_$sdate$.nc")
  selected_dates <- dates_reanalysis
  forecast <- FALSE
} else if (dataset == "chirps") {
  path_dataset <- paste0(local_path,"/esarchive/obs/ucsb/chirps-v2/monthly_mean/$var$/$var$_$sdate$.nc")
  selected_dates <- dates_reanalysis
} else {
  print("Dataset not known")
}

# load data:
data <- Start(dat = path_dataset,
              var = var,
              sdate = selected_dates,
              time = "all",
              split_multiselected_dims = TRUE,
              latitude = startR::values(list(lats_min, lats_max)),
              latitude_reorder = Sort(decreasing = TRUE),
              longitude = startR::values(list(lons_min, lons_max)),
              longitude_reorder = CircularSort(-180, 180),
              synonims = list(latitude = c("lat", "latitude"),
                              longitude = c("lon", "longitude")),
              return_vars = list(latitude = "dat",
                                 longitude = "dat",
                                 time = c("sdate")),
              retrieve = TRUE)
# data<-readRDS("/home/sbelman/Documents/BRD/SouthAfrica/climate/data/data.rds")

# transform units
if (var %in% c("tas", "tasmax", "tasmin", "tdps")) {
  tas <- data - 273.15
  attr(tas, "Variables")$common$tas$units <- "C"
  tas <- ClimProjDiags::Subset(tas,
                               along = c("dat", "var"),
                               indices = list(1, 1),
                               drop = "selected")
} else if (var == "prlr") {
  prlr <- data * 3600 * 24  * 1000
  attr(prlr, "Variables")$common$prlr$units <- "mm"
  prlr <- ClimProjDiags::Subset(prlr,
                                along = c("dat", "var"),
                                indices = list(1, 1),
                                drop = "selected")
}

#convert to raster
print("convert tas to raster")
data <- tas

tas_raster <- NULL

# Assign latitude and longitude coordinates from data
lon <- attr(data, "Variables")[["dat1"]][["longitude"]]
lat <- attr(data, "Variables")[["dat1"]][["latitude"]]

# Set time variables
nmonths <- dim(data)[["month"]]
ndays <- dim(data)[["time"]]
nyears <- dim(data)[["syear"]]
ntimes <- nyears * nmonths * ndays

# Transform multidimensional array into list of rasters
out <- NULL
for (x in 1:ntimes){

  # Calculate i and j from x
  d <- ((x - 1) %/% (nyears * nmonths)) + 1
  m <- ((x - 1) %% nmonths) + 1
  y <- (((x - 1) %/% nmonths)) %% nyears + 1

  # Extract date
  b <- as.character(attr(data, "Variables")[["common"]][["time"]][[x]])
  # Extract date
  # b <- substr(data$attrs$Dates[x], 1, 10)
  print(b)

  if (is.na(b)) {
    print(paste(x, b))
    next
  }
  # Extract array of values for year-month-day combination
  a <-  data[m, y, d, , ] #month, year, day

  # Convert to a raster
  out[[b]] <- raster::raster(a,
                             xmn = min(lon), xmx = max(lon),
                             ymn = min(lat), ymx = max(lat))

  # Assign CRS
  crs(out[[b]]) <- "+init=epsg:4326"
}

# Convert list of rasters into stack
tas_raster <- raster::stack(out)

# 
# # tasmin ------------------------------------------------------------------
print("step: tasmin run")
var <- "tasmin"
# obtain path:
if (dataset == "era5land") {
  ff <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps") {"_f1h"} else {0}
  gg <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps") {"_mean"} else {0}
  path_dataset <- paste0(local_path,"/esarchive/recon/ecmwf/era5land/daily", gg, "/$var$", ff, "/$var$_$sdate$.nc")
  selected_dates <- dates_reanalysis
  forecast <- FALSE
} else if (dataset == "chirps") {
  path_dataset <- paste0(local_path,"/esarchive/obs/ucsb/chirps-v2/monthly_mean/$var$/$var$_$sdate$.nc")
  selected_dates <- dates_reanalysis
} else {
  print("Dataset not known")
}

# load data:
data <- Start(dat = path_dataset,
              var = var,
              sdate = selected_dates,
              time = "all",
              split_multiselected_dims = TRUE,
              latitude = startR::values(list(lats_min, lats_max)),
              latitude_reorder = Sort(decreasing = TRUE),
              longitude = startR::values(list(lons_min, lons_max)),
              longitude_reorder = CircularSort(-180, 180),
              synonims = list(latitude = c("lat", "latitude"),
                              longitude = c("lon", "longitude")),
              return_vars = list(latitude = "dat",
                                 longitude = "dat",
                                 time = c("sdate")),
              retrieve = TRUE)

# transform units
if (var %in% c("tas", "tasmax", "tasmin", "tdps")) {
  tas <- data - 273.15
  attr(tas, "Variables")$common$tas$units <- "C"
  tas <- ClimProjDiags::Subset(tas,
                               along = c("dat", "var"),
                               indices = list(1, 1),
                               drop = "selected")
} else if (var == "prlr") {
  prlr <- data * 3600 * 24  * 1000
  attr(prlr, "Variables")$common$prlr$units <- "mm"
  prlr <- ClimProjDiags::Subset(prlr,
                                along = c("dat", "var"),
                                indices = list(1, 1),
                                drop = "selected")
}


# convert to raster
data <- tas

tasmin_raster <- NULL

# Assign latitude and longitude coordinates from data
lon <- attr(data, "Variables")[["dat1"]][["longitude"]]
lat <- attr(data, "Variables")[["dat1"]][["latitude"]]

# Set time variables
nmonths <- dim(data)[["month"]]
ndays <- dim(data)[["time"]]
nyears <- dim(data)[["syear"]]
ntimes <- nyears * nmonths * ndays

# Transform multidimensional array into list of rasters
out <- NULL
for (x in 1:ntimes){

  # Calculate i and j from x
  d <- ((x - 1) %/% (nyears * nmonths)) + 1
  m <- ((x - 1) %% nmonths) + 1
  y <- (((x - 1) %/% nmonths)) %% nyears + 1

  # Extract date
  b <- as.character(attr(data, "Variables")[["common"]][["time"]][[x]])

  if (is.na(b)) {
    print(paste(x, b))
    next
  }
  # Extract array of values for year-month-day combination
  a <-  data[m, y, d, , ] #month, year, day

  # Convert to a raster
  out[[b]] <- raster::raster(a,
                             xmn = min(lon), xmx = max(lon),
                             ymn = min(lat), ymx = max(lat))
  # Assign CRS
  crs(out[[b]]) <- "+init=epsg:4326"
}

# Convert list of rasters into stack
tasmin_raster <- raster::stack(out)

# 
# # tasmax ------------------------------------------------------------------
print("step: tasmax run")
var <- "tasmax"
# obtain path:
if (dataset == "era5land") {
  ff <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps") {"_f1h"} else {0}
  gg <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps") {"_mean"} else {0}
  path_dataset <- paste0(local_path,"/esarchive/recon/ecmwf/era5land/daily", gg, "/$var$", ff, "/$var$_$sdate$.nc")
  selected_dates <- dates_reanalysis
  forecast <- FALSE
} else if (dataset == "chirps") {
  path_dataset <- paste0(local_path,"/esarchive/obs/ucsb/chirps-v2/monthly_mean/$var$/$var$_$sdate$.nc")
  selected_dates <- dates_reanalysis
} else {
  print("Dataset not known")
}

# load data:
data <- Start(dat = path_dataset,
              var = var,
              sdate = selected_dates,
              time = "all",
              split_multiselected_dims = TRUE,
              latitude = startR::values(list(lats_min, lats_max)),
              latitude_reorder = Sort(decreasing = TRUE),
              longitude = startR::values(list(lons_min, lons_max)),
              longitude_reorder = CircularSort(-180, 180),
              synonims = list(latitude = c("lat", "latitude"),
                              longitude = c("lon", "longitude")),
              return_vars = list(latitude = "dat",
                                 longitude = "dat",
                                 time = c("sdate")),
              retrieve = TRUE)


# transform units
if (var %in% c("tas", "tasmax", "tasmin", "tdps")) {
  tas <- data - 273.15
  attr(tas, "Variables")$common$tas$units <- "C"
  tas <- ClimProjDiags::Subset(tas,
                               along = c("dat", "var"),
                               indices = list(1, 1),
                               drop = "selected")
} else if (var == "prlr") {
  prlr <- data * 3600 * 24  * 1000
  attr(prlr, "Variables")$common$prlr$units <- "mm"
  prlr <- ClimProjDiags::Subset(prlr,
                                along = c("dat", "var"),
                                indices = list(1, 1),
                                drop = "selected")
}

# convert to raster
print("convert tasmax to raster")
data <- tas

tasmax_raster <- NULL

# Assign latitude and longitude coordinates from data
lon <- attr(data, "Variables")[["dat1"]][["longitude"]]
lat <- attr(data, "Variables")[["dat1"]][["latitude"]]

# Set time variables
nmonths <- dim(data)[["month"]]
ndays <- dim(data)[["time"]]
nyears <- dim(data)[["syear"]]
ntimes <- nyears * nmonths * ndays

# Transform multidimensional array into list of rasters
out <- NULL
for (x in 1:ntimes){

  # Calculate i and j from x
  d <- ((x - 1) %/% (nyears * nmonths)) + 1
  m <- ((x - 1) %% nmonths) + 1
  y <- (((x - 1) %/% nmonths)) %% nyears + 1

  # Extract date
  b <- as.character(attr(data, "Variables")[["common"]][["time"]][[x]])

  if (is.na(b)) {
    print(paste(x, b))
    next
  }
  # Extract array of values for year-month-day combination
  a <-  data[m, y, d, , ] #month, year, day

  # Convert to a raster
  out[[b]] <- raster::raster(a,
                             xmn = min(lon), xmx = max(lon),
                             ymn = min(lat), ymx = max(lat))

  # Assign CRS
  crs(out[[b]]) <- "+init=epsg:4326"
}

# Convert list of rasters into stack
tasmax_raster <- raster::stack(out)


# prlr --------------------------------------------------------------------
print("step: prlr run")
var <- "prlr"
# obtain path:
if (dataset == "era5land") {
  ff <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps") {"_f1h"} else {0}
  gg <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps") {"_mean"} else {0}
  path_dataset <- paste0(local_path,"/esarchive/recon/ecmwf/era5land/daily", gg, "/$var$", ff, "/$var$_$sdate$.nc")
  selected_dates <- dates_reanalysis
  forecast <- FALSE
} else if (dataset == "chirps") {
  path_dataset <- paste0(local_path,"/esarchive/obs/ucsb/chirps-v2/monthly_mean/$var$/$var$_$sdate$.nc")
  selected_dates <- dates_reanalysis
} else {
  print("Dataset not known")
}

# load data:
data <- Start(dat = path_dataset,
              var = var,
              sdate = selected_dates,
              time = "all",
              split_multiselected_dims = TRUE,
              latitude = startR::values(list(lats_min, lats_max)),
              latitude_reorder = Sort(decreasing = TRUE),
              longitude = startR::values(list(lons_min, lons_max)),
              longitude_reorder = CircularSort(-180, 180),
              synonims = list(latitude = c("lat", "latitude"),
                              longitude = c("lon", "longitude")),
              return_vars = list(latitude = "dat",
                                 longitude = "dat",
                                 time = c("sdate")),
              retrieve = TRUE)

# transform units
if (var %in% c("tas", "tasmax", "tasmin", "tdps")) {
  tas <- data - 273.15
  attr(tas, "Variables")$common$tas$units <- "C"
  tas <- ClimProjDiags::Subset(tas,
                               along = c("dat", "var"),
                               indices = list(1, 1),
                               drop = "selected")
} else if (var == "prlr") {
  prlr <- data * 3600 * 24 * 1000
  attr(prlr, "Variables")$common$prlr$units <- "mm"
  prlr <- ClimProjDiags::Subset(prlr,
                                along = c("dat", "var"),
                                indices = list(1, 1),
                                drop = "selected")
}

# convert to raster
print("convert prlr to raster")
data <- prlr

prlr_raster <- NULL

# Assign latitude and longitude coordinates from data
lon <- attr(data, "Variables")[["dat1"]][["longitude"]]
lat <- attr(data, "Variables")[["dat1"]][["latitude"]]

# Set time variables
nmonths <- dim(data)[["month"]]
ndays <- dim(data)[["time"]]
nyears <- dim(data)[["syear"]]
ntimes <- nyears * nmonths * ndays

# Transform multidimensional array into list of rasters
out <- NULL
for (x in 1:ntimes) {

  # Calculate i and j from x
  d <- ((x - 1) %/% (nyears * nmonths)) + 1
  m <- ((x - 1) %% nmonths) + 1
  y <- (((x - 1) %/% nmonths)) %% nyears + 1

  # Extract date
  b <- as.character(attr(data, "Variables")[["common"]][["time"]][[x]])

  if (is.na(b)) {
    print(paste(x, b))
    next
  }
  # Extract array of values for year-month-day combination
  a <-  data[m, y, d, , ] #month, year, day

  # Convert to a raster
  out[[b]] <- raster::raster(a,
                             xmn = min(lon), xmx = max(lon),
                             ymn = min(lat), ymx = max(lat))

  # Assign CRS
  crs(out[[b]]) <- "+init=epsg:4326"
}

# Convert list of rasters into stack
prlr_raster <- raster::stack(out)
#
# # save rasters
print("save rasters")
saveRDS(tas_raster, file=paste0("./",country,"/climate/data/tas_raster.rds"))
saveRDS(tasmin_raster, file=paste0("./",country,"/climate/data/tasmin_raster.rds"))
saveRDS(tasmax_raster, file=paste0("./",country,"/climate/data/tasmax_raster.rds"))
saveRDS(prlr_raster, file=paste0("./",country,"/climate/data/prlr_raster.rds"))

# # read rasters
# print("read rasters")
prlr_raster<- readRDS(paste0("./",country,"/climate/data/prlr_raster.rds"))
tasmin_raster<- readRDS(paste0("./",country,"/climate/data/tasmin_raster.rds"))
tasmax_raster<- readRDS(paste0("./",country,"/climate/data/tasmax_raster.rds"))
tas_raster<- readRDS(paste0("./",country,"/climate/data/tas_raster.rds"))

nlayers(prlr_raster)
nlayers(tasmin_raster)
nlayers(tasmax_raster)
nlayers(tas_raster)

# obtain data for spatial unit: adm1 -------------------------------------------
print("obtain data for spatial units")
adm1 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_1.shp")
adm1_id <- adm1$GID_1
adm1_name <- adm1$NAME_1
adm1_rows <- nrow(adm1)

daily_adm1 <- data.frame()
for (l in 1:raster::nlayers(tasmin_raster)) {
  # for (l in 6938) {
  print(paste("Raster Layer",l,"admin1"))
  name <- names(tas_raster[[l]])
  new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm1_rows),
                       month = rep(substr(name, 7, 8), times = adm1_rows),
                       year = rep(substr(name, 2, 5), times = adm1_rows),
                       adm1_id = adm1_id,
                       adm1_name = adm1_name,
                       tas = exactextractr::exact_extract(tas_raster[[l]], adm1, "mean"),
                       tasmin = exactextractr::exact_extract(tasmin_raster[[l]], adm1, "mean"),
                       tasmax = exactextractr::exact_extract(tasmax_raster[[l]], adm1, "mean"),
                       prlr = exactextractr::exact_extract(prlr_raster[[l]], adm1, "mean"))
  daily_adm1 <- rbind(daily_adm1, new_df)
}


# obtain data for spatial unit: adm2 -------------------------------------------

adm2 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_2.shp")
adm2_id <- adm2$GID_2
adm2_name <- adm2$NAME_2
adm2_rows <- nrow(adm2)

daily_adm2 <- data.frame()
for (l in 1:raster::nlayers(tasmin_raster)) {
  print(paste("Raster Layer",l,"admin2"))
  
  name <- names(tas_raster[[l]])
  new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm2_rows),
                       month = rep(substr(name, 7, 8), times = adm2_rows),
                       year = rep(substr(name, 2, 5), times = adm2_rows),
                       adm2_id = adm2_id,
                       adm2_name = adm2_name,
                       tas = exactextractr::exact_extract(tas_raster[[l]], adm2, "mean"),
                       tasmin = exactextractr::exact_extract(tasmin_raster[[l]], adm2, "mean"),
                       tasmax = exactextractr::exact_extract(tasmax_raster[[l]], adm2, "mean"),
                       prlr = exactextractr::exact_extract(prlr_raster[[l]], adm2, "mean"))
  daily_adm2 <- rbind(daily_adm2, new_df)
}

# Weekly agregation ------------------------------------------------------------
daily_adm1 <- fread(paste0("./",country,"/climate/data/climate_adm1_daily.csv"))
daily_adm2 <- fread(paste0("./",country,"/climate/data/climate_adm2_daily.csv"))
print("aggregate")
#adm1
daily_adm1$date <- as.Date(paste0(daily_adm1$day, "-", daily_adm1$month, "-", daily_adm1$year), format = "%d-%m-%Y")
daily_adm1$epiweek <- epiweek(daily_adm1$date)
daily_adm1$weekyear <- paste0(daily_adm1$epiweek, "-", daily_adm1$year)

# Aggregate weekly data
weekly_adm1 <- daily_adm1 %>%
  group_by(weekyear, adm1_id, adm1_name) %>%
  summarise(
    tas = mean(tas, na.rm = TRUE),
    tasmin = mean(tasmin, na.rm = TRUE),
    tasmax = mean(tasmax, na.rm = TRUE),
    prlr_mean = mean(prlr, na.rm = TRUE),
    prlr_max = max(prlr, na.rm = TRUE),
    prlr_sum = sum(prlr, na.rm = TRUE)
  )

# adm2
daily_adm2$date <- as.Date(paste0(daily_adm2$day, "-", daily_adm2$month, "-", daily_adm2$year), format = "%d-%m-%Y")
daily_adm2$epiweek <- epiweek(daily_adm2$date)
daily_adm2$weekyear <- paste0(daily_adm2$epiweek, "-", daily_adm2$year)

# weekly_adm2 <- aggregate(list(tas = daily_adm2$tas,
#                               tasmin = daily_adm2$tasmin,
#                               tasmax = daily_adm2$tasmax,
#                               prlr = daily_adm2$prlr),
#                          by = list(epiweek = daily_adm2$weekyear,
#                                    adm2_id = daily_adm2$adm2_id,
#                                    adm2_name = daily_adm2$adm2_name),
#                          FUN = mean, na.rm = TRUE)
# Aggregate weekly data
weekly_adm2 <- daily_adm2 %>%
  group_by(weekyear, adm2_id, adm2_name) %>%
  summarise(
    tas = mean(tas, na.rm = TRUE),
    tasmin = mean(tasmin, na.rm = TRUE),
    tasmax = mean(tasmax, na.rm = TRUE),
    prlr_mean = mean(prlr, na.rm = TRUE),
    prlr_max = max(prlr, na.rm = TRUE),
    prlr_sum = sum(prlr, na.rm = TRUE)
  )

# Save files -------------------------------------------------------------------
print("save")
write.csv(daily_adm1, paste0("./",country,"/climate/data/climate_adm1_daily.csv"), row.names = FALSE)
write.csv(daily_adm2,  paste0("./",country,"/climate/data/climate_adm2_daily.csv"), row.names  =  FALSE)
write.csv(weekly_adm1, paste0("./",country,"/climate/data/climate_adm1_weekly.csv"), row.names = FALSE)
write.csv(weekly_adm2,  paste0("./",country,"/climate/data/climate_adm2_weekly.csv"), row.names  =  FALSE)

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