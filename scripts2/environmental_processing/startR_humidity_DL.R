## ---------------------------
##
## Script name: 02_humidity 2
##
## Purpose of script: Getting the climate variables from era5land
## Author: Daniela Lührsen
## Date Created: 2023-08-01
## Email: daniela.luhrsen@bsc.es
##
## ---------------------------

## load up the packages
packages <- c("startR", "exactextractr", "sf", "progress", "lubridate","CSTools","raster")
install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, require, character.only = TRUE)

## ---------------------------
print("define variables")

# define spatial and temporal extent--------------------------------------------
# local_root="/home/sbelman/Documents/"
local_root=""

country = "SouthAfrica"
lats_min <- -35
lats_max <- -20
lons_min <- 15
lons_max <- 33

start_month <- 1
end_month <- 12
start_year <- 2005
end_year <- 2023

dataset <- "era5land" #era5land, chrips

# obtain correct data format ---------------------------------------------------

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



# load tas ----------------------------------------------------------------
print("tas")
var <- "tas"

# obtain path:
if (dataset == "era5land") {
  ff <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps") {"_f1h"} else {0}
  gg <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps") {"_mean"} else {0}
  path_dataset <- paste0(local_root,"/esarchive/recon/ecmwf/era5land/daily", gg, "/$var$", ff, "/$var$_$sdate$.nc")
  selected_dates <- dates_reanalysis
  forecast <- FALSE
} else if (dataset == "chirps") {
  path_dataset <- paste0(local_root,"/esarchive/obs/ucsb/chirps-v2/monthly_mean/$var$/$var$_$sdate$.nc")
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
  prlr <- data * 3600 * 24 * 30.44 * 1000
  attr(prlr, "Variables")$common$prlr$units <- "mm"
  prlr <- ClimProjDiags::Subset(prlr,
                                along = c("dat", "var"),
                                indices = list(1, 1),
                                drop = "selected")
} else if (var == "ps") {
  ps <- data / 100 #convert Pa to hPa
  attr(ps, "Variables")$common$ps$units <- "hPa"
  ps <- ClimProjDiags::Subset(ps, along = c("dat", "var"),
                              indices = list(1, 1),
                              drop = "selected")
}


# convert to raster
print("tas to raster")
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
tas_raster <- raster::stack(out)


# load ps ----------------------------------------------------------------
# print("surface pressure (known as ps)")
# var <- "ps"
# 
# # obtain path:
# if (dataset == "era5land") {
#   ff <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps" | var == "ps") {"_f1h"} else {0}
#   gg <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps" | var == "ps") {"_mean"} else {0}
#   path_dataset <- paste0("/esarchive/recon/ecmwf/era5land/daily", gg, "/$var$", ff, "/$var$_$sdate$.nc")
#   selected_dates <- dates_reanalysis
#   forecast <- FALSE
# } else if (dataset == "chirps") {
#   path_dataset <- "/esarchive/obs/ucsb/chirps-v2/monthly_mean/$var$/$var$_$sdate$.nc"
#   selected_dates <- dates_reanalysis
# } else {
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
# 
# # transform units
# if (var %in% c("tas", "tasmax", "tasmin", "tdps")) {
#   tas <- data - 273.15
#   attr(tas, "Variables")$common$tas$units <- "C"
#   tas <- ClimProjDiags::Subset(tas,
#                                along = c("dat", "var"),
#                                indices = list(1, 1),
#                                drop = "selected")
# } else if (var == "prlr") {
#   prlr <- data * 3600 * 24 * 30.44 * 1000
#   attr(prlr, "Variables")$common$prlr$units <- "mm"
#   prlr <- ClimProjDiags::Subset(prlr,
#                                 along = c("dat", "var"),
#                                 indices = list(1, 1),
#                                 drop = "selected")
# } else if (var == "ps") {
#   ps <- data / 100 #convert Pa to hPa
#   attr(ps, "Variables")$common$ps$units <- "hPa"
#   ps <- ClimProjDiags::Subset(ps, along = c("dat", "var"),
#                               indices = list(1, 1),
#                               drop = "selected")
# }
# 
# 
# # convert to raster
# print("convert to raster")
# data <- ps
# 
# ps_raster <- NULL
# 
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
# for (x in 1:ntimes) {
#   
#   # Calculate i and j from x
#   d <- ((x - 1) %/% (nyears * nmonths)) + 1
#   m <- ((x - 1) %% nmonths) + 1
#   y <- (((x - 1) %/% nmonths)) %% nyears + 1
#   
#   # Extract date
#   b <- as.character(attr(data, "Variables")[["common"]][["time"]][[x]])
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
# ps_raster <- raster::stack(out)


# load tdps ---------------------------------------------------------------
print("tdps load (dewpoint temp)")
var <- "tdps"

# obtain path:
if (dataset == "era5land") {
  ff <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps") {"_f1h"} else {0}
  gg <- if (var == "tasmax" | var == "tasmin") {""} else if (var == "prlr" | var == "tas" | var == "tdps") {"_mean"} else {0}
  path_dataset <- paste0(local_root,"/esarchive/recon/ecmwf/era5land/daily", gg, "/$var$", ff, "/$var$_$sdate$.nc")
  selected_dates <- dates_reanalysis
  forecast <- FALSE
} else if (dataset == "chirps") {
  path_dataset <- paste0(local_root,"/esarchive/obs/ucsb/chirps-v2/monthly_mean/$var$/$var$_$sdate$.nc")
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
  prlr <- data * 3600 * 24 * 30.44 * 1000
  attr(prlr, "Variables")$common$prlr$units <- "mm"
  prlr <- ClimProjDiags::Subset(prlr,
                                along = c("dat", "var"),
                                indices = list(1, 1),
                                drop = "selected")
} else if (var == "ps") {
  ps <- data / 100 #convert Pa to hPa
  attr(ps, "Variables")$common$ps$units <- "hPa"
  ps <- ClimProjDiags::Subset(ps, along = c("dat", "var"),
                              indices = list(1, 1),
                              drop = "selected")
}

# convert to raster
print("convert to raster")
data <- tas

tdps_raster <- NULL

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
tdps_raster <- raster::stack(out)

# calculate hurs1 (Only hurs1 is used, I'm leaving the others for future reference) -----------
print("calculate hurs")
calculate_hurs1 <- function(tas, tdps) {
  result <- 100 * exp((2500000 / 461.5) * ((1 / (tas + 273.15) - 1 / (tdps + 273.15))))
  return(result)
}

pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = nlayers(tas_raster)
)
hurs1_raster <- stack()
for (day in 1:nlayers(tas_raster)) {
  tas_layer <- tas_raster[[day]]
  tdps_layer <- tdps_raster[[day]]
  print(nlayers(tas_raster))
  print(nlayers(tdps_raster))
  if(nlayers(tas_raster)!=nlayers(tdps_raster)){print("Raster NLayers unequal!")}
  # Apply the function to the layers
  result_layer <- calculate_hurs1(tas_layer, tdps_layer)
  names(result_layer) <- names(tas_raster[[day]])
  
  # Append the result to the result stack
  hurs1_raster <- addLayer(hurs1_raster, result_layer)
  pb$tick()
}

# # calculate hurs2 (Only hurs1 is used, I"m leaving the others for future reference) -----------
# 
# calculate_hurs2 <- function(tas, tdps) {
#   result <- 100 * ((10^(7.5 * tdps / (237.7 + tdps))) / (10^(7.5 * tas / (237.7 + tas))))
#   return(result)
# }
# 
# hurs2_raster <- stack()
# for (day in 1:nlayers(tas_raster)) {
#   tas_layer <- tas_raster[[day]]
#   tdps_layer <- tdps_raster[[day]]
#   
#   # Apply the function to the layers
#   result_layer <- calculate_hurs2(tas_layer, tdps_layer)
#   names(result_layer) <- names(tas_raster[[day]])
#   
#   # Append the result to the result stack
#   hurs2_raster <- addLayer(hurs2_raster, result_layer)
#   pb$tick()
# }
# 
# # calculate hurs3 (Only hurs1 is used, I"m leaving the others for future reference) -----------
# 
# calculate_hurs3 <- function(tas, tdps) {
#   result <- 100 * exp((17.625 * tdps) / (243.04 + tdps)) / exp((17.625 * tas) / (243.04 + tas))
#   return(result)
# }
# pb <- progress_bar$new(
#   format = "[:bar] :percent ETA: :eta",
#   total = nlayers(tas_raster)
# )
# hurs3_raster <- stack()
# for (day in 1:nlayers(tas_raster)) {
#   tas_layer <- tas_raster[[day]]
#   tdps_layer <- tdps_raster[[day]]
#   
#   # Apply the function to the layers
#   result_layer <- calculate_hurs3(tas_layer, tdps_layer)
#   names(result_layer) <- names(tas_raster[[day]])
#   
#   # Append the result to the result stack
#   hurs3_raster <- addLayer(hurs3_raster, result_layer)
#   pb$tick()
# }



# absolute humidity -------------------------------------------------------
print("calculate absolute humidity")
calculate_absh <- function(hurs, tas) {  #https://carnotcycle.wordpress.com/2012/08/04/how-to-convert-relative-humidity-to-absolute-humidity/
  result <- (6.112 * exp((17.67 * tas) / (tas + 243.5)) * hurs * 18.02) / ((273.15 + tas) * 100 * 0.08314)
  return(result)
} # grams/m3

pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = raster::nlayers(tas_raster)
)

absh_raster <- raster::stack()
for (day in 1:raster::nlayers(tas_raster)) {
  tas_layer <- tas_raster[[day]]
  hurs_layer <- hurs1_raster[[day]]
  
  # Apply the function to the layers
  result_layer <- calculate_absh(hurs_layer, tas_layer)
  names(result_layer) <- names(tas_raster[[day]])
  
  # Append the result to the result stack
  absh_raster <- raster::addLayer(absh_raster, result_layer)
  pb$tick()
}


# # specific humidity  -------------------------------------------------------
# print("calculate specific humidity (huss)")
# calculate_huss <- function(tdps, ps) { #https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
#   e <- 6.112 * exp((17.67 * tdps) / (tdps + 243.5))
#   result <- (0.622 * e) / (ps  - (0.378 * e))
#   return(result)
# } # kg/kg
# 
# huss_raster <- raster::stack()
# for (day in 1:raster::nlayers(tdps_raster)) {
#   ps_layer <- ps_raster[[day]]
#   tdps_layer <- tdps_raster[[day]]
#   
#   # Apply the function to the layers
#   result_layer <- calculate_huss(tdps_layer, ps_layer)
#   names(result_layer) <- names(tdps_raster[[day]])
#   
#   # Append the result to the result stack
#   huss_raster <- raster::addLayer(huss_raster, result_layer)
#   print(day)
# }

# obtain data for spatial unit: adm1 -------------------------------------------
print("obtain data and aggregate for spatial units")
# load shapefile and get info
adm1 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_1.shp")
adm1_id <- adm1$GID_1
adm1_name <- adm1$NAME_1
adm1_rows <- nrow(adm1)

# create a dataframe with the results of exactextract for each spatial unit
# loop over the layers to get a value for each week
daily_adm1 <- data.frame()
for (l in 1:raster::nlayers(hurs1_raster)){
  name <- names(hurs1_raster[[l]])
  new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm1_rows),
                       month = rep(substr(name, 7, 8), times = adm1_rows),
                       year = rep(substr(name, 2, 5), times = adm1_rows),
                       adm1_id = adm1_id,
                       adm1_name = adm1_name,
                       hurs = exactextractr::exact_extract(hurs1_raster[[l]], adm1, "mean"),
                       absh = exactextractr::exact_extract(absh_raster[[l]], adm1, "mean"))
                       # huss = exactextractr::exact_extract(huss_raster[[l]], adm1, "mean"))
  daily_adm1 <- rbind(daily_adm1, new_df)
}


# obtain data for spatial unit: adm2 --------------------------------------------------------------------
print("aggregate to admin2")
# load shapefile and get info
adm2 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_2.shp")
adm2_id <- adm2$GID_2
adm2_name <- adm2$NAME_2
adm2_rows <- nrow(adm2)

# create a dataframe with the results of exactextract to get the values for each spatial unit
# loop over the layers to get a value for each week
daily_adm2 <- data.frame()
for (l in 1:raster::nlayers(hurs1_raster)){
  name <- names(hurs1_raster[[l]])
  new_df <- data.frame(day = rep(substr(name, 10, 11), times = adm2_rows),
                       month = rep(substr(name, 7, 8), times = adm2_rows),
                       year = rep(substr(name, 2, 5), times = adm2_rows),
                       adm2_id = adm2_id,
                       adm2_name = adm2_name,
                       hurs = exactextractr::exact_extract(hurs1_raster[[l]], adm2, "mean"),
                       absh = exactextractr::exact_extract(absh_raster[[l]], adm2, "mean"))
                       # huss = exactextractr::exact_extract(huss_raster[[l]], adm2, "mean"))
  daily_adm2 <- rbind(daily_adm2, new_df)
}


# Weekly agregation -------------------------------------------------------
# adm1
daily_adm1$date <- as.Date(paste0(daily_adm1$day, "-", daily_adm1$month, "-", daily_adm1$year), format = "%d-%m-%Y")
daily_adm1$epiweek <- lubridate::epiweek(daily_adm1$date)
daily_adm1$weekyear <- paste0(daily_adm1$epiweek, "-", daily_adm1$year)

weekly_adm1 <- aggregate(list(absh = daily_adm1$absh,
                              hurs = daily_adm1$hurs),
                              # huss = daily_adm1$huss),
                         by = list(epiweek = daily_adm1$weekyear,
                                   adm1_id = daily_adm1$adm1_id,
                                   adm1_name = daily_adm1$adm1_name),
                         FUN = mean, na.rm = TRUE)

# adm2
daily_adm2$date <- as.Date(paste0(daily_adm2$day, "-", daily_adm2$month, "-", daily_adm2$year), format = "%d-%m-%Y")
daily_adm2$epiweek <- lubridate::epiweek(daily_adm2$date)
daily_adm2$weekyear <- paste0(daily_adm2$epiweek, "-", daily_adm2$year)

weekly_adm2 <- aggregate(list(absh = daily_adm2$absh,
                              hurs = daily_adm2$hurs),
                              # huss = daily_adm2$huss),
                         by = list(epiweek = daily_adm2$weekyear,
                                   adm2_id = daily_adm2$adm2_id,
                                   adm2_name = daily_adm2$adm2_name),
                         FUN = mean, na.rm = TRUE)


# save data ---------------------------------------------------------------
print("save data")
write.csv(weekly_adm1, paste0("./",country,"/climate/data/hum_adm1_weekly.csv"), row.names = FALSE)
write.csv(weekly_adm2,  paste0("./",country,"/climate/data/hum_adm2_weekly.csv"), row.names  =  FALSE)
write.csv(daily_adm1, paste0("./",country,"/climate/data/hum_adm1_daily.csv"), row.names = FALSE)
write.csv(daily_adm2,  paste0("./",country,"/climate/data/hum_adm2_daily.csv"), row.names  =  FALSE)
