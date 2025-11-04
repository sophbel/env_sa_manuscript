## ---------------------------
##
## Script name: colombia_spi
##
## Purpose of script: Calculate SPI and SPEI for colombia
## Author: Daniela Lührsen
## Date Created: 2024-06-11
## Email: daniela.luhrsen@bsc.es
##
## ---------------------------

## load up the packages
packages <- c("startR", "s2dv", "CSTools", "multiApply", "ClimProjDiags",
              "SPEI", "zoo", "TLMoments", "lmomco", "lmom", "lubridate")
install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, require, character.only = TRUE)

## ---------------------------

# example SPI real data from reanlaysis:

# select parameters:
# local_root="/home/sbelman/Documents/"
local_root=""


country = "SouthAfrica"
lats_min <- -35
lats_max <- -20
lons_min <- 15
lons_max <- 33

start_month <- 1
end_month <- 12
start_year <- 2004
end_year <- 2023


# obtain dates for startR call:
for (mm in start_month:end_month){
  if (mm == start_month) {
    dates <- as.Date(paste0(start_year:end_year, "-", start_month, "-1"))
  } else {
    dates <- append(dates, as.Date(paste0(start_year:end_year, "-", mm, "-1")))
  }
}
dates <- dates[order(dates)]
dates <- array(substr(gsub("-", "", as.character(dates)), 1, 6),
               dim = c(time = (end_month - start_month + 1),
                       syear = (end_year - start_year + 1)))

# precipitation ----
# obtain path for startR call:
path_dataset <- paste0(local_root,"/esarchive/recon/ecmwf/era5land/monthly_mean/$var$_f1h/$var$_$sdate$.nc")

# startR call:
data_prlr <- startR::Start(dat = path_dataset,
                           var = "prlr",
                           sdate = dates,
                           split_multiselected_dims = TRUE,
                           latitude = startR::values(list(lats_min, lats_max)),
                           latitude_reorder = startR::Sort(decreasing = TRUE),
                           longitude = startR::values(list(lons_min, lons_max)),
                           longitude_reorder = startR::CircularSort(-180, 180),
                           synonims = list(latitude = c("lat", "latitude"),
                                           longitude = c("lon", "longitude")),
                           return_vars = list(latitude = "dat",
                                              longitude = "dat",
                                              time = c("sdate")),
                           retrieve = TRUE)
data_prlr <- data_prlr * 3600 * 24 * 30.44 * 1000
attr(data_prlr, "Variables")$common$prlr$units <- "mm"
data_prlr <- CSTools::as.s2dv_cube(data_prlr)
if (!("ensemble" %in% names(data_prlr$data))) {
  data_prlr$data <- s2dv::InsertDim(data_prlr$data, pos =  5, len = 1, name = "ensemble")
}
saveRDS(data_prlr, file=paste0("./",country,"/climate/drought/",country,"_prlr.rds"))

# tasmax ----
path_dataset <- paste0(local_root,"/esarchive/recon/ecmwf/era5land/monthly_mean/$var$_f24h/$var$_$sdate$.nc")

# startR call:
data_tasmax <- startR::Start(dat = path_dataset,
                             var = "tasmax",
                             sdate = dates,
                             split_multiselected_dims = TRUE,
                             latitude = startR::values(list(lats_min, lats_max)),
                             latitude_reorder = startR::Sort(decreasing = TRUE),
                             longitude = startR::values(list(lons_min, lons_max)),
                             longitude_reorder = startR::CircularSort(-180, 180),
                             synonims = list(latitude = c("lat", "latitude"),
                                             longitude = c("lon", "longitude")),
                             return_vars = list(latitude = "dat",
                                                longitude = "dat",
                                                time = c("sdate")),
                             retrieve = TRUE)
data_tasmax <- data_tasmax - 273.15
attr(data_tasmax, "Variables")$common$tasmax$units <- "C"
data_tasmax <- CSTools::as.s2dv_cube(data_tasmax)
if (!("ensemble" %in% names(data_prlr$data))) {
  data_tasmax$data <- s2dv::InsertDim(data_tasmax$data, pos =  5, len = 1, name = "ensemble")
}
saveRDS(data_prlr, file=paste0("./",country,"/climate/drought/",country,"_tasmax.rds"))

# tasmin ----
path_dataset <- paste0(local_root,"/esarchive/recon/ecmwf/era5land/monthly_mean/$var$_f24h/$var$_$sdate$.nc")

# startR call:
data_tasmin <- startR::Start(dat = path_dataset,
                             var = "tasmin",
                             sdate = dates,
                             split_multiselected_dims = TRUE,
                             latitude = startR::values(list(lats_min, lats_max)),
                             latitude_reorder = startR::Sort(decreasing = TRUE),
                             longitude = startR::values(list(lons_min, lons_max)),
                             longitude_reorder = startR::CircularSort(-180, 180),
                             synonims = list(latitude = c("lat", "latitude"),
                                             longitude = c("lon", "longitude")),
                             return_vars = list(latitude = "dat",
                                                longitude = "dat",
                                                time = c("sdate")),
                             retrieve = TRUE)
data_tasmin <- data_tasmin - 273.15
attr(data_tasmin, "Variables")$common$tasmin$units <- "C"
data_tasmin <- CSTools::as.s2dv_cube(data_tasmin)
if (!("ensemble" %in% names(data_tasmin$data))) {
  data_tasmin$data <- s2dv::InsertDim(data_tasmin$data, pos =  5, len = 1, name = "ensemble")
}
saveRDS(data_prlr, file=paste0("./",country,"/climate/drought/",country,"_tasmin.rds"))
# SPI ----
data_prlr <- readRDS(paste0("./",country,"/climate/drought/",country,"_prlr.rds"))

source("https://earth.bsc.es/gitlab/es/csindicators/-/raw/develop-SPEI/R/zzz.R")
source("https://earth.bsc.es/gitlab/es/csindicators/-/raw/develop-SPEI/R/PeriodPET.R")
source("https://earth.bsc.es/gitlab/es/csindicators/-/raw/develop-SPEI/R/PeriodSPEI.R")
# source("colombia-drought/PeriodSPEI.R")


spi1 <- CST_PeriodSPEI(exp = list(pr = data_prlr), pet_exp = array(0, dim = dim(data_prlr$data)), accum = 1, standardization = TRUE)
spi3 <- CST_PeriodSPEI(exp = list(pr = data_prlr), pet_exp = array(0, dim = dim(data_prlr$data)), accum = 3, standardization = TRUE)
spi6 <- CST_PeriodSPEI(exp = list(pr = data_prlr), pet_exp = array(0, dim = dim(data_prlr$data)), accum = 6, standardization = TRUE)
spi12 <- CST_PeriodSPEI(exp = list(pr = data_prlr), pet_exp = array(0, dim = dim(data_prlr$data)), accum = 12, standardization = TRUE)



spei1 <- CST_PeriodSPEI(exp = list(tmax = data_tasmax, tmin = data_tasmin, pr = data_prlr), accum = 1, standardization = TRUE)
spei3 <- CST_PeriodSPEI(exp = list(tmax = data_tasmax, tmin = data_tasmin, pr = data_prlr), accum = 3, standardization = TRUE)
spei6 <- CST_PeriodSPEI(exp = list(tmax = data_tasmax, tmin = data_tasmin, pr = data_prlr), accum = 6, standardization = TRUE)
spei12 <- CST_PeriodSPEI(exp = list(tmax = data_tasmax, tmin = data_tasmin, pr = data_prlr), accum = 12, standardization = TRUE)




# convert to raster -------------------------------------------------------
names <- c("spi1", "spi3", "spi6", "spi12", "spei1", "spei3", "spei6", "spei12")

for (c in 1: length(names)){
  
  if (c == 1) {
    data <- spi1
  } else if (c == 2) {
    data <- spi3
  } else if (c == 3) {
    data <- spi6
  } else if (c == 4) {
    data <- spi12
  } else if (c == 5) {
    data <- spei1
  } else if (c == 6) {
    data <- spei3
  } else if (c == 7) {
    data <- spei6
  } else if (c == 8) {
    data <- spei12
  }
  
  # Assign latitude and longitude coordinates from data
  lon <- data$coords$longitude
  lat <- data$coords$latitude
  
  # Set time variables
  nmonths <- data$dims[["time"]]
  nyears <- data$dims[["syear"]]
  ntimes <- nyears * nmonths
  
  # Transform multidimensional array into list of rasters
  out <- NULL
  for (x in 1:ntimes){
    
    # Calculate i and j from x
    m <- ((x - 1) %% nmonths) + 1
    y <- (((x - 1) %/% nmonths)) %% nyears + 1
    
    # Extract date
    b <- substr(data$attrs$Dates[x], 1, 10)
    
    if (is.na(b)) {
      print(paste(x, b))
      next
    }
    # Extract array of values for year-month-day combination
    a <-  data$data[1, 1, m, y, 1, , ] #month, year, day
    
    # Convert to a raster
    out[[b]] <- raster::raster(a,
                               xmn = min(lon), xmx = max(lon),
                               ymn = min(lat), ymx = max(lat))
    # Assign CRS
    raster::crs(out[[b]]) <- "+init=epsg:4326"
    print(x)
  }
  
  # Convert list of rasters into stack
  assign(paste0(names[c], ".raster"), NULL)
  assign(paste0(names[c], ".raster"), raster::stack(out))
  print(c)
}


# adm1 --------------------------------------------------------------------

adm1 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_1.shp")
adm1_id <- adm1$GID_1
adm1_name <- adm1$NL_NAME_1


ntimestep <- raster::nlayers(spi1.raster) 
monthly_adm1 <- data.frame()
for (l in 1:ntimestep){
  name <- names(spi1.raster[[l]])
  name <- substr(name, 2, 11)
  new_df <- data.frame(month = rep(substr(name, 6, 7), times = length(adm1_id)),
                       year = rep(substr(name, 1, 4), times = length(adm1_id)),
                       adm1_id = adm1_id,
                       adm1_name = adm1_name,
                       spi1 = exactextractr::exact_extract(spi1.raster[[l]], adm1, "mean"),
                       spi3 = exactextractr::exact_extract(spi3.raster[[l]], adm1, "mean"),
                       spi6 = exactextractr::exact_extract(spi6.raster[[l]], adm1, "mean"),
                       spi12 = exactextractr::exact_extract(spi12.raster[[l]], adm1, "mean"),
                       spei1 = exactextractr::exact_extract(spei1.raster[[l]], adm1, "mean"),
                       spei3 = exactextractr::exact_extract(spei3.raster[[l]], adm1, "mean"),
                       spei6 = exactextractr::exact_extract(spei6.raster[[l]], adm1, "mean"),
                       spei12 = exactextractr::exact_extract(spei12.raster[[l]], adm1, "mean"))
  monthly_adm1 <- rbind(monthly_adm1, new_df)
  print(l)
}

# adm2 --------------------------------------------------------------------

adm2 <- sf::st_read("./SouthAfrica/shps/gadm41_ZAF_2.shp")
adm2_id <- adm2$GID_2
adm2_name <- adm2$NL_NAME_2


ntimestep <- raster::nlayers(spi1.raster) 
monthly_adm2 <- data.frame()
for (l in 1:ntimestep){
  name <- names(spi1.raster[[l]])
  name <- substr(name, 2, 11)
  new_df <- data.frame(month = rep(substr(name, 6, 7), times = length(adm2_id)),
                       year = rep(substr(name, 1, 4), times = length(adm2_id)),
                       adm2_id = adm2_id,
                       adm2_name = adm2_name,
                       spi1 = exactextractr::exact_extract(spi1.raster[[l]], adm2, "mean"),
                       spi3 = exactextractr::exact_extract(spi3.raster[[l]], adm2, "mean"),
                       spi6 = exactextractr::exact_extract(spi6.raster[[l]], adm2, "mean"),
                       spi12 = exactextractr::exact_extract(spi12.raster[[l]], adm2, "mean"),
                       spei1 = exactextractr::exact_extract(spei1.raster[[l]], adm2, "mean"),
                       spei3 = exactextractr::exact_extract(spei3.raster[[l]], adm2, "mean"),
                       spei6 = exactextractr::exact_extract(spei6.raster[[l]], adm2, "mean"),
                       spei12 = exactextractr::exact_extract(spei12.raster[[l]], adm2, "mean"))
  monthly_adm2 <- rbind(monthly_adm2, new_df)
  print(l)
}



# save -------------------------------------------------------------------------
write.csv(monthly_adm2, paste0("./",country,"/climate/drought/adm2_droughtindices_monthly.csv"), row.names = FALSE)
write.csv(monthly_adm1, paste0("./",country,"/climate/drought/adm1_droughtindices_monthly.csv"), row.names = FALSE)


# test -------------------------------------------------------------------------

#plot the data
# test <- spei1$data[1,1,1,1,1,,]
# lat <- spei1$coords$latitude
# lon <- spei1$coords$longitude
# red_to_green <- rgb(seq(1, 0, length = 100), seq(0, 1, length = 100), 0)
# png(filename = "test_spei1.png", width = 800, height = 600, res = 100)
# image( test, col= red_to_green,  main = "spei1", xlab = "Longitude", ylab = "Latitude")
# dev.off()


# png("test_plot.png")
# sp::plot(spi3.raster[[3]])
# dev.off()

### END
# 
# ##plot test
# monthly_adm2$month<-as.numeric(monthly_adm2$month)
# monthly_adm2$year<-as.numeric(monthly_adm2$year)
# 
# monthly_adm2$date <- as.Date(paste(monthly_adm2$year, monthly_adm2$month, "15", sep = "-"), "%Y-%m-%d")
# 
# ggplot(monthly_adm2)+
#   geom_line(aes(x=date,y=spei12))
