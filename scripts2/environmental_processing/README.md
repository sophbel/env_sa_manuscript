## Environmental processing ##
Each of the startR scripts download multilayered netcdf files, truncate them to South Africa specifically, and process the datasets to the correct concentrations or units. All the datasets are from CAMS with the exception of the air quality datasets.

*startR_spi_DL.R*

The spi script calculates the standardised precipitation and evapotranspiration index (SPEI) and the standardised precipitation index (SPI). These are both drought indicators. 

*startR_climate_DL.R*

Pulls the precipitation and temperature data.

*startR_humidity_DL.R*

Pulls the specific humidity and adjusts to include the absolute and relative humidity values. 

*startR_sfcWind_DL.R*

Pulls the windspeed files.

*startR_airquality_cams.R*

The airquality specific script converts all units to ug/m3 adjusting for the density of the molecules in question. It pulls the data for cams, MERRA, MERRA-2.

## Air Quality Sensitivity ##
*01_explore_aq_reanalysis.R*

The script processes observational data and the reanalysis data from all other sources (ECMWF and NASA) as well as bias corrected and gapless models as described in the manuscript. 
It plots the outcomes and calculates the performance of the reanalysis models given the observational data available. 
