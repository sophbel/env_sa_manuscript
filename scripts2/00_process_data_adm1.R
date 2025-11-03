################################################################################
###### THIS SCRIPT GENERATES WEEKLY AND MONTHLY DATA AT PROVINCE LEVEL IN SOUTH AFRICA
################################################################################

################################################################################
###### LOAD LIBRARIES AND SOURCE CODE #####
################################################################################

## load libraries and set functions
'%notin%'<-Negate('%in%')
library(lubridate)
library(raster)       # For handling raster data
library(sf)           # For handling vector data (shapefiles)
library(exactextractr) # For extracting raster values based on polygons
library(dplyr)
library(tidyr)
library(data.table)
library(tsModel)
library(ggplot2)
library(INLA)
library(spdep)
library(ISOweek)

source("/home/sbelman/Documents/env_sa_manuscript/scripts2/0_source_functions.R")

## merge province data with shape file area ids
## dist_id, prov_id, dist_name, prov_name
## load disease data
data<-fread("/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_point.csv")
unknowns<-data[which(data$laneid!=""),]
## load shape for adm2 so that I can index and aggreagte the populations to province level
shp_adm1 <-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_1.shp")
shp_adm2 <-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_2.shp")
adm_idx <- shp_adm2[,c("GID_1","GID_2","NAME_1")]
adm_idx <- st_drop_geometry(adm_idx)

## load population size raster from landscan
landscan_raster <- raster("/home/sbelman/Documents/env_sa_manuscript/input_datasets/sociodemographic/landscan-global-2022.tif")
population_data<- fread("/home/sbelman/Documents/env_sa_manuscript/input_datasets/sociodemographic/ZAF_worldpop_2000_2023.csv")

###### TABLE OF CONTENTS #######
############# LINE 543: CLIMATE DATA DEVELOPMENT
############# LINE 905: MERGE CLIMATE AND DISEASE DATA

################################################################################
####### POPULATION SIZE AND DENSITY #####
################################################################################  
# Extract population data for provinces
# district_population <- exact_extract(landscan_raster, shp, fun = "sum")
# Add population data to shapefile
# shp_adm2$population_size <- district_population

# Calculate the area in square kilometers (assuming the shapefile is projected in meters)
shp_adm1 <- shp_adm1 %>%
  mutate(area_km2 = st_area(.) / 1e6)  # Convert from square meters to square kilometers
shp_adm2 <- shp_adm2 %>%
  mutate(area_km2 = st_area(.) / 1e6)  # Convert from square meters to square kilometers

## check match with worldpop names
## adjust the unmatching names
population_data$NAME_2[which(population_data$NAME_2=="O.R.Tambo")] <- "OR Tambo"
population_data$NAME_2[which(population_data$NAME_2=="Uthukela")] <- "uThukela"
population_data$NAME_2[which(population_data$NAME_2=="Uthungulu")] <- "uThungulu"
population_data$NAME_2[which(population_data$NAME_2=="Umzinyathi")] <- "uMzinyathi"
population_data$NAME_2[which(population_data$NAME_2=="Bojanala")] <- "Bojanala Platinum"
population_data$NAME_2[which(population_data$NAME_2=="Umgungundlovu")] <- "uMgungundlovu"
population_data$NAME_2[which(population_data$NAME_2=="Thabo Mofutsanyane")] <- "Thabo Mofutsanyana"
population_data$NAME_2[which(population_data$NAME_2=="Cacadu")] <- "Sarah Baartman"
population_data$NAME_2[which(population_data$NAME_2=="Sisonke")] <- "Harry Gwala"
population_data$NAME_2[which(population_data$NAME_2=="Umkhanyakude")] <- "uMkhanyakude"

pop_name_vec<-unique(population_data$NAME_2)
dist_vec<-unique(data$district)
table(dist_vec%in%pop_name_vec)
dist_vec[dist_vec %notin% pop_name_vec]
pop_name_vec[pop_name_vec %notin% dist_vec]
## rename columns in population data
colnames(population_data)[which(colnames(population_data)=="NAME_2")] <- "district"
colnames(population_data)[which(colnames(population_data)=="GID2")] <- "GID_2"
colnames(population_data)[which(colnames(population_data)=="population")] <- "population_size"

## join shape and population size
population_data <-left_join(population_data, shp_adm2[,c("GID_1","GID_2","area_km2")], by="GID_2")

## turn this into populaiton by year per province
population_data_adm1 <- population_data %>%
  group_by(GID_1,year) %>%
  summarise(population_size = sum(population_size))
population_data_adm1 <-left_join(population_data_adm1, shp_adm1[,c("GID_1","area_km2")], by="GID_1")
population_data_adm1$population_density <- as.numeric(population_data_adm1$population_size/population_data_adm1$area_km2)
## over-write district populations
population_data<-st_drop_geometry(population_data_adm1)
# population_data<-population_data[,-c("geometry")]

## give shape adm1 province name
shp_adm1$province <- shp_adm1$NAME_1

################################################################################
######### MERGE DISEASE AND SHAPEFILE #####
################################################################################   

## merge both with renames etc
## using a single population size for each district
# tmp<-left_join(data,shp[,c("district","GID_2","NAME_1","GID_1","population_size","population_density")],by="district")

## using yearly population sizes from world pop
tmp1<-left_join(data,shp_adm1[,c("province","NAME_1","GID_1")],by="province")
tmp1<-st_drop_geometry(tmp1)
## confirm that only one GID_2 matches each district name
# tmp2<-as.matrix(table(tmp$district,tmp$GID_2))
# tmp2[which(tmp2==0)]<-NA
# isnum<-!is.na(tmp2) ## convert to true false matrix
# rowSums(isnum,na.rm=T) ## all should be 1 
data2<-tmp1[,-c("geometry")]

################################################################################
######  DEFINE FACTORS FOR MODELS ######
################################################################################    
## ensure date is date
data2 <- data2 %>%
  mutate(date = as.Date(date))
## transform serotypes into counts of pcv7, pcv13, and nvt serotypes
# include Vaccine Status
unknown <-c(names(table(data2$serotype))[grep("/",names(table(data2$serotype)))],"NEG38","NEG42","ND","INS","POOL G","POOL G/8","9ALVN","25AF38","","UNTYPABLE")
pcv7<-c("4","6B","9V","14","18C","19F","23F","6E(6B)")
pcv13<-c("1","5","7F","3","6A","19A")
data2$vaccine_status_phen<-ifelse(data2$serotype%in%pcv7,"PCV7",ifelse(
  data2$serotype%in%pcv13,"PCV13",
  ifelse(data2$serotype%in%unknown,"unknown","NVT")
))
##############################################################################  
###### ADD IN GENOMIC DATA ######
##############################################################################  
## include GPS data######
# lane ids that are in GPS
gps1<-fread("/home/sbelman/Documents/BRD/SouthAfrica/genomes/gps_data/gps1_merged_nicd-only.csv")
gps1<-data.frame(gps1)
gps2<-fread("/home/sbelman/Documents/BRD/SouthAfrica/genomes/gps_data/gps2_merged_nicd-only.csv")
gps2<-data.frame(gps2)
gps2.1 <- fread("/home/sbelman/Documents/BRD/SouthAfrica/genomes/gps_data/nicd_2020to2023.csv")
gps2.1<-data.frame(gps2.1)

data2$year[which(data2$laneid%in%gps2.1$Lane_id_y)]

# bind lane IDs
gps_lids<-c(gps1$Lane_id,gps2$Lane_id_x, gps2.1$Lane_id_x)
drug_res<-colnames(gps1)[grepl("SIR",colnames(gps1))][!grepl("colour",colnames(gps1)[grepl("SIR",colnames(gps1))])]
col_names<-c("Lane_id","Region","City","Country","Manifestation","HIV_status","In_silico_serotype","GPSC", drug_res,"PBP1A_2B_2X__autocolour")
colnames(gps2)[which(colnames(gps2)=="Lane_id_x")]<-"Lane_id"
colnames(gps2.1)[which(colnames(gps2.1)=="Lane_id_x")]<-"Lane_id"
colnames(gps2.1)[which(colnames(gps2.1)=="Clinical_manifestation")]<-"Manifestation"

gps1<-gps1[col_names]
gps2<-gps2[col_names]
gps2.1<-gps2.1[col_names]

gps12<-rbind(gps1,gps2,gps2.1)
colnames(gps12)[which(colnames(gps12)=="Lane_id")]<-"laneid"
## find which are in the epi data
data3<-left_join(data2,gps12,by="laneid")
## serotypes
data3$vaccine_status_insil<-ifelse(data3$In_silico_serotype%in%pcv7,"PCV7",ifelse(
  data3$In_silico_serotype%in%pcv13,"PCV13",
  ifelse(data3$In_silico_serotype%in%unknown,"unknown","NVT")
))

## drop unknown district and province
data3 <- subset(data3, data3$province!="unknown")

## INCLUDE EGGNOG AND PANAROO OUTPUTS
# if(COGS==TRUE){
cog_meta<-fread(file="/home/sbelman/Documents/BRD/SouthAfrica/genomes/eggnog/output/cog_meta.csv")
cols<-c("laneid",colnames(cog_meta)[grepl("COG",colnames(cog_meta))])
cog_meta<-cog_meta[,cols,with = FALSE]
## select those with <10000 total
cols<-c("laneid",names(colSums(cog_meta[,-"laneid"])[which(colSums(cog_meta[,-"laneid"])>10000)]))
cog_meta<-cog_meta[,cols,with = FALSE]

# merge with data3
data3<-left_join(data3,cog_meta,by="laneid")
# }else{

# write.table(data3, file="/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_point_base.csv",quote=FALSE,row.names = FALSE, col.names = TRUE,sep=",")
# }
################################################################################
######  CREATE DAILY AGGREGATED DATA FRAME WITH 9 (province=9) NUMBERS FOR EACH DATE ######
################################################################################ 
data2<- fread(file="/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_point_base.csv",quote=FALSE, header = TRUE)


# select which GPSCs and R to include
gpsc_df<-data.table(table(data2$GPSC))
gpsc_vec<-as.vector(gpsc_df[which(gpsc_df$N>100),"V1"])$V1
resistance_columns<-colnames(data2)[grepl("WGS",colnames(data2))]
cog_columns <- colnames(data2)[grepl("COG",colnames(data2))]
# ############DAILY DATA##################
      # Complete the sequence of dates for each province level (GID_1)
        data_summarized <- data2 %>%
          group_by(GID_1,province,date) %>%
          summarize(disease = sum(disease, na.rm = TRUE),
                    sequenced = sum(sequenced, na.rm = TRUE),
                    pcv7 = sum(vaccine_status_phen=="PCV7",na.rm=TRUE),  # summarise phenotypic serotype data
                    pcv13 = sum(vaccine_status_phen=="PCV13",na.rm=TRUE),
                    nvt = sum(vaccine_status_phen=="NVT",na.rm=TRUE),
                    unknown_serotype = sum(vaccine_status_phen=="unknown",na.rm=TRUE),
                    # Creating counts for each GPSC type
                    GPSC1_count = sum(GPSC == "1", na.rm = TRUE),
                    GPSC10_count = sum(GPSC == "10", na.rm = TRUE),
                    GPSC13_count = sum(GPSC == "13", na.rm = TRUE),
                    GPSC14_count = sum(GPSC == "14", na.rm = TRUE),
                    GPSC17_count = sum(GPSC == "17", na.rm = TRUE),
                    GPSC2_count = sum(GPSC == "2", na.rm = TRUE),
                    GPSC21_count = sum(GPSC == "21", na.rm = TRUE),
                    GPSC26_count = sum(GPSC == "26", na.rm = TRUE),
                    GPSC3_count = sum(GPSC == "3", na.rm = TRUE),
                    GPSC5_count = sum(GPSC == "5", na.rm = TRUE),
                    GPSC56_count = sum(GPSC == "56", na.rm = TRUE),
                    GPSC41_count = sum(GPSC == "41", na.rm = TRUE),
                    GPSC8_count = sum(GPSC == "8", na.rm = TRUE),
                    
                    # Counts for  resistance
                    across(all_of(resistance_columns),
                           list(R_count = ~ sum(. == "R", na.rm = TRUE)),
                           .names = "{.col}_{.fn}"),
                    # Counts for COG groups
                    across(all_of(cog_columns),
                           list(count = ~ sum(., na.rm = TRUE)),
                           .names = "{.col}_{.fn}"),

                    .groups = 'drop')
        # ## create a location index to include population sizes etc back in
        # index_locations<-unique(data2[,c("GID_2","year","population_size","population_density","GID_1","province")])
        # data_summarized<-left_join(data_summarized,index_locations,by=c("GID_2","year"))
        data_summarized$date<-as.Date(data_summarized$date)

        # Get the full date range across all data
        full_date_range <- seq.Date(min(data_summarized$date), max(data_summarized$date), by = "day")
        # Get the list of all unique GID_1 and hospital values
        all_provinces <- unique(data_summarized$GID_1)
        all_provinces<-all_provinces[which(!is.na(all_provinces))]
        # Create a complete grid of all dates, districts, and hospitals
        complete_grid <- expand.grid(date = full_date_range, GID_1 = all_provinces)
        
        # Merge this grid with summarized data and fill NAs with 0's
        # Create a lookup table for GID_1 and province
        province_lookup <- data_summarized %>%
          dplyr::select(GID_1, province) %>%
          distinct()
        # Add province column to complete_grid_week
        complete_grid <- complete_grid %>%
          left_join(province_lookup, by = "GID_1")
        df_complete_day <- complete_grid %>%
          left_join(data_summarized, by = c("GID_1", "province","date"))
        # Identify the columns to replace NAs (all except location/time vars)
        non_location_columns <- setdiff(colnames(df_complete_day), c("GID_1", "province", "date"))
        replace_list <- as.list(setNames(rep(0, length(non_location_columns)), non_location_columns))
        df_complete_day <- df_complete_day %>%
          mutate(across(all_of(non_location_columns), ~replace_na(., 0))) %>%
          ungroup()
        
        
        
        ## add in a vaccine column to account for
        df_complete_day$vaccinated <- ifelse(year(df_complete_day$date)<2009,0,1)
        df_complete_day$year <- year(df_complete_day$date)
        ## create a location index to include population sizes etc back in
        df_complete_day<-left_join(df_complete_day,population_data[,c("GID_1","year","population_size")],by=c("GID_1","year"))

        # glance at data
        # ggplot(df_complete_day) + geom_line(aes(x=date,y=population_size, group=province))


############WEEKLY DATA##################
data_summarized_week <- data2 %>%
  mutate(week = floor_date(date, unit = "week")) %>%  # Convert date to start of the week being a Sunday
  group_by(GID_1, province, week) %>%  # Group by GID_1, province, and weekly date
  summarize(disease = sum(disease, na.rm = TRUE),
            sequenced = sum(sequenced, na.rm = TRUE),
            pcv7 = sum(vaccine_status_phen == "PCV7", na.rm = TRUE),  # Count occurrences of PCV7
            pcv13 = sum(vaccine_status_phen == "PCV13", na.rm = TRUE),  # Count occurrences of PCV13
            nvt = sum(vaccine_status_phen == "NVT", na.rm = TRUE),  # Count occurrences of NVT
            unknown_serotype = sum(vaccine_status_phen == "unknown", na.rm = TRUE),  # Count occurrences of unknown serotype
            
            # Counts for Gauteng and Western Cape specifically
            gauteng_count = sum(province == "Gauteng", na.rm=T),
            westerncape_count = sum(province == "Western Cape", na.rm = T),
            
            # Creating counts for ages <1, <5, 5-21,, >65, >80
            age_lt1 = sum(ageyears <1, na.rm=T),
            age_lt6 = sum(ageyears <=5, na.rm=T),
            age_6t14 = sum(ageyears>5 & ageyears<=14,na.rm=T),
            age_15t24 = sum(ageyears>14 & ageyears<=24,na.rm=T),
            age_15t64 = sum(ageyears>14 & ageyears<65,na.rm=T),
            age_18t64 = sum(ageyears>=18 & ageyears<65,na.rm=T),
            age_gt65 = sum(ageyears >65, na.rm=T),
            age_gt80 = sum(ageyears >80, na.rm=T),
            
            # Sex count
            female_count = sum(sex == "FEMALE", na.rm = TRUE),
            male_count = sum(sex == "MALE", na.rm = TRUE),
            unknownsex_count = sum(sex == "UNKNOWN", na.rm = TRUE),
            
            # Age and Sex Count
            age_lt1_fem = sum(ageyears <1 & sex=="FEMALE", na.rm=T),
            age_lt5_fem = sum(ageyears <=5 & sex=="FEMALE", na.rm=T),
            age_gt65_fem = sum(ageyears > 65 & sex=="FEMALE", na.rm=T),
            age_gt80_fem = sum(ageyears > 80 & sex=="FEMALE", na.rm=T),
            age_lt1_mal = sum(ageyears <1 & sex=="MALE", na.rm=T),
            age_lt5_mal = sum(ageyears <=5 & sex=="MALE", na.rm=T),
            age_gt65_mal = sum(ageyears > 65 & sex=="MALE", na.rm=T),
            age_gt80_mal = sum(ageyears > 80 & sex=="MALE", na.rm=T),
            
            # Smoking count
            smoker_count = sum(smoker == "YES", na.rm = TRUE),
            nonsmoker_count = sum(smoker == "NO", na.rm = TRUE),
            unknownsmoker_count = sum(smoker == "", na.rm = TRUE),
            
            # Pregnancy count
            pregnant_count = sum(pregnant == "Y", na.rm = TRUE),
            notpregnant_count = sum(pregnant == "N", na.rm = TRUE),
            unknownpregnant_count = sum(pregnant %notin% c("Y","N"), na.rm = TRUE),
            
            ### mening, bacteremia, other count
            mening_count = sum(spec_diagnosis == "Meningitis", na.rm = TRUE),
            bact_count = sum(spec_diagnosis == "Bacteremia without focus", na.rm = TRUE),
            other_count = sum(spec_diagnosis == "Other", na.rm = TRUE),
            
            # Creating counts for each GPSC type
            GPSC1_count = sum(GPSC == "1", na.rm = TRUE),
            GPSC10_count = sum(GPSC == "10", na.rm = TRUE),
            GPSC13_count = sum(GPSC == "13", na.rm = TRUE),
            GPSC14_count = sum(GPSC == "14", na.rm = TRUE),
            GPSC17_count = sum(GPSC == "17", na.rm = TRUE),
            GPSC2_count = sum(GPSC == "2", na.rm = TRUE),
            GPSC21_count = sum(GPSC == "21", na.rm = TRUE),
            GPSC26_count = sum(GPSC == "26", na.rm = TRUE),
            GPSC3_count = sum(GPSC == "3", na.rm = TRUE),
            GPSC5_count = sum(GPSC == "5", na.rm = TRUE),
            GPSC56_count = sum(GPSC == "56", na.rm = TRUE),
            GPSC41_count = sum(GPSC == "41", na.rm = TRUE),
            GPSC8_count = sum(GPSC == "8", na.rm = TRUE),
            
            # Creating counts for each serotype
            # PCV7
            Sero4_count = sum(serotype == "4", na.rm = TRUE),
            Sero6B_count = sum(serotype == "6B", na.rm = TRUE),
            Sero9V_count = sum(serotype == "9V", na.rm = TRUE),
            Sero14_count = sum(serotype == "14", na.rm = TRUE),
            Sero18C_count = sum(serotype == "18C", na.rm = TRUE),
            Sero19F_count = sum(serotype == "19F", na.rm = TRUE),
            Sero23F_count = sum(serotype == "23F", na.rm = TRUE),
            # PCV13
            Sero1_count = sum(serotype == "1", na.rm = TRUE),
            Sero3_count = sum(serotype == "3", na.rm = TRUE),
            # Sero5_count = sum(serotype == "5", na.rm = TRUE), ## none in this dataset
            Sero6A_count = sum(serotype == "6A", na.rm = TRUE),
            Sero7F_count = sum(serotype == "7F", na.rm = TRUE),
            Sero19A_count = sum(serotype == "19A", na.rm = TRUE),
            # NVT
            Sero12F_count = sum(serotype == "12F", na.rm = TRUE),
            Sero15A_count = sum(serotype == "15A", na.rm = TRUE),
            Sero15BC_count = sum(serotype == "15B/C", na.rm = TRUE), 
            Sero16F_count = sum(serotype == "16F", na.rm = TRUE),
            Sero8_count = sum(serotype == "8", na.rm = TRUE),
            Sero9N_count = sum(serotype == "9N", na.rm = TRUE),
            
            ## in silico serotypes
            # Sero1_count = sum(In_silico_serotype == "1", na.rm = TRUE),
            # Sero10A_count = sum(In_silico_serotype == "10A", na.rm = TRUE),
            # Sero12F_count = sum(In_silico_serotype == "12F", na.rm = TRUE),
            # Sero13_count = sum(In_silico_serotype == "13", na.rm = TRUE),
            # Sero14_count = sum(In_silico_serotype == "14", na.rm = TRUE),
            # Sero15A_count = sum(In_silico_serotype == "15A", na.rm = TRUE),
            # Sero15B_count = sum(In_silico_serotype == "15B", na.rm = TRUE),
            # Sero15BC_count = sum(In_silico_serotype == "15B/15C", na.rm = TRUE),
            # Sero16F_count = sum(In_silico_serotype == "16F", na.rm = TRUE),
            # Sero17F_count = sum(In_silico_serotype == "17F", na.rm = TRUE),
            # Sero18C_count = sum(In_silico_serotype == "18C", na.rm = TRUE),
            # Sero19A_count = sum(In_silico_serotype == "19A", na.rm = TRUE),
            # Sero19F_count = sum(In_silico_serotype == "19F", na.rm = TRUE),
            # Sero22F_count = sum(In_silico_serotype == "22F", na.rm = TRUE),
            # Sero23A_count = sum(In_silico_serotype == "23A", na.rm = TRUE),
            # Sero23F_count = sum(In_silico_serotype == "23F", na.rm = TRUE),
            # Sero3_count = sum(In_silico_serotype == "3", na.rm = TRUE),
            # Sero34_count = sum(In_silico_serotype == "34", na.rm = TRUE),
            # Sero35B_count = sum(In_silico_serotype == "35B", na.rm = TRUE),
            # Sero4_count = sum(In_silico_serotype == "4", na.rm = TRUE),
            # Sero5_count = sum(In_silico_serotype == "5", na.rm = TRUE),
            # Sero6A_count = sum(In_silico_serotype == "6A", na.rm = TRUE),
            # Sero6B_count = sum(In_silico_serotype == "6B", na.rm = TRUE),
            # Sero6C_count = sum(In_silico_serotype == "6C", na.rm = TRUE),
            # Sero7C_count = sum(In_silico_serotype == "7C", na.rm = TRUE),
            # Sero7F_count = sum(In_silico_serotype == "7F", na.rm = TRUE),
            # Sero8_count = sum(In_silico_serotype == "8", na.rm = TRUE),
            # Sero9N_count = sum(In_silico_serotype == "9N", na.rm = TRUE),
            # Sero9V_count = sum(In_silico_serotype == "9V", na.rm = TRUE),
            
            # Counts for  resistance 
            across(all_of(resistance_columns), 
                   list(R_count = ~ sum(. == "R", na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            # Counts for COG groups
            across(all_of(cog_columns),
                   list(count = ~ sum(., na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            .groups = 'drop')

# ggplot(data_summarized_week)+
# geom_line(aes(x=week,y=GPSC3_count))
# ggplot(data_summarized_week)+
#   geom_line(aes(x=week,y=Sero8_count))
# ggplot(data_summarized_week)+
#   geom_line(aes(x=week,y=disease))
##### cleaning weekly epidemiological data
#### complete weekly date range
data_summarized_week$week<-as.Date(data_summarized_week$week)

# Get the full date range across all data
full_date_range_week <- seq.Date(min(data_summarized_week$week), max(data_summarized_week$week), by = "week")
# Get the list of all unique GID_1 and hospital values
all_provinces <- unique(data_summarized_week$GID_1)
all_provinces<-all_provinces[which(!is.na(all_provinces))]
# Create a complete grid of all dates, provinces, and hospitals
complete_grid_week <- expand.grid(week = full_date_range_week, GID_1 = all_provinces)

# Merge this grid with summarized data and fill NAs with 0's
# Create a lookup table for GID_1 and province
province_lookup <- data_summarized_week %>%
  dplyr::select(GID_1, province) %>%
  distinct()
# Add province column to complete_grid_week
complete_grid_week <- complete_grid_week %>%
  left_join(province_lookup, by = "GID_1")
df_complete_week <- complete_grid_week %>%
  left_join(data_summarized_week, by = c("GID_1", "province","week"))
# Identify the columns to replace NAs (all except location/time vars)
non_location_columns <- setdiff(colnames(df_complete_week), c("GID_1", "province", "week"))
replace_list <- as.list(setNames(rep(0, length(non_location_columns)), non_location_columns))
df_complete_week <- df_complete_week %>%
  mutate(across(all_of(non_location_columns), ~replace_na(., 0))) %>%
  ungroup()

# df_complete_week$population<-df_complete_week$population_size

df_complete_week$year <- lubridate::year(df_complete_week$week)
df_complete_week<-left_join(df_complete_week,population_data[,c("GID_1","year","population_size", "population_density")],by=c("year","GID_1"))

## add in a vaccine column to account for
df_complete_week$vaccinated <- ifelse(year(df_complete_week$week)<2009,0,1)


# check population sizes
# ggplot(df_complete_week)+geom_line(aes(x=year,y=population_size,group=province))
# ggplot(df_complete_week)+geom_line(aes(x=week,y=population_size,group=province))

### Calculate the proportions of each COG group as well as the counts
# df_complete_week<-data.frame(df_complete_week)
# cog_names<-grep("COG",colnames(df_complete_week),value=TRUE)
# df_complete_week$total_cog_count <- rowSums(df_complete_week[, cog_names], na.rm = TRUE)
# 
# for (col in cog_names) {
#   prop_col <- sub("_count$", "_prop", col)  # Replace "_count" with "_prop"
#   df_complete_week[[prop_col]] <- df_complete_week[[col]] / df_complete_week$total_cog_count  # Calculate the proportion
# }

# df_complete_week <- df_complete_week %>%
# mutate(epiweek = epiweek(week),               # Extract epidemiological week
#        year = year(week),                     # Extract year
#        weekyear = paste0(epiweek, "-", year))  # Combine epiweek and year as "week-year"

# ggplot(df_complete_week) +
#   geom_line(aes(x=week,y=age_gt80_fem)) + ggplot(df_complete_week) +
#   geom_line(aes(x=week,y=age_gt80_mal))

############MONTHLY DATA##################
## serotypes
sero_df<-data.table(table(data2$In_silico_serotype))
sero_vec<-as.vector(sero_df[which(sero_df$N>35),"V1"])$V1
pcv7vec<-c("23F","6B","14","19F","18C","9V","4")
pcv13vec<-c("3","6A","1","19A","5","7F")
nvts<-sero_vec[sero_vec%notin%c(pcv7vec,pcv13vec)]

## ages
# Check if ageyears is 0 where agemonths is less than 12
age_check <- data2 %>%
  filter(agemonths > 12) %>%
  dplyr::select(agemonths,ageyears)

## summarise
data_summarized_month <- data2 %>%
  mutate(month = floor_date(date, unit = "month")) %>%  # Convert date to start of the month being a Sunday
  group_by(GID_1, province, month) %>%  # Group by GID_1, province, and monthly date
  summarize(disease = sum(disease, na.rm = TRUE),
            sequenced = sum(sequenced, na.rm = TRUE),
            pcv7 = sum(vaccine_status_phen == "PCV7", na.rm = TRUE),  # Count occurrences of PCV7
            pcv13 = sum(vaccine_status_phen == "PCV13", na.rm = TRUE),  # Count occurrences of PCV13
            nvt = sum(vaccine_status_phen == "NVT", na.rm = TRUE),  # Count occurrences of NVT
            unknown_serotype = sum(vaccine_status_phen == "unknown", na.rm = TRUE),  # Count occurrences of unknown serotype
            
            # Counts for Gauteng and Western Cape specifically
            gauteng_count = sum(province == "Gauteng", na.rm=T),
            westerncape_count = sum(province == "Western Cape", na.rm = T),
            
            # Creating counts for ages <1, <5, 5-21,, >65, >80
            age_lt1 = sum(ageyears <1, na.rm=T),
            age_lt6 = sum(ageyears <=5, na.rm=T),
            age_6t14 = sum(ageyears>5 & ageyears<=14,na.rm=T),
            age_15t24 = sum(ageyears>14 & ageyears<=24,na.rm=T),
            age_15t64 = sum(ageyears>14 & ageyears<65,na.rm=T),
            age_18t64 = sum(ageyears>=18 & ageyears<65,na.rm=T),
            age_gt65 = sum(ageyears >65, na.rm=T),
            age_gt80 = sum(ageyears >80, na.rm=T),
            
            # Sex count
            female_count = sum(sex == "FEMALE", na.rm = TRUE),
            male_count = sum(sex == "MALE", na.rm = TRUE),
            unknownsex_count = sum(sex == "UNKNOWN", na.rm = TRUE),
            
            # Age and Sex Count
            age_lt1_fem = sum(ageyears <1 & sex=="FEMALE", na.rm=T),
            age_lt5_fem = sum(ageyears <=5 & sex=="FEMALE", na.rm=T),
            age_gt65_fem = sum(ageyears > 65 & sex=="FEMALE", na.rm=T),
            age_gt80_fem = sum(ageyears > 80 & sex=="FEMALE", na.rm=T),
            age_lt1_mal = sum(ageyears <1 & sex=="MALE", na.rm=T),
            age_lt5_mal = sum(ageyears <=5 & sex=="MALE", na.rm=T),
            age_gt65_mal = sum(ageyears > 65 & sex=="MALE", na.rm=T),
            age_gt80_mal = sum(ageyears > 80 & sex=="MALE", na.rm=T),
            
            # Smoking count
            smoker_count = sum(smoker == "YES", na.rm = TRUE),
            nonsmoker_count = sum(smoker == "NO", na.rm = TRUE),
            unknownsmoker_count = sum(smoker == "", na.rm = TRUE),
            
            # Pregnancy count
            pregnant_count = sum(pregnant == "Y", na.rm = TRUE),
            notpregnant_count = sum(pregnant == "N", na.rm = TRUE),
            unknownpregnant_count = sum(pregnant %notin% c("Y","N"), na.rm = TRUE),
            
            ### mening, bacteremia, other count
            mening_count = sum(spec_diagnosis == "Meningitis", na.rm = TRUE),
            bact_count = sum(spec_diagnosis == "Bacteremia without focus", na.rm = TRUE),
            other_count = sum(spec_diagnosis == "Other", na.rm = TRUE),
            
            # Creating counts for each GPSC type
            GPSC1_count = sum(GPSC == "1", na.rm = TRUE),
            GPSC10_count = sum(GPSC == "10", na.rm = TRUE),
            GPSC13_count = sum(GPSC == "13", na.rm = TRUE),
            GPSC14_count = sum(GPSC == "14", na.rm = TRUE),
            GPSC17_count = sum(GPSC == "17", na.rm = TRUE),
            GPSC2_count = sum(GPSC == "2", na.rm = TRUE),
            GPSC21_count = sum(GPSC == "21", na.rm = TRUE),
            GPSC26_count = sum(GPSC == "26", na.rm = TRUE),
            GPSC3_count = sum(GPSC == "3", na.rm = TRUE),
            GPSC5_count = sum(GPSC == "5", na.rm = TRUE),
            GPSC56_count = sum(GPSC == "56", na.rm = TRUE),
            GPSC41_count = sum(GPSC == "41", na.rm = TRUE),
            GPSC8_count = sum(GPSC == "8", na.rm = TRUE),
            
            # Creating counts for each serotype
            # PCV7
            Sero4_count = sum(serotype == "4", na.rm = TRUE),
            Sero6B_count = sum(serotype == "6B", na.rm = TRUE),
            Sero9V_count = sum(serotype == "9V", na.rm = TRUE),
            Sero14_count = sum(serotype == "14", na.rm = TRUE),
            Sero18C_count = sum(serotype == "18C", na.rm = TRUE),
            Sero19F_count = sum(serotype == "19F", na.rm = TRUE),
            Sero23F_count = sum(serotype == "23F", na.rm = TRUE),
            # PCV13
            Sero1_count = sum(serotype == "1", na.rm = TRUE),
            Sero3_count = sum(serotype == "3", na.rm = TRUE),
            # Sero5_count = sum(serotype == "5", na.rm = TRUE), ## none in this dataset
            Sero6A_count = sum(serotype == "6A", na.rm = TRUE),
            Sero7F_count = sum(serotype == "7F", na.rm = TRUE),
            Sero19A_count = sum(serotype == "19A", na.rm = TRUE),
            # NVT
            Sero12F_count = sum(serotype == "12F", na.rm = TRUE),
            Sero15A_count = sum(serotype == "15A", na.rm = TRUE),
            Sero15BC_count = sum(serotype == "15B/C", na.rm = TRUE), 
            Sero16F_count = sum(serotype == "16F", na.rm = TRUE),
            Sero8_count = sum(serotype == "8", na.rm = TRUE),
            Sero9N_count = sum(serotype == "9N", na.rm = TRUE),
            
            # Counts for  resistance 
            across(all_of(resistance_columns), 
                   list(R_count = ~ sum(. == "R", na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            # Counts for COG groups
            across(all_of(cog_columns), 
                   list(count = ~ sum(., na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            .groups = 'drop') 

## include the proportion of each GPSC per the number sequenced each month
data_summarized_month <- data_summarized_month %>%
  mutate(year = year(month)) %>%  
  group_by(year, month) %>%
  summarise(month_seq = sum(sequenced), .groups = 'drop') %>%
  left_join(data_summarized_month, by = c("month")) %>%  # Merge back to get month_seq
  mutate(across(starts_with("GPSC"), ~ . / month_seq, .names = "{.col}_prop")) %>%
  rename_with(~ gsub("_count_prop$", "_prop", .), ends_with("_count_prop")) # Calculate proportions


data_summarized_month$month<-as.Date(data_summarized_month$month)

# Get the full date range across all data
full_date_range_month <- seq.Date(min(data_summarized_month$month), max(data_summarized_month$month), by = "month")
# Get the list of all unique GID_1 and hospital values
all_provinces <- unique(data_summarized_month$GID_1)
all_provinces<-all_provinces[which(!is.na(all_provinces))]

# Create a complete grid of all dates, provinces, and hospitals
complete_grid_month <- expand.grid(month = full_date_range_month, GID_1 = all_provinces)
gpsc_cols<-colnames(data_summarized_month)[grep("GPSC",colnames(data_summarized_month))]
sero_cols<-colnames(data_summarized_month)[grep("Sero",colnames(data_summarized_month))]

# Merge this copmlete grid with your summarized data
province_lookup <- data_summarized_month %>%
  dplyr::select(GID_1, province) %>%
  distinct()
# Add province column to complete_grid_week
complete_grid_month <- complete_grid_month %>%
  left_join(province_lookup, by = "GID_1")
df_complete_month <- complete_grid_month %>%
  left_join(data_summarized_month, by = c("GID_1", "province","month"))
# Identify the columns to replace NAs (all except location/time vars)
non_location_columns <- setdiff(colnames(df_complete_month), c("GID_1", "province", "month"))
replace_list <- as.list(setNames(rep(0, length(non_location_columns)), non_location_columns))
df_complete_month <- df_complete_month %>%
  mutate(across(all_of(non_location_columns), ~replace_na(., 0))) %>%
  ungroup()

## set month and year and join
df_complete_month$year <- as.integer(lubridate::year(df_complete_month$month))
df_complete_month<-left_join(df_complete_month,population_data[,c("year","GID_1","population_size", "population_density")],by=c("year", "GID_1"))
## add in a vaccine column to account for
df_complete_month$vaccinated <- ifelse(year(df_complete_month$month)<2009,0,1)

### check how data looks
## check population sizes
# ggplot(population_data)+geom_line(aes(x=year,y=population_size,group=GID_1))
# ggplot(df_complete_month)+geom_line(aes(x=month,y=population_size,group=province))
# # Convert data to long format
# data_long <- df_complete_month %>%
#   dplyr::select(month, starts_with("GPSC"), -matches("_count$")) %>%  # Select only proportion columns
#   pivot_longer(cols = starts_with("GPSC"),
#                names_to = "GPSC_type",
#                values_to = "proportion")
# # Plot
# ggplot(data_long, aes(x = month, y = proportion, fill = GPSC_type, group = GPSC_type)) +
#   geom_col(position = "fill", alpha = 0.8) +  # Stacked area chart
#   scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") + # Format x-axis
#   labs(title = "GPSC Proportions Over Time",
#        x = "Month",
#        y = "Proportion",
#        color = "GPSC Type") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
################################################################################
###### SAVE DATA FRAMES AS AGGREGATED AT WEEKLY, MONTHLY, OR LEFT AS POINTS ######
################################################################################  
## remove before 2005 due to climate data going from then
## drop everything before 2005 because the climatic data only starts in 2005
df_complete_week <- subset(df_complete_week,year(df_complete_week$week)>2004)
df_complete_month <- subset(df_complete_month,year(df_complete_month$month)>2004)
df_complete_day <- subset(df_complete_day,year(df_complete_day$date)>2004)

## set area id province, year, and month as a factor
# df_complete$id_u<-as.numeric(as.factor(df_complete$GID_1))


## plot cases
d<-ggplot()+
  geom_line(data=df_complete_day,aes(x=date,y=disease,group=GID_1),color="purple")
w<-ggplot()+
  geom_line(data=df_complete_week,aes(x=week,y=disease,group=GID_1),color="blue")

m<-ggplot()+
  geom_line(data=df_complete_month,aes(x=month,y=disease,group=GID_1),color="red")
library(patchwork)
d+w+m

# Ensure year values match key types before applying mapping
## add factor indices for INLA
df_complete_day <- re_id2(df_complete_day, "day")
df_complete_week <- re_id2(df_complete_week, "week")
df_complete_month <- re_id2(df_complete_month, "month")


################################################################################
######  SOCIODEMOGRAPHIC DATA MERGE ######
################################################################################  
# df_complete_week<-fread("/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_weekly_province.csv")
# df_complete_month<-fread("/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_monthly_province.csv")

sociodemographic_df <- readRDS(file="/home/sbelman/Documents/env_sa_manuscript/input_datasets/sociodemographic/statsZA_sociodemographic_data.rds")
sociodemograph_sub <- sociodemographic_df[,c("year","month","NAME_1","GID_1","flu_positivity","rsv_positivity",
                                             "Incidence_in_youth_15_24", "Total_people_living_with_HIV","ART_coverage","Total_on_ART", "Dependency_ratio",
                                             "Aging_index","percent_children_under15","percent_elderly","in_migration_5year","net_migration_5year","out_migration_5year",
                                             "PCV_coverage_3dose","WUENIC_coverage_3dose" ,"hiv_prevalence_adults_15_49","hiv_prevalence_youth_15_24","hiv_prevalence_totalpop")]


df_complete_week <- left_join(df_complete_week,sociodemograph_sub,by=c("GID_1","year","month"))
df_complete_month <- left_join(df_complete_month,sociodemograph_sub,by=c("GID_1","year","month"))
df_complete_day <- left_join(df_complete_day,sociodemograph_sub,by=c("GID_1","year","month"))

# ## save as base data frame
write.table(df_complete_day, file="/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_daily_province.csv",quote=FALSE,row.names = FALSE, col.names = TRUE,sep=",")
write.table(df_complete_week, file="/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_weekly_province.csv",quote=FALSE,row.names = FALSE, col.names = TRUE,sep=",")
write.table(df_complete_month, file="/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_monthly_province.csv",quote=FALSE,row.names = FALSE, col.names = TRUE,sep=",")

################################################################################
###### CLIMATE DATA MERGE  ######
################################################################################   

################DAILY DATA MERGE ######################
## load climatic and humidity and merge by the day
clim_adm1_daily<-fread("/home/sbelman/Documents/BRD/SouthAfrica/climate/data/climate_adm1_daily.csv")
hum_adm1_daily<-fread("/home/sbelman/Documents/BRD/SouthAfrica/climate/data/hum_adm1_daily.csv")
aq_adm1_daily<-fread("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/reanalysis/aq_adm1_daily_cams.csv")
# aq_adm1_daily<-fread("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/reanalysis/aq_adm1_daily_merra.csv")
sfcWind_adm1_daily <- fread("/home/sbelman/Documents/BRD/SouthAfrica/climate/data/sfcWind_adm1_daily.csv")

# sort dates
clim_adm1_daily$date <- as.Date(clim_adm1_daily$date)
hum_adm1_daily$date <- as.Date(hum_adm1_daily$date)
aq_adm1_daily$date <- as.Date(aq_adm1_daily$date)
sfcWind_adm1_daily$date <- as.Date(sfcWind_adm1_daily$date)

####### check  data ##########
# clim_adm1_daily$date <- as.Date(clim_adm1_daily$date)
# test <- clim_adm1_daily %>%
#   group_by(date) %>%
#   summarise(prlr = mean(prlr,na.rm=T))
# ggplot(test)+
#   geom_line(aes(x=date,y=prlr))
#############

aq_adm1_daily$date<-as.Date(aq_adm1_daily$date)
## adjust the aq data to ug/m3
aq_adm1_daily$pm2p5 <- aq_adm1_daily$pm2p5 * 1e9
## pm10 for cams from kg m-3  to ug/m3
aq_adm1_daily$pm10 <- aq_adm1_daily$pm10 * 1e9
## pm10 for cams from kg m-1  to ug/m3
air_density_reg <- 1.225
air_density_so2 <- 2.63  # SO2 density in kg/m³ (standard at sea level, 15°C)
air_density_o3 <- 2.14  # O3 density in kg/m³ (standard at sea level, 15°C)
aq_adm1_daily$o3 <- aq_adm1_daily$o3 * air_density_o3 * 1e9
aq_adm1_daily$so2 <- aq_adm1_daily$so2 * air_density_so2 * 1e9

### convert anything >100ug/m3 t0 100ug/m3
aq_adm1_daily$o3[which(aq_adm1_daily$o3>400)] <- 400
aq_adm1_daily$so2[which(aq_adm1_daily$so2>400)] <- 400
aq_adm1_daily$pm10[which(aq_adm1_daily$pm10>200)] <- 200
aq_adm1_daily$pm2p5[which(aq_adm1_daily$pm2p5>200)] <-200

## merge humiidity and other meteorology data
hum_adm1_daily$date <- as.Date(hum_adm1_daily$date)
tmp<-hum_adm1_daily[,c("date","adm1_id","hurs","absh")]
rm(hum_adm1_daily)
clim_hum<-left_join(clim_adm1_daily,tmp,by=c("date","adm1_id"))
rm(tmp)

### merge humiidity and climate with air quality
clim_hum$date<-as.Date(clim_hum$date)
tmpaq<-aq_adm1_daily[,c("date","adm1_id","pm2p5","pm10","o3","so2")]
clim_hum_aq<-left_join(clim_hum,tmpaq,by=c("date","adm1_id"))
rm(tmpaq)

## merge wind with the rest
sfcWind_adm1_daily$date<-as.Date(sfcWind_adm1_daily$date)
clim_hum_aq_wind <- left_join(clim_hum_aq,sfcWind_adm1_daily[,c("date","adm1_id","adm1_name","sfcWind")], by=c("date","adm1_id","adm1_name"))
rm(clim_hum_aq)

## load spi drought indicators monthly
spi_adm1_monthly<-fread("/home/sbelman/Documents/BRD/SouthAfrica/climate/drought/adm1_droughtindices_monthly.csv")
spi_adm1_monthly<- data.frame(spi_adm1_monthly)
## remove infinities## ## 
spi_vec<-c(paste0(rep("spi",4),c(1,3,6,12)),paste0(rep("spei",4),c(1,3,6,12)))
adm_vec<-unique(spi_adm1_monthly$adm1_id)
for(dd in 1:length(spi_vec)){
  spi_idx<-as.numeric(which(colnames(spi_adm1_monthly)%in%spi_vec[dd]))
  for(rr in 1:length(adm_vec)){
    for(mm in 1:12){	
      maxspi1 <- max(spi_adm1_monthly[which(spi_adm1_monthly$month == mm & 
                                              spi_adm1_monthly$adm1_id == adm_vec[rr] & 
                                              is.finite(spi_adm1_monthly[,spi_idx])),spi_idx])
      minspi1 <- min(spi_adm1_monthly[which(spi_adm1_monthly$month == mm & 
                                              spi_adm1_monthly$adm1_id == adm_vec[rr] & 
                                              is.finite(spi_adm1_monthly[,spi_idx])),spi_idx])
      spi_adm1_monthly[which(spi_adm1_monthly$month == mm & 
                               spi_adm1_monthly$adm1_id == adm_vec[rr] & 
                               spi_adm1_monthly[,spi_idx] == Inf),spi_idx] <- maxspi1
      spi_adm1_monthly[which(spi_adm1_monthly$month == mm & 
                               spi_adm1_monthly$adm1_id == adm_vec[rr] & 
                               spi_adm1_monthly[,spi_idx] == -Inf),spi_idx] <- minspi1
    }
  }
  print(table(is.infinite(spi_adm1_monthly[,spi_idx])))
}
## ## ## ## ## 
spi_adm1_monthly<-data.table(spi_adm1_monthly)
tmp<-spi_adm1_monthly[,-c("adm1_name")]
climate_data <- left_join(clim_hum_aq_wind,tmp,by=c("adm1_id","month","year"))
## change the name so it matches
colnames(climate_data)[which(colnames(climate_data)=="adm1_id")]<-"GID_1"
# save
write.table(climate_data,file="/home/sbelman/Documents/env_sa_manuscript/input_datasets/climate/env_adm1_daily.csv", quote=FALSE , row.names = FALSE, col.names = TRUE,sep=",")

##########WEEKLY DATA MERGE ############
## load climatic and humidity and merge by the day
clim_adm1_weekly2<-fread("/home/sbelman/Documents/BRD/SouthAfrica/climate/data/climate_adm1_weekly.csv")
sfcWind_adm1_weekly <-fread("/home/sbelman/Documents/BRD/SouthAfrica/climate/data/sfcWind_adm1_weekly.csv")
clim_adm1_weekly <- left_join(clim_adm1_weekly2,sfcWind_adm1_weekly,by=c("adm1_id", "adm1_name","weekyear"))

## split the epiweek column to create a date column with the week
clim_adm1_weekly[,c("week_num","year") := tstrsplit(weekyear, "-", fixed = TRUE)]
clim_adm1_weekly$week_num <- as.numeric(clim_adm1_weekly$week_num)
clim_adm1_weekly$year <- as.numeric(clim_adm1_weekly$year)

## weekly 
clim_adm1_weekly[, week := floor_date(make_date(year) + weeks(week_num - 1), "week")]

# # Check for duplicates
duplicates <- clim_adm1_weekly %>%
  group_by(adm1_id, week) #%>%
# filter(n() > 1)
# clim_adm1_weekly[which(clim_adm1_weekly$week=="2008-12-28" & clim_adm1_weekly$adm1_id=="ZAF.1.1_1"),]

### take the means of the duplicate weeks
clim_adm1_weekly <- clim_adm1_weekly %>%
  group_by(week,adm1_id, adm1_name) %>%
  summarise(
    tas = mean(tas, na.rm = TRUE),  # Mean of tas for duplicates
    tasmin = mean(tasmin, na.rm = TRUE),   # Mean of tasmin for duplicate
    tasmax = mean(tasmax, na.rm = TRUE),   # Mean of tasmax for duplicate
    prlrmean = mean(prlr_mean, na.rm = TRUE),   # Mean of prlr for duplicate
    prlrsum = sum(prlr_sum, na.rm = TRUE),   # Mean of prlr for duplicate
    prlrmax= mean(prlr_max, na.rm = TRUE),   # Mean of prlr for duplicate
    sfcWind = mean(sfcWind, na.rm = TRUE)
  ) %>%
  ungroup()
clim_adm1_weekly<-unique(clim_adm1_weekly[,c("adm1_id","adm1_name","tas","tasmin","tasmax","prlrsum","prlrmax","prlrmean","sfcWind","week")])

####humidity
hum_adm1_weekly<-fread("/home/sbelman/Documents/BRD/SouthAfrica/climate/data/hum_adm1_weekly.csv")
## split the epiweek column to create a date column with the week
hum_adm1_weekly[,c("week_num","year") := tstrsplit(epiweek, "-", fixed = TRUE)]
hum_adm1_weekly$week_num <- as.numeric(hum_adm1_weekly$week_num)
hum_adm1_weekly$year <- as.numeric(hum_adm1_weekly$year)
hum_adm1_weekly[, week := floor_date(make_date(year) + weeks(week_num - 1), "week")]
# # Check for duplicates
# duplicates <- hum_adm1_weekly %>%
#   group_by(adm1_id, week) %>%
#   filter(n() > 1)
# 
### take the means of the duplicate weeks
hum_adm1_weekly <- hum_adm1_weekly %>%
  group_by(week, adm1_id) %>%
  mutate(
    absh = mean(absh, na.rm = TRUE),  # Mean of absh for duplicates
    hurs = mean(hurs, na.rm = TRUE)   # Mean of hurs for duplicate
  ) %>%
  ungroup()
# hum_adm1_weekly[which(hum_adm1_weekly$week=="2008-12-28" & hum_adm1_weekly$adm1_id=="ZAF.1.1_1"),]
hum_adm1_weekly<-unique(hum_adm1_weekly[,c("week","adm1_id","hurs","absh")])
hum_adm1_weekly<-hum_adm1_weekly[,c("week","adm1_id","hurs","absh")]
clim_hum<-left_join(clim_adm1_weekly,hum_adm1_weekly,by=c("week","adm1_id"))
clim_hum$week<-as.Date(clim_hum$week)

### air quality
aq_adm1_weekly<-fread("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/reanalysis/aq_adm1_weekly_cams.csv")
## ## adjust the aq data to ug/m3
aq_adm1_weekly$pm2p5 <- aq_adm1_weekly$pm2p5 * 1e9
## pm10 for cams from kg m-3  to ug/m3
aq_adm1_weekly$pm10 <- aq_adm1_weekly$pm10 * 1e9
## pm10 for cams from kg m-1  to ug/m3
air_density_reg <- 1.225
air_density_so2 <- 2.63  # SO2 density in kg/m³ (standard at sea level, 15°C)
air_density_o3 <- 2.14  # O3 density in kg/m³ (standard at sea level, 15°C)
aq_adm1_weekly$o3 <- aq_adm1_weekly$o3 * air_density_o3 * 1e9
aq_adm1_weekly$so2 <- aq_adm1_weekly$so2 * air_density_so2 * 1e9
### convert anything >100ug/m3 t0 100ug/m3
aq_adm1_weekly$o3[which(aq_adm1_weekly$o3>400)] <- 400
aq_adm1_weekly$so2[which(aq_adm1_weekly$so2>400)] <- 400
aq_adm1_weekly$pm10[which(aq_adm1_weekly$pm10>200)] <- 200
aq_adm1_weekly$pm2p5[which(aq_adm1_weekly$pm2p5>200)] <-200



## split the epiweek column to create a date column with the week
aq_adm1_weekly[,c("week_num","year") := tstrsplit(epiweek, "-", fixed = TRUE)]
aq_adm1_weekly$week_num <- as.numeric(aq_adm1_weekly$week_num)
aq_adm1_weekly$year <- as.numeric(aq_adm1_weekly$year)
aq_adm1_weekly[, week := floor_date(make_date(year) + weeks(week_num - 1), "week")]
# # Check for duplicates
duplicates <- aq_adm1_weekly %>%
  group_by(adm1_id, week) %>%
  filter(n() > 1)

### take the means of the duplicate weeks
aq_adm1_weekly <- aq_adm1_weekly %>%
  group_by(week, adm1_id) %>%
  mutate(
    pm2p5 = mean(pm2p5, na.rm = TRUE),  # Mean of tas for duplicates
    pm10 = mean(pm10, na.rm = TRUE),   # Mean of tasmin for duplicate
    o3 = mean(o3, na.rm = TRUE),   # Mean of prlr for duplicate
    so2 = mean(so2, na.rm = TRUE)   # Mean of prlr for duplicate
    
  ) %>%
  ungroup()
aq_adm1_weekly<-unique(aq_adm1_weekly[,c("adm1_id","adm1_name","pm2p5","pm10","o3","so2","week")])


aq_adm1_weekly$week<-as.Date(aq_adm1_weekly$week)
tmpaq<-aq_adm1_weekly[,c("week","adm1_id","pm2p5","pm10","o3","so2")]
clim_hum_aq<-left_join(clim_hum,tmpaq,by=c("week","adm1_id"))
## load spi drought indicators monthly
spi_adm1_monthly<-fread("/home/sbelman/Documents/BRD/SouthAfrica/climate/drought/adm1_droughtindices_monthly.csv")
spi_adm1_monthly<- data.frame(spi_adm1_monthly)
## remove infinities## ## 
spi_vec<-c(paste0(rep("spi",4),c(1,3,6,12)),paste0(rep("spei",4),c(1,3,6,12)))
adm_vec<-unique(spi_adm1_monthly$adm1_id)
for(dd in 1:length(spi_vec)){
  spi_idx<-as.numeric(which(colnames(spi_adm1_monthly)%in%spi_vec[dd]))
  for(rr in 1:length(adm_vec)){
    for(mm in 1:12){	
      maxspi1 <- max(spi_adm1_monthly[which(spi_adm1_monthly$month == mm & 
                                              spi_adm1_monthly$adm1_id == adm_vec[rr] & 
                                              is.finite(spi_adm1_monthly[,spi_idx])),spi_idx])
      minspi1 <- min(spi_adm1_monthly[which(spi_adm1_monthly$month == mm & 
                                              spi_adm1_monthly$adm1_id == adm_vec[rr] & 
                                              is.finite(spi_adm1_monthly[,spi_idx])),spi_idx])
      spi_adm1_monthly[which(spi_adm1_monthly$month == mm & 
                               spi_adm1_monthly$adm1_id == adm_vec[rr] & 
                               spi_adm1_monthly[,spi_idx] == Inf),spi_idx] <- maxspi1
      spi_adm1_monthly[which(spi_adm1_monthly$month == mm & 
                               spi_adm1_monthly$adm1_id == adm_vec[rr] & 
                               spi_adm1_monthly[,spi_idx] == -Inf),spi_idx] <- minspi1
    }
  }
  print(table(is.infinite(spi_adm1_monthly[,spi_idx])))
}
## ## ## ## ## 
spi_adm1_monthly<-data.table(spi_adm1_monthly)
tmp<-spi_adm1_monthly[,-c("adm1_name")]
# add a month column
clim_hum_aq$month <- month(clim_hum_aq$week)
clim_hum_aq$year <- year(clim_hum_aq$week)

climate_data_weekly <- left_join(clim_hum_aq,tmp,by=c("adm1_id","month","year"))
## change the name so it matches
colnames(climate_data_weekly)[which(colnames(climate_data_weekly)=="adm1_id")]<-"GID_1"

# save
write.table(climate_data_weekly,file="/home/sbelman/Documents/env_sa_manuscript/input_datasets/climate/env_adm1_weekly.csv", quote=FALSE , row.names = FALSE, col.names = TRUE,sep=",")

######## MONTHLY DATA MERGE ####################
## load climatic and humidity and merge by the day
clim_adm1_daily<-fread("/home/sbelman/Documents/BRD/SouthAfrica/climate/data/climate_adm1_daily.csv")
hum_adm1_daily<-fread("/home/sbelman/Documents/BRD/SouthAfrica/climate/data/hum_adm1_daily.csv")
aq_adm1_daily<-fread("/home/sbelman/Documents/BRD/SouthAfrica/airquality/data/reanalysis/aq_adm1_daily_cams.csv")
sfcWind_adm1_daily <- fread("/home/sbelman/Documents/BRD/SouthAfrica/climate/data/sfcWind_adm1_daily.csv")

# sort dates
clim_adm1_daily$date <- as.Date(clim_adm1_daily$date)
hum_adm1_daily$date <- as.Date(hum_adm1_daily$date)
aq_adm1_daily$date <- as.Date(aq_adm1_daily$date)
sfcWind_adm1_daily$date <- as.Date(sfcWind_adm1_daily$date)

## adjust the aq data to ug/m3
aq_adm1_daily$pm2p5 <- aq_adm1_daily$pm2p5 * 1e9
## pm10 for cams from kg m-3  to ug/m3
aq_adm1_daily$pm10 <- aq_adm1_daily$pm10 * 1e9
## pm10 for cams from kg m-3  to ug/m3
air_density_reg <- 1.225
air_density_so2 <- 2.63  # SO2 density in kg/m³ (standard at sea level, 15°C)
air_density_o3 <- 2.14  # O3 density in kg/m³ (standard at sea level, 15°C)
aq_adm1_daily$o3 <- aq_adm1_daily$o3 * air_density_o3 * 1e9
aq_adm1_daily$so2 <- aq_adm1_daily$so2 * air_density_so2 * 1e9

### convert anything >100ug/m3 t0 100ug/m3
aq_adm1_daily$o3[which(aq_adm1_daily$o3>400)] <- 400
aq_adm1_daily$so2[which(aq_adm1_daily$so2>400)] <- 400
aq_adm1_daily$pm10[which(aq_adm1_daily$pm10>200)] <- 200
aq_adm1_daily$pm2p5[which(aq_adm1_daily$pm2p5>200)] <-200

### merge all the data sets
tmp<-hum_adm1_daily[,c("date","adm1_id","hurs","absh")]
clim_hum<-left_join(clim_adm1_daily,tmp,by=c("date","adm1_id"))
clim_hum$date<-as.Date(clim_hum$date)
tmpaq<-aq_adm1_daily[,c("date","adm1_id","pm2p5","pm10","o3","so2")]
clim_hum_aq<-left_join(clim_hum,tmpaq,by=c("date","adm1_id"))

sfcWind_adm1_daily$date <- as.Date(sfcWind_adm1_daily$date)
clim_hum_aq_wind<-left_join(clim_hum_aq,sfcWind_adm1_daily,by=c("date","adm1_id","adm1_name"))
clim_hum_aq_wind<-clim_hum_aq_wind[,c("date","adm1_id","adm1_name","tas","tasmin","tasmax","prlr","hurs","absh","pm2p5","pm10","o3","so2", "sfcWind")]

### aggregate them monthly
clim_hum_aq_monthly<-clim_hum_aq_wind %>%
  mutate(year_month = floor_date(date, unit = "month")) %>%  # Convert date to start of the month 
  group_by(adm1_id, adm1_name, year_month) %>%  # Group by GID_1, province, and weekly date
  summarise(
    tas = mean(tas, na.rm = TRUE),  
    tasmin = mean(tasmin, na.rm = TRUE), 
    tasmax = mean(tasmax, na.rm = TRUE), 
    prlrmean = mean(prlr, na.rm = TRUE), 
    prlrsum = sum(prlr, na.rm = TRUE), 
    prlrmax = max(prlr, na.rm = TRUE), 
    sfcWind = mean(sfcWind, na.rm = TRUE), 
    hurs = mean(hurs, na.rm = TRUE), 
    absh = mean(absh, na.rm = TRUE),   
    pm2p5 = mean(pm2p5, na.rm = TRUE),  
    pm10 = mean(pm10, na.rm = TRUE),
    so2 = mean(so2, na.rm = TRUE),
    o3 = mean(o3, na.rm = TRUE)
    
  ) %>%
  ungroup()        
clim_hum_aq_monthly$month<-month(clim_hum_aq_monthly$year_month)
clim_hum_aq_monthly$year<-year(clim_hum_aq_monthly$year_month)

# ggplot(clim_hum_aq_monthly)+
# geom_line(aes(x=year_month,y=pm10,group=adm1_id))
## load spi drought indicators monthly
spi_adm1_monthly<-fread("/home/sbelman/Documents/BRD/SouthAfrica/climate/drought/adm1_droughtindices_monthly.csv")
spi_adm1_monthly<- data.frame(spi_adm1_monthly)
## remove infinities## ## 
spi_vec<-c(paste0(rep("spi",4),c(1,3,6,12)),paste0(rep("spei",4),c(1,3,6,12)))
adm_vec<-unique(spi_adm1_monthly$adm1_id)
for(dd in 1:length(spi_vec)){
  spi_idx<-as.numeric(which(colnames(spi_adm1_monthly)%in%spi_vec[dd]))
  for(rr in 1:length(adm_vec)){
    for(mm in 1:12){	
      maxspi1 <- max(spi_adm1_monthly[which(spi_adm1_monthly$month == mm & 
                                              spi_adm1_monthly$adm1_id == adm_vec[rr] & 
                                              is.finite(spi_adm1_monthly[,spi_idx])),spi_idx])
      minspi1 <- min(spi_adm1_monthly[which(spi_adm1_monthly$month == mm & 
                                              spi_adm1_monthly$adm1_id == adm_vec[rr] & 
                                              is.finite(spi_adm1_monthly[,spi_idx])),spi_idx])
      spi_adm1_monthly[which(spi_adm1_monthly$month == mm & 
                               spi_adm1_monthly$adm1_id == adm_vec[rr] & 
                               spi_adm1_monthly[,spi_idx] == Inf),spi_idx] <- maxspi1
      spi_adm1_monthly[which(spi_adm1_monthly$month == mm & 
                               spi_adm1_monthly$adm1_id == adm_vec[rr] & 
                               spi_adm1_monthly[,spi_idx] == -Inf),spi_idx] <- minspi1
    }
  }
  print(table(is.infinite(spi_adm1_monthly[,spi_idx])))
}
## ## ## ## ## 
spi_adm1_monthly<-data.table(spi_adm1_monthly)
tmp<-spi_adm1_monthly[,-c("adm1_name")]
climate_data_monthly <- left_join(clim_hum_aq_monthly,tmp,by=c("adm1_id","month","year"))
## change the name so it matches
colnames(climate_data_monthly)[which(colnames(climate_data_monthly)=="adm1_id")]<-"GID_1"
write.table(climate_data_monthly,file="/home/sbelman/Documents/env_sa_manuscript/input_datasets/climate/env_adm1_monthly.csv", quote=FALSE , row.names = FALSE, col.names = TRUE,sep=",")


################################################################################
### JOIN CLIMATE AND DISEASE DATA ######
################################################################################  

#########CHECK SHAPE FILE INDEX ORDER #########
shp<-st_read("/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/gadm41_namematch_ZAF_1.shp")
adm_idx <- shp[,c("GID_1","NAME_1")]
adm_idx <- st_drop_geometry(adm_idx)
idx <- unique(shp$GID_1)
## save adjacency matrix
nb <- poly2nb(shp_adm1)
nb2INLA(file="/home/sbelman/Documents/env_sa_manuscript/input_datasets/shps/sa_adjacency_map_adm1.adj", nb)

seq_month <- fread("/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/number_sequenced_monthly.csv")

########DAILY DATA#############
disease_data<-fread("/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_daily_province.csv")
climate_data<-fread("/home/sbelman/Documents/env_sa_manuscript/input_datasets/climate/env_adm1_daily.csv")
# ## merge climate and disease data
sa_adm1_daily<-left_join(climate_data,disease_data,by=c("date","year","month","GID_1"))
# 
# ## reorder for idx
sa_adm1_daily <- sa_adm1_daily %>%
  mutate(GID_1 = factor(GID_1, levels = idx)) %>%
  arrange(date, GID_1)
# 
# 
# ######### CHECKS BEFORE SAVE
# ## check SPATIAL match  - if TRUE they match
###### reorder so that the indexes all match the order
sa_adm1_daily$id_u <- as.numeric(as.factor(sa_adm1_daily$GID_1))
all(diff(unique(sa_adm1_daily$id_u)) == 1)
if(all(diff(match(unique(sa_adm1_daily$GID_1),adm_idx$GID_1))==1)){print("the indexes match the order")}

## check YEAR match  - if TRUE they match
table_y <- table(sa_adm1_daily$year, sa_adm1_daily$id_y)
table_m <- table(sa_adm1_daily$month, sa_adm1_daily$id_m)
if(all(table_m[row(table_m) != col(table_m)] == 0)){print("All off diagonals are zero for the months")}
if(all(table_y[row(table_y) != col(table_y)] == 0)){print("All off diagonals are zero for the years")}
# 
# ## add monthly sequenced
sa_adm1_daily <- left_join(sa_adm1_daily, seq_month, by= c("year","month"))
# 
# ######## SAVE
write.table(sa_adm1_daily,file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_daily.csv", quote=FALSE , row.names = FALSE, col.names = TRUE,sep=",")
# 
# 
########WEEKLY DATA############
disease_data<-fread("/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_weekly_province.csv")
climate_data<-fread("/home/sbelman/Documents/env_sa_manuscript/input_datasets/climate/env_adm1_weekly.csv")
## merge climate and disease data
sa_adm1_weekly<-left_join(climate_data,disease_data,by=c("week","month","GID_1","year"))


sa_adm1_weekly <- subset(sa_adm1_weekly,sa_adm1_weekly$year>2004)
if("date"%notin%colnames(sa_adm1_weekly)){
  print("Change week to date")
  df<- sa_adm1_weekly
  df <- df %>%
    mutate(date = as.Date(week))}else{
      df<-sa_adm1_weekly
    }



## reorder for idx
df <- df %>%
  mutate(GID_1 = factor(GID_1, levels = idx)) %>%
  arrange(week, GID_1)



######### CHECKS BEFORE SAVE
## check SPATIAL match  - if TRUE they match
###### reorder so that the indexes all match the order
df$id_u <- as.numeric(as.factor(df$GID_1))
all(diff(unique(df$id_u)) == 1)
if(all(diff(match(unique(df$GID_1),adm_idx$GID_1))==1)){print("the indexes match the order")}

## check YEAR match  - if TRUE they match
table_y <- table(df$year, df$id_y)
table_m <- table(df$month, df$id_m)
if(all(table_m[row(table_m) != col(table_m)] == 0)){print("All off diagonals are zero for the months")}
if(all(table_y[row(table_y) != col(table_y)] == 0)){print("All off diagonals are zero for the years")}

## add monthly sequenced
df2 <- left_join(df, seq_month, by= c("year","month"))

########save
write.table(df2,file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly.csv", quote=FALSE , row.names = FALSE, col.names = TRUE,sep=",")

########MONTHLY DATA#############
disease_data<-fread("/home/sbelman/Documents/env_sa_manuscript/input_datasets/disease/SA_disease_monthly_province.csv")
climate_data<-fread("/home/sbelman/Documents/env_sa_manuscript/input_datasets/climate/env_adm1_monthly.csv")
colnames(climate_data)[colnames(climate_data)=="year_month"]<-"date"
table(climate_data$date%in%disease_data$date)
## merge climate and disease data
sa_adm1_monthly<-left_join(climate_data,disease_data,by=c("date","GID_1","year","month"))
if("date"%notin%colnames(sa_adm1_monthly)){
  print("Change monthh to date")
  df <- sa_adm1_monthly
  df <- df %>%
    mutate(date = as.Date(month)) }else{
      df <- sa_adm1_monthly } 

## reorder for idx
df <- df %>%
  mutate(GID_1 = factor(GID_1, levels = idx)) %>%
  arrange(date, GID_1)

######### CHECKS BEFORE SAVE
## check SPATIAL match  - if TRUE they match
###### reorder so that the indexes all match the order
df$id_u <- as.numeric(as.factor(df$GID_1))
all(diff(unique(df$id_u)) == 1)
if(all(diff(match(unique(df$GID_1),adm_idx$GID_1))==1)){print("the indexes match the order")}

## check YEAR match  - if TRUE they match
table_y <- table(df$year, df$id_y)
table_m <- table(df$month, df$id_m)
if(all(table_m[row(table_m) != col(table_m)] == 0)){print("All off diagonals are zero for the months")}
if(all(table_y[row(table_y) != col(table_y)] == 0)){print("All off diagonals are zero for the years")}

########save
write.table(df,file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_monthly.csv", quote=FALSE , row.names = FALSE, col.names = TRUE,sep=",")



################################################################################
### TEST DISTRIBUTION OF VARIABLES AND SCALE IF NON NORMAL ######
################################################################################ 

sa_adm1_weekly<-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly.csv")
all_names<- colnames(sa_adm1_weekly)
cov_names<-c(grep("tas",all_names,value=TRUE),grep("sp",all_names,value=TRUE),"hurs","absh","prlrsum","prlrmax","prlrmean","pm2p5","pm10", "o3", "so2","sfcWind")

# test distribution before scaling
# test_distribution<-function(data,cov){
#   p <- data %>%
#     ggplot()+
#       geom_point(aes(x=.data[[cov]],y=disease))+
#       ggtitle(cov)+
#       theme_bw()
#   return(p)
# }
# dist_list<-lapply(cov_names, function(x) test_distribution(sa_adm1_weekly,x))
# grid.arrange(grobs=dist_list)

## test normality
# result <- t(sapply(cov_names, function(cov) c(p_value = as.numeric(shapiro.test(sample(sa_adm1_weekly[[cov]], 5000))$p.value), var_name = cov)))

# st<-shapiro.test(sample(sa_adm1_weekly$spi3,5000))
# st<-shapiro.test(sample(sa_adm1_weekly$pm2p5,5000))


## scale prlr and absh
sa_adm1_weekly_sc<-sa_adm1_weekly
scale_covs<-c("hurs","absh","prlrsum","prlrmax","prlrmean" )
# Scale the specified columns using .SD and lapply
sa_adm1_weekly_sc[, (scale_covs) := lapply(.SD, scale), .SDcols = scale_covs]
write.table(sa_adm1_weekly_sc,file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly_sc.csv", quote=FALSE , row.names = FALSE, col.names = TRUE,sep=",")

##################### MONTHLY SCALING
sa_adm1_monthly<-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_monthly.csv")
all_names<- colnames(sa_adm1_monthly)
cov_names<-c(grep("tas",all_names,value=TRUE),grep("sp",all_names,value=TRUE),"hurs","absh","prlrsum","prlrmax","prlrmean","pm2p5","pm10" ,"o3","so2")
# dist_list<-lapply(cov_names, function(x) test_distribution(sa_adm1_monthly,x))
# grid.arrange(grobs=dist_list)

## scale prlr and absh
sa_adm1_monthly_sc<-sa_adm1_monthly
scale_covs<-c("hurs","absh","prlrsum","prlrmax","prlrmean")
sa_adm1_monthly_sc[, (scale_covs) := lapply(.SD, scale), .SDcols = scale_covs]
write.table(sa_adm1_monthly_sc,file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_monthly_sc.csv", quote=FALSE , row.names = FALSE, col.names = TRUE,sep=",")
# dist_list<-lapply(cov_names, function(x) test_distribution(sa_adm1_monthly_sc,x))
# grid.arrange(grobs=dist_list)



################################################################################
### LAG ENVIRONMENTAL VARIABLES ######,"o3","so2"
################################################################################       
sa_adm1_weekly<-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly.csv")
all_names<- colnames(sa_adm1_weekly)
cov_names_lag<-c(grep("tas",all_names,value=TRUE),grep("sp",all_names,value=TRUE),"hurs","absh","prlrsum","prlrmax","prlrmean","pm2p5","pm10","o3","so2","sfcWind",grep("monthly_count",all_names,value=TRUE))
sa_adm1_weekly_lag<-lag_dataframe(sa_adm1_weekly,cov_names_lag, 12)
colnames(sa_adm1_weekly_lag)[which(colnames(sa_adm1_weekly_lag)%in%cov_names_lag)] <- paste0(colnames(sa_adm1_weekly_lag)[which(colnames(sa_adm1_weekly_lag)%in%cov_names_lag)],"_lag0")
write.table(sa_adm1_weekly_lag,file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly_lag.csv", quote=FALSE , row.names = FALSE, col.names = TRUE,sep=",")

sa_adm1_weekly_sc<-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly_sc.csv")
all_names<- colnames(sa_adm1_weekly_sc)
cov_names_lag<-c(grep("tas",all_names,value=TRUE),grep("sp",all_names,value=TRUE),"hurs","absh","prlrsum","prlrmax","prlrmean","pm2p5","pm10","o3","so2","sfcWind",grep("monthly_count",all_names,value=TRUE))
sa_adm1_weekly_lag_sc<-lag_dataframe(sa_adm1_weekly_sc,cov_names_lag, 12)
colnames(sa_adm1_weekly_lag_sc)[which(colnames(sa_adm1_weekly_lag_sc)%in%cov_names_lag)] <- paste0(colnames(sa_adm1_weekly_lag_sc)[which(colnames(sa_adm1_weekly_lag_sc)%in%cov_names_lag)],"_lag0")
write.table(sa_adm1_weekly_lag_sc,file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly_lag_sc.csv", quote=FALSE , row.names = FALSE, col.names = TRUE,sep=",")

sa_adm1_monthly_sc<-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_monthly_sc.csv")
all_names<- colnames(sa_adm1_monthly_sc)
cov_names_lag<-c(grep("tas",all_names,value=TRUE),grep("sp",all_names,value=TRUE),"hurs","absh","prlrsum","prlrmax","prlrmean","pm2p5","pm10","o3","so2","sfcWind",grep("monthly_count",all_names,value=TRUE))
sa_adm1_monthly_lag_sc<-lag_dataframe(sa_adm1_monthly_sc,cov_names_lag, 6)
colnames(sa_adm1_monthly_lag_sc)[which(colnames(sa_adm1_monthly_lag_sc)%in%cov_names_lag)] <- paste0(colnames(sa_adm1_monthly_lag_sc)[which(colnames(sa_adm1_monthly_lag_sc)%in%cov_names_lag)],"_lag0")
write.table(sa_adm1_monthly_lag_sc,file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_monthly_lag_sc.csv", quote=FALSE , row.names = FALSE, col.names = TRUE,sep=",")

sa_adm1_monthly<-fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_monthly.csv")
all_names<- colnames(sa_adm1_monthly)
cov_names_lag<-c(grep("tas",all_names,value=TRUE),grep("sp",all_names,value=TRUE),"hurs","absh","prlrsum","prlrmax","prlrmean","pm2p5","pm10","o3","so2","sfcWind",grep("monthly_count",all_names,value=TRUE))
sa_adm1_monthly_lag<-lag_dataframe(sa_adm1_monthly,cov_names_lag, 6)
colnames(sa_adm1_monthly_lag)[which(colnames(sa_adm1_monthly_lag)%in%cov_names_lag)] <- paste0(colnames(sa_adm1_monthly_lag)[which(colnames(sa_adm1_monthly_lag)%in%cov_names_lag)],"_lag0")
write.table(sa_adm1_monthly_lag,file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_monthly_lag.csv", quote=FALSE , row.names = FALSE, col.names = TRUE,sep=",")


### check the order according to the shape files
idx1 <- unique(sa_adm1_weekly$GID_1)
idx2 <- unique(sa_adm1_weekly_lag_sc$GID_1)
idx3 <- unique(sa_adm1_monthly$GID_1)
idx4 <- unique(sa_adm1_monthly_lag_sc$GID_1)
idx5 <- unique(shp$GID_1)

(all(diff(match(idx5, idx2)) == 1))
(all(diff(match(idx5, idx3)) == 1))
(all(diff(match(idx5, idx4)) == 1))
(all(diff(match(idx5, idx1)) == 1))

idx1 <- unique(sa_adm1_weekly$id_u)
idx2 <- unique(sa_adm1_weekly_lag_sc$id_u)
idx3 <- unique(sa_adm1_monthly$id_u)
idx4 <- unique(sa_adm1_monthly_lag_sc$id_u)

(all(diff(match(idx1, idx2)) == 1))
(all(diff(match(idx1, idx3)) == 1))
(all(diff(match(idx1, idx4)) == 1))
(all(diff(match(idx1, idx2)) == 1))



# #### sanity check plots
# # disease counts by age
# sum<-sa_adm1_monthly_lag%>%
#   group_by(date) %>%
#   summarise(
#     age_lt1=sum(age_lt1,na.rm=T),
#     age_lt5=sum(age_lt5,na.rm=T),
#     age_5t14=sum(age_5t14,na.rm=T), 
#     age_14t24=sum(age_14t24,na.rm=T), 
#     age_gt65=sum(age_gt65,na.rm=T), 
#     age_gt80=sum(age_gt80,na.rm=T) 
#   )
# (ggplot(sum)+
#     geom_line(aes(x=date,y=age_lt1))+
#     ggtitle("Age <1")+
#     ylab("Count")+
#     ylim(0,200) + 
#     theme(axis.text = element_text(size=12), axis.title = element_text(size=12))+
#     theme_bw()) +
#   (ggplot(sum)+
#      geom_line(aes(x=date,y=age_lt5))+
#      ggtitle("Age <5")+
#      ylab("Count")+
#      ylim(0,200) + 
#      theme(axis.text = element_text(size=12), axis.title = element_text(size=12))+
#      theme_bw()) +
#   (ggplot(sum)+
#      geom_line(aes(x=date,y=age_5t14))+
#      ggtitle("Age 5-14")+
#      ylab("Count")+
#      ylim(0,200) + 
#      theme(axis.text = element_text(size=12), axis.title = element_text(size=12))+
#      theme_bw()) +
#   (ggplot(sum)+
#      geom_line(aes(x=date,y=age_14t24))+
#      ggtitle("Age 14-24")+
#      ylab("Count")+
#      ylim(0,200) + 
#      theme_bw()) +
#   (ggplot(sum)+
#      geom_line(aes(x=date,y=age_gt65))+
#      ggtitle("Age >65")+
#      ylab("Count")+
#      ylim(0,200) + 
#      theme(axis.text = element_text(size=12), axis.title = element_text(size=12))+
#      theme_bw()) +
#   (ggplot(sum)+
#      geom_line(aes(x=date,y=age_gt80))+
#      ggtitle("Age >80")+
#      ylab("Count")+
#      ylim(0,200) + 
#      theme(axis.text = element_text(size=12), axis.title = element_text(size=12))+
#      theme_bw()) 
# 
# ## disease counts by province
# sum<-sa_adm1_monthly_lag%>%
#   group_by(date, NAME_1) %>%
#   summarise(
#     counts = sum(disease, na.rm =T)
#   )
# ggplot(sum) + 
#   geom_line(aes(x=date,y=counts, group=NAME_1))+
#   ggtitle("")+
#   ylab("Count")+
#   theme_bw() +
#   theme(axis.text = element_text(size=12), axis.title = element_text(size=12))+
#   facet_wrap(NAME_1~.)
# 
# ggplot(sa_adm1_monthly_lag_sc)+
#   geom_line(aes(x=date,y=GPSC13_prop))
# 
# 
# ############################### 
# ### explore weekly ###
# 
# sa_adm1_weekly_lag_sc <- fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly_lag_sc.csv")
# table(table(sa_adm1_weekly_lag_sc$date))
# df_sum <- sa_adm1_weekly_lag_sc %>%
#   group_by(GID_1,adm1_name,date) %>%
#   summarise(count = sum(disease,na.rm=T))
# ggplot(df_sum) +
#   geom_line(aes(x=date,y=count,group=adm1_name,color=adm1_name,group=GID_1))+
#   theme_bw()+
#   theme(legend.position = "none")+
#   facet_grid(GID_1 ~. )
# 
# sa_adm1_weekly_lag <- fread(file="/home/sbelman/Documents/env_sa_manuscript/dataframes/sa_adm1_weekly_lag.csv")
# 
# df_sum <- sa_adm1_weekly_lag %>%
#   group_by(GID_1,adm1_name,date) %>%
#   summarise(count = mean(pm2p5_lag0,na.rm=T))
# ggplot(df_sum) +
#   geom_line(aes(x=date,y=count,group=adm1_name,color=adm1_name,group=GID_1))+
#   theme_bw()+
#   theme(legend.position = "none")+
#   facet_grid(GID_1 ~. )
