########### LOAD LIBRARIES #############3
library(data.table)
library(tidyverse)
library(lubridate)
library(tidyr)
library(dplyr)
library(tsModel)
library(sf)
library(tmaptools)
library(tidygeocoder)
library(osmdata)
library(gridExtra)
library(readxl)
library(readr)
# setwd("/home/sbelman/Documents/env_sa_manuscript/")
################ RENAME AND SELECT COLUMNS REQUIRED FROM GERMS-SA ##########
'%notin%'<-Negate('%in%')
# setwd("")
### read in raw data file
setwd("/home/sbelman/Documents/BRD/SouthAfrica/disease/germs/")
data<-fread("/home/sbelman/Documents/BRD/SouthAfrica/disease/germs/GERMS_SP_2003-2023_5Jul2024.csv")
data<-data.frame(data)
data<-data[,c("Lane_id","YEAR","COLLECTDTE",colnames(data)[grep("HOSPITAL",colnames(data))],"SSEROTYPE","AGEMONTHS","AGEYEARS","SPECIMENTYPE","SPECDIAG","SPEN", "PREGNA", "SMOKER", "SEX")]
colnames(data)<-c("laneid","year","date","hospital_name","province","district","subdistrict","serotype","agemonths","ageyears","spec_type","spec_diagnosis","spec_penR","pregnant","smoker","sex")
data$district[data$district==""]<-"unknown"
data$province[data$province==""]<-"unknown"
data$hospital_name[data$hospital_name==""]<-"unknown"

#################### accessions ##################
#### exploring which datasets had accession numbers including GPS project from ENA, Monocle itself, and published paper Lekhuleni, et al., 2023. 
#### accessions are included in the supplementary material.

# data<-fread("input_datasets/disease/SA_disease_point.csv")
# # only 2005-2020 of data
# dtmpac <- data[which(!is.na(data$laneid)& data$year<2020),]
# 
# # 1. Fetch the entire sequencing run report for the GPS project from ENA
# url <- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB3084&result=read_run&fields=study_accession,sample_accession,run_accession,country,fastq_ftp"
# gps_metadata <- read_tsv(url)
# germs_sa_pneumo <- gps_metadata %>% 
#   filter(grepl("South Africa", country, ignore.case = TRUE))
# 
# ### read in genomic accessions to check they match
# accs <- read_excel("./input_datasets/germs/41467_2024_52459_MOESM10_ESM.xlsx")
# accs <- accs[,c("Isolate ID","Isolate accession numbers")]
# colnames(accs) <- c("laneid","accession","year")
# 
# ### read in from monocle
# monc <- fread("./input_datasets/germs/monocle-metadata_all.csv")
# monc <- monc[,c("Lane_id","ERS","ERR")]
# colnames(monc) <- c("laneid","ERS","ERR")
# 
# # which are present
# missingaccs <- dtmpac[which(dtmpac$laneid%notin%monc$laneid),]
# newmiss <- missingaccs[which(missingaccs$laneid%notin%accs$`Isolate ID`),]
# data <- data |> subset(district!="")

############### ASSIGN DISEASE COLUMN ######################
data$disease<-1
# Extract epi week and epi year
data$date<-as.Date(data$date)
data$epi_week <- epiweek(data$date)
data$epi_year <- epiyear(data$date)

# Separate district names
data<- data %>% separate(district, into= "district", sep= " District Municipality") %>%
  separate(district, into="district", sep = " Metropolitan Municipality")

# Add a 1 where there is a sequence and 0 where there is not
data$sequenced <- ifelse(data$laneid=="",0,1)
data$laneid[which(data$laneid=="")]<-NA
# Save processed file including disease
# write.table(data,"input_datasets/disease/SA_disease_point.csv",quote=FALSE,row.names = FALSE, col.names = TRUE,sep=",")

# data5 <- data[,-c("laneid")]
# write.table(data5,"input_datasets/disease/SA_disease_point_share.csv",quote=TRUE,row.names = FALSE, col.names = TRUE,sep=",")

################################################################################
## SAVE NAME MATCHED SHAPE FILES 
################################################################################
# LOAD SHAPE FILES
data<-fread("input_datasets/disease/SA_disease_point_share.csv")

######## district level
shp<-st_read("input_datasets/shps/gadm41_ZAF_2.shp")
shp_name_vec<-shp$NAME_2
dist_vec<-unique(data$district)
table(dist_vec%in%shp_name_vec)
dist_vec[dist_vec %notin% shp_name_vec]
shp_name_vec[shp_name_vec %notin% dist_vec]

## reassign the names to match between the data and the shape file
shp <- shp %>%
  mutate(district = case_when(
    NAME_2 == "O.R.Tambo" ~ "OR Tambo",
    NAME_2 == "Bojanala" ~ "Bojanala Platinum",
    NAME_2 == "Thabo Mofutsanyane" ~ "Thabo Mofutsanyana",
    NAME_2 == "Uthungulu" ~ "uThungulu",
    NAME_2 == "Uthukela" ~ "uThukela",
    NAME_2 == "Umzinyathi" ~ "uMzinyathi",
    NAME_2 == "Umgungundlovu" ~ "uMgungundlovu",
    NAME_2 == "Umkhanyakude" ~ "uMkhanyakude",
    NAME_2 == "Sisonke" ~ "Harry Gwala",
    NAME_2 == "Cacadu" ~ "Sarah Baartman",
    TRUE ~ NAME_2  # retain original value if no match
  )) 
st_write(shp, "input_datasets/shps/gadm41_namematch_ZAF_2.shp", delete_dsn=FALSE)


######## repeat for province level 
shp1<-st_read("input_datasets/shps/gadm41_ZAF_1.shp")
shp_name_vec<-shp1$NAME_1
prov_vec<-unique(data$province)
table(prov_vec%in%shp_name_vec)
prov_vec[prov_vec %notin% shp_name_vec]
shp_name_vec[shp_name_vec %notin% prov_vec]

st_write(shp1, "input_datasets/shps/gadm41_namematch_ZAF_1.shp", delete_dsn=FALSE)



################################################################################
### USE COMBINATION OF OPEN STREET MAPS AND OTHER GEOCODE TO ASSIGN 
### LATS AND LONGS TO HOSPITALS WHERE POSSIBLE AND DISTRICTS WHERE NOT
################################################################################
tmp <- data
      # test geocoding the hospitals
      # tmp$hospitalA<-paste0(tmp$hospital_name,", ",tmp$district, ", District, South Africa")
      # tmpa<-sort(unique(tmp$hospitalA))
### save data necessary
hospitals<-unique(tmp[,c("hospital_name","district","province")])

# Geocoding using tidygeocoder
geocoded_hospitals <- hospitals %>%
  mutate(address = paste(hospital_name,",", district,",", province, ",South Africa")) %>%
  geocode(address, method = 'osm', full_results = TRUE)

## save province level geocode
saveRDS(geocoded_hospitals, file="input_datasets/germs/geocoded_hospitals.rds" )


##### prepare district level geocode
## subset those which are missing
nolat<-geocoded_hospitals[is.na(geocoded_hospitals$lat),]

# for those with no data geocode to district level
inp<-paste(nolat$district,",", nolat$province, ", South Africa")
inp[inp=="OR Tambo , Eastern Cape , South Africa"]<- "OR Tambo District Municipality , Eastern Cape , South Africa"
geocoded_districts <-  geocode_OSM(inp,  as.data.frame = TRUE)
# subset query and coordinates assign hospital name to NA
geocoded_districts<-geocoded_districts[c("query","lat","lon")]

### test the number of districts
test<-geocoded_districts
test <- test %>%
  separate(query, into = c("district", "province", "country"), sep = " , ") %>%
  dplyr::select(-country)

hosp_ind<-unique(as.data.frame(nolat[c("district","province","hospital_name")]))
test<-unique(test)
test1<-hosp_ind%>%
  left_join(test,by=c("district","province"))
geocoded_districts_new<-test1

# subset query and coordinates from hospital data frame
geocoded_hospitals <- geocoded_hospitals[c("hospital_name","district","province","lat","long")]
## add geocoded districts to this data frame
nas<-geocoded_hospitals[is.na(geocoded_hospitals$lat),]
hosp_vec<-unique(nas$hospital_name)
for(i in 1:length(hosp_vec)){
  nas[nas$hospital_name==hosp_vec[i],c("lat","long")] <- geocoded_districts_new[geocoded_districts_new$hospital_name==hosp_vec[i],c("lat","lon")]
}
nas[is.na(nas$lat),"lat"]<--26.134789
nas[is.na(nas$long),"long"]<-28.240528

## merge all lats and longs
geocoded_hospitals_dists<-rbind(nas,geocoded_hospitals[!is.na(geocoded_hospitals$lat),])
saveRDS(geocoded_hospitals_dists, file="input_datasets/germs/geocoded_hospitals_districts.rds" )



################################################################################
### SUMMARISE DISTRICT COUNCIL PROJECTIONS FROM MID-YEAR POPULATION ESTIMATES FROM STATISTICS SOUTH AFRICA
################################################################################

population_statsZA <- read.xlsx2("input_datasets/sociodemographic/District_Council_projection_by_sex_and_age_2002-2024.xlsx",  sheetIndex = 1 , header = TRUE)
population_long <- population_statsZA %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "year",
    values_to = "population_size"
  ) %>%
  mutate(
    year = str_remove(year, "^X"),
    population_size = as.numeric(population_size)
  )
colnames(population_long)<- c("district","sex","age","year","population_size")
population_long <- population_long %>%
  mutate(
    district_clean = district %>%
      # 1. Strip prefixes and suffixes
      str_remove("^[A-Z]{2,3}\\s*-\\s*") %>%         # remove province code prefix
      str_remove_all("\\(DC[0-9]+\\)") %>%           # remove (DCxx)
      str_remove_all("\\(MAN\\)") %>%                # remove (MAN) from Mangaung
      str_remove_all("District Municipality") %>%    # drop "District Municipality"
      str_remove_all("Metropolitan Municipality") %>%# drop "Metropolitan Municipality"
      str_remove_all("DM") %>%                       # drop "DM"
      str_squish() %>%                               # trim extra spaces
      # 2. Handle known mismatches with landscan
      str_replace("^Garden Route$", "Eden") %>%
      str_replace("^ZF Mgcawu$", "Siyanda") %>%
      str_replace("^King Cetshwayo$", "uThungulu")
  )

population_long <- subset(population_long, population_long$district_clean !="")
population_long$district_clean[which(population_long$district_clean=="MP - Ehlanzeni")] <- "Ehlanzeni"
population_long$district_clean[which(population_long$district_clean=="Thabo Mofutsanyane")] <- "Thabo Mofutsanyana"

population_data <- population_long %>%
  group_by(district_clean,year) %>%
  summarise(population_size = sum(population_size, na.rm = T))

## merge with shp for the GID values
colnames(population_data)[which(colnames(population_data)=="district_clean")] <- "district"
population_data <- left_join(population_data,st_drop_geometry(shp[,c("district","GID_1","GID_2")]), by = "district")
population_data$population_size <- as.numeric(population_data$population_size)
population_data$year <- as.numeric(population_data$year)
population_data <- data.table(population_data)

write.table(population_data, file = "input_datasets/sociodemographic/population_2000_2023_statsZA.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)



################################################################################
### PLOTS OF WHERE SEQUENCED ISOLATES ARE LOCATED 
################################################################################
############### sequence map plots ##############
# ## plot those we have seqs for by district
# data_seq<-data[which(data$sequenced==1),]
# weekly_data <- data_seq %>%
#   group_by(date,epi_year, epi_week, district) %>%
#   summarise(weekly_value = sum(disease), .groups = 'drop')
# 
# ggplot(weekly_data)+
#   geom_line(aes(x=date,y=weekly_value))+
#   theme_classic()+
#   xlab("Year")+
#   ylab("Weekly Case Count")
# 
# #### sum across districts over all time
# dist_count <- data_seq %>%
#   group_by(district) %>%
#   summarise(dist_count = sum(disease), .groups = 'drop')
# shp_dist_seq <- left_join(shp, dist_count,by="district")
# table(shp_dist_seq$dist_count!=0)
# seq_plot <- ggplot(shp_dist_seq) +
#   geom_sf(aes(fill=dist_count), color=NA) +
#   scale_fill_gradient(low="#789EBF",high="#AA336A")+
#   theme_bw() +
#   theme(legend.position = "bottom", legend.text = element_text(size = 10, angle = 45, hjust = 1))+
#   labs(fill="") +
#   ggtitle("Genomic Sequence Count")
# 
# dis_count <- data %>%
#   group_by(district)%>% 
#   summarise(disease_count = sum(disease), .groups = 'drop')
# shp_dist_disease <- left_join(shp, dis_count,by="district")
# case_plot <- ggplot(shp_dist_disease) +
#   geom_sf(aes(fill=disease_count), color=NA) +
#   scale_fill_gradient(low="#789EBF",high="#AA336A")+
#   theme_bw() +
#   theme(legend.position = "bottom", legend.text = element_text(size = 10, angle = 45, hjust = 1))+
#   labs(fill="") +
#   ggtitle("Disease Count")
# 
# case_plot + seq_plot
#   
# ### plot the cases weekly per district
# ggplot(weekly_data) +
#   # ggplot(subset(weekly_data,weekly_data$district=="City of Johannesburg" & weekly_data$epi_year==2005)) +
#   geom_line(aes(x = date, y = weekly_value, group = district, color = district)) +
#   theme_classic() +
#   theme(legend.position = "none",axis.text = element_text(size=12),axis.title = element_text(size=12), axis.text.x = element_text(angle=90)) +
#   xlab("Date")+
#   ylab("Weekly Case Count")+
#   facet_wrap(~district, nrow=8)
# # ggsave("./data_summary/district_weekly.png",width=10,height=8)
# 
# 
# ##### plot the seasonality colored by year
# ##### plot the number from each province sequenced over time
# data_seq <- subset(data_seq, data_seq$province!="unknown")
# df<-as.data.frame.matrix( table(data_seq$year,data_seq$province))
# df$year<-rownames(df)
# # Reshape the data from wide to long format
# df_long <- df %>%
#   pivot_longer(cols = -year, # All columns except 'year'
#                names_to = "province", # New column for provinces
#                values_to = "count") # New column for the counts
# df_long$year<-as.numeric(df_long$year)
# ggplot(df_long)+geom_line(aes(x=year,y=count,color=province))+theme_classic() +
#   theme(axis.text = element_text(size=15),axis.title = element_text(size=15)) +
#   labs(color="Province") +
#   xlab("Year")

################################################################################
#### PLOT HOSPITAL LOCATIONS ON MAP
################################################################################
# geocoded_hospitals_dists<-readRDS(file="/home/sbelman/Documents/BRD/SouthAfrica/disease/germs/geocoded_hospitals_districts.rds" )
# 
# # TURN TO AN SF
# points_sf <- st_as_sf(geocoded_hospitals_dists, coords = c("long", "lat"), crs = 4326)
# # Check the CRS of the shapefile and points
# correct<-st_crs(shp)
# st_crs(points_sf)
# # If needed, transform the CRS of points to match the shapefile
# points_sf <- st_transform(points_sf, correct)
# # st_crs(points_sf)
# 
# # PLOT SHP
# pdist<-ggplot(shp)+
#   geom_sf()+
#   # geom_sf(data = points_sf,  size = 3) +
#   geom_sf(data = points_sf, aes(color = factor(district)), size = 3) +
#   theme_classic()+
#   theme(legend.position = "none")+
#   labs(color="Districts")
# pprov<-ggplot(shp)+
#   geom_sf()+
#   # geom_sf(data = points_sf,  size = 3) +
#   geom_sf(data = points_sf, aes(color = factor(province)), size = 3) +
#   theme_classic()+
#   # theme(legend.position = "none")+
#   labs(color="Province")
# 
#  pprov
# total_points<-dim(points_sf)[1]
# number_hosps<-length(unique(data$hospital_name))
# 
# ggsave(pdist, file="./data_summary/district_map.png",width=10,height=6)
# ggsave(pprov, file="./data_summary/province_map.png",width=10,height=10)
# 
# # ggsave("./data_summary/district_weekly.png",width=10,height=8)