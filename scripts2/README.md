## Script Library ##
## Data Processing ##
*000_summarise_data.R* 

This summarises the line-list invasive pneumococcal disease data from GERMS-SA into a cohesive spatiotemporal data frame including geocoding hospitals and adding genomic data to epidemiological data. 

*0_process_data.R* & *0_process_data_adm1.R* 

These bring together the epidemiological, genomic, sociodemographic, environmental datasets together into either weekly or monthly resolution. They are for district level (adm2) and province level (adm1) respectively.

*02_write_formulas.R* 

Creates lists of formulas to use in different model configurations.
## Testing base models (no environmental covariates) ##
*03_base_model.R* 

This reads in formulas from a list in the previous script and runs base models. There are 2 flags which can be set in lines 15 and 16 respecitively specifying the time and space resolution whereby time can either be "weekly" or "monthly" and space is either "adm1" or "adm2". This dictates which dataframes are read in and where model outputs are saved.

*03_base_ARTPCV.R* 

This runs more base models testing different integrations of the antiretroviral therapy (ART) and pneumococcal conjugate vaccine (PCV) and their effects. Tested formulations including coverage (%) for each nationally, coverage (%) at the province level for ART, categorical covariates for vaccine period including 3 periods (prePCV, postPCV7, and postPCV13) and including only 2 periods (prePCV and postPCV). This script also includes supplementary plots incorporating the interannual random effect with thse interventions included as covariates as well as the fixed effects.
## Running Models with Environmental Covariates ##
These scripts include a base model with seasonal, spatial, interannual (province replicate) random effects, as well as vaccination period and population density covariates. They iteratively run across environmental factors of choice and save the goodness-of-fit metrics, summary of the model outputs, and the fits. The bivariable script includes an additional loop which allows the addition of a second environmental covariate.
Both scripts allow modifications to run models including invasive disease cases from 2005-2019 and 2005-2023, district or province level, and weekly or monthly models.

*04_run_univariable_dlnm_GPSC.R* This includes a single environmental factor.

*04.1_run_bivariable_dlnm_GPSC.R* This script allows inclusion of two environmental factors. Some exploratory plots are included.
## Running Models with Sensitivity Modifications ##
*06.1_run_univariable_dlnm_GPSC_singleprovinces.R*

This is the same as the above run_univariable_dlnm_GPSC script but allows a subset by just the Gauteng province or Western Cape province either at adm1 (province) or adm2 (district) level.

*06.3_run_univariable_dlnm_GPSC_PCVSensitivity.R*

This script stratifies the interaction model with each GPSC to run only on 2005-2008 (prePCV) and 2009-2019 or 2009-2023 to determine whether the exposure response curves are a result of the vaccine intervention or something that persists across both periods for the GPSCs or serotypes. This script includes summary plots of the pre-PCV and post-PCV periods.
## Running Models with outcome specific disease types ##
*06_run_univariable_dlnm_multioutcomes.R*

This script allows the outcome to be modified given demographic, or serotype variables specified in a vector. Includes space time modifications and only including 2005-2019 or including all of 2005-2023. 

