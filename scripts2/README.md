## Script Library ##
## Data Processing ##
1) 000_summarise_data.R: This summarises the line-list invasive pneumococcal disease data from GERMS-SA into a cohesive spatiotemporal data frame including geocoding hospitals and adding genomic data to epidemiological data. 
2) 0_process_data.R and 0_process_data_adm1.R: These bring together the epidemiological, genomic, sociodemographic, environmental datasets together into either weekly or monthly resolution. They are for district level (adm2) and province level (adm1) respectively.
3) 02_write_formulas.R: Creates lists of formulas to use in different model configurations.
## Testing base models (no environmental covariates) ##
5) 03_base_model.R: This reads in formulas from a list in the previous script and runs base models. There are 2 flags which can be set in lines 15 and 16 respecitively specifying the time and space resolution whereby time can either be "weekly" or "monthly" and space is either "adm1" or "adm2". This dictates which dataframes are read in and where model outputs are saved.
6) 03_base_ARTPCV.R: This runs more base models testing different integrations of the antiretroviral therapy (ART) and pneumococcal conjugate vaccine (PCV) and their effects. Tested formulations including coverage (%) for each nationally, coverage (%) at the province level for ART, categorical covariates for vaccine period including 3 periods (prePCV, postPCV7, and postPCV13) and including only 2 periods (prePCV and postPCV). This script also includes supplementary plots incorporating the interannual random effect with thse interventions included as covariates as well as the fixed effects.
## Running Models with Environmental Covariates ##
7) 04_run_univariable_dlnm_GPSC.R & 04.1_run_bivariable_dlnm_GPSC.R: These scripts include a base model with seasonal, spatial, interannual (province replicate) random effects, as well as vaccination period and population density covariates. They iteratively run across environmental factors of choice and save the goodness-of-fit metrics, summary of the model outputs, and the fits. The bivariable script includes an additional loop which allows the addition of a second environmental covariate.
Both scripts allow modifications to run models including invasive disease cases from 2005-2019 and 2005-2023, district or province level, and weekly or monthly models.
## Running Models with Sensitivity Modifications ##

