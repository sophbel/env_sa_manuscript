# Pneumococcal population structure influences the effects of air pollution on invasive disease risk in South Africa
 *Streptococcus pneumoniae* is highly diverse, comprising over 100 serotypes and hundreds of genomic lineages amidst widespread vaccination. While it can cause invasive pneumococcal disease (IPD) which exhibits pronounced seasonal spikes, the interplay between pneumococcal diversity and environmental drivers remains unexplored. Here, we analyzed 59,017 IPD cases over 19 years from South Africa, incorporating 4,350 genome-sequenced isolates, using Bayesian spatio-temporal models to link environmental exposures and pneumococcal diversity. Cumulatively, across an 8-week period, moderate relative humidity (33-49%) and cold minimum temperatures (4-10°C) increased IPD risk by 5% and 4%, respectively. Conversely, warm maximum temperatures (27-38°C) were associated with up to a 10% increased risk within a week of exposure. There was a positive association between air pollution (PM<sub>2.5</sub>) and IPD although it varied by age, disease presentation, and most notably serotype and lineage. Specifically, the lag time between PM<sub>2.5</sub> exposure and disease onset varied by serotype with only serotypes 4, 8, and 23F conferring an immediate IPD risk. High prevalence of GPSC21 lineage (serotype 19F) also modified the pollution response, shifting the lag structure to produce immediate risk of disease following high PM<sub>2.5</sub> exposure. Our results demonstrate that pneumococcal population structure shapes air quality risk which in turn can shape the fitness landscape of microbial populations. Integration of these data may inform public health policy.

 ## Folders ##
```./scripts2/``` folder contains the scripts required for this code to run - including the data processing and models.

```./scripts2/environmental_processing``` contains the scripts to process the climatic and air pollution datasets including comparisons of the observational and reanalysis air pollution datasets.

```./dataframes/``` add the harmonized datasets from Zenodo to this folders (e.g. sa_adm2_weekly_lag.csv)

```./input_datasets/shps/``` add the shapefiles from Zenodo to this folders (e.g. gadm41_ZAF_xx.shx)

 ## Code Availability ##
 Purpose of manuscript is investigating the impact of environmental factors (meteorological and air pollution), while accounting for underlying microbial diversity in invasive pneumococcal disease. 

All code for analysis is included in this repository and the associated packages can be installed with the yaml file.

__R Libraries__

```
if (!require("yaml")) install.packages("yaml")
if (!require("remotes")) install.packages("remotes")

# Read the YAML file
env <- yaml::read_yaml("r_environment.yml")

# Install the packages via 'remotes' to ensure version matching
invisible(lapply(names(env$packages), function(pkg) {
  version <- env$packages[[pkg]]
  message(paste("Installing", pkg, "version", version))
  remotes::install_version(pkg, version = version, upgrade = "never")
}))
```

## Data Availability ##
__Data and Models__

Available on Zenodo (https://doi.org/10.5281/zenodo.17877254)

Saved as three separate tar.gz files which can be unzipped locally into the appropriate folders at:

**input files:** input_datasets.tar.gz

**main harmonized dataframes:** dataframes.tar.gz

**model outputs:** models.tar.gz

__Genomic Data__

The microbial genomic data used in this study were previously published in Lekhuleni et al., 2024 (doi: https://doi.org/10.1038/s41467-024-52459-3) and the accession numbers are included in Supplementary Table S17. This project is part of the Global Pneumococcal Sequencing Project (https://www.pneumogen.net/gps/)


