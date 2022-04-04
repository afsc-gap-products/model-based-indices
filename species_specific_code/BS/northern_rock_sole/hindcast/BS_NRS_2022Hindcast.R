##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Hindcast Index Calculation, EBS northern rock sole
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov), modifying author
##                Jason Conner (jason.conner@noaa.gov), original author
## Output POC:    Carey McGilliard (carey.mcgilliard AT noaa.gov)
## Description:   2022 hindcast index standardization for northern rock sole
##                     (Lepidopsetta polyxystra) in the Eastern Bering Sea
## NOTES: 
##                Make sure R and package versions are consistent with the 
##                versions speficied in the 2022 Terms of Reference (TOR)
## https://docs.google.com/document/d/1t-pIruLZ-F_iNzCysWLH8cdsM0gZpuUb/edit
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
finalanalysis <- T # this will make the model work faster while troubleshooting

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# Import packages ---------------------------------------------------------

library(coldpool)
library(tidyverse)
library(VAST)
library(sf)
library(scales)
library(DHARMa)
library(Matrix)

# source("http://www.math.ntnu.no/inla/givemeINLA.R")  
# remotes::install_github("James-Thorson-NOAA/VAST", ref="3.8.2") 
# remotes::install_github("nwfsc-assess/geostatistical_delta-GLMM", ref="3.3.0") 
# remotes::install_github("James-Thorson-NOAA/FishStatsUtils@2.10.2")

# install.packages(pkgs = "https://cran.r-project.org/bin/windows/contrib/4.2/Matrix_1.4-0.zip", repos = NULL)
# install.packages(pkgs = "https://cran.r-project.org/src/contrib/Archive/TMB/TMB_1.7.22.tar.gz", repos = NULL)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   System preferences ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
R_version <- "R version 4.0.2 (2020-06-22)"
VAST_cpp_version <- "VAST_v13_1_0"
pck_version <- c("VAST" = "3.9.0", 
                 "FishStatsUtils" = "2.11.0", 
                 "Matrix" = "1.4-0", 
                 "TMB" = "1.7.22", 
                 "DHARMa" = "0.4.5")

{
  if(sessionInfo()$R.version$version.string == R_version) 
    message(paste0(sessionInfo()$R.version$version.string, 
                   " is consistent with the 2022 TOR."))
  
  if(!sessionInfo()$R.version$version.string == R_version) 
    message(paste0(sessionInfo()$R.version$version.string, 
                   " is NOT consistent with the 2022 TOR. ",
                   "Please update R version to ", R_version))
  
  for (pck in 1:length(pck_version)) {
    temp_version <- packageVersion(pkg = names(pck_version)[pck])
    
    if(temp_version == pck_version[pck])
      message(paste0("The version of the '", names(pck_version)[pck], 
                     "' package (", temp_version, ") is consistent",
                     " with the 2022 TOR."))
    
    if(!temp_version == pck_version[pck])
      message(paste0("The version of the '", names(pck_version)[pck], 
                     "' package (", temp_version, ") is NOT consistent",
                     " with the 2022 TOR. Please update the '", 
                     names(pck_version)[pck], "' package to ", pck_version[pck]))
  }
  rm(pck, temp_version)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Result dir ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
species_name <- "northern_rock_sole"
folder <- paste0("species_specific_code/BS/", species_name, "/hindcast/")
if(!dir.exists(folder)) dir.create(folder)
folder <- paste0("species_specific_code/BS/", species_name,"/hindcast/results/")
if(!dir.exists(folder)) dir.create(folder)
rm(folder)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load VAST data ----
##   For northern rock sole, only include records from 1996 on &
##   exclude the northern bering sea
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_geostat <- readRDS(file = paste0("species_specific_code/BS/", 
                                      species_name, 
                                      "/hindcast/data/EBS_NBS_Index.RDS"))

data_geostat <- subset(x = data_geostat, 
                       subset = YEAR >= 1996 & REGION == "EBS")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Covariate data ----
##   Cold pool index from coldpool package
##   Area (square km) <= 2 degrees Celcius, scaled
##   Clunky addition: 2020 covariate data needs to be specified even though
##   no data were collected that year. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
cpi <- scale(coldpool:::cold_pool_index$AREA_LTE2_KM2)
covariate_data <- data.frame(Year = c(coldpool:::cold_pool_index$YEAR, 2020),
                             Lat = mean(data_geostat$START_LATITUDE),
                             Lon = mean(data_geostat$START_LONGITUDE), 
                             AREA_LTE2_KM2 = c(cpi, 0))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set VAST settings ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
settings <- FishStatsUtils::make_settings( 
  n_x = 750,
  Region = "Eastern_Bering_Sea",
  purpose = "index2",
  fine_scale = TRUE,
  strata.limits = data.frame(STRATA = as.factor('All_areas')),
  ObsModel = c(2,1),
  FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", 
                  "Omega2"="IID", "Epsilon2"="IID"),
  RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0),
  OverdispersionConfig = c("Eta1"=0, "Eta2"=0),
  Options = c("Calculate_Range" = TRUE, 
              "Calculate_effective_area" = TRUE, 
              "treat_nonencounter_as_zero" = FALSE ),
  use_anisotropy = TRUE,
  Version = VAST_cpp_version,
  max_cells = 2000,
  knot_method = "grid",
  bias.correct = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Run VAST model ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
fit <- FishStatsUtils::fit_model( 
  "settings" = settings, 
  "Lat_i" = data_geostat[, "START_LATITUDE"], 
  "Lon_i" = data_geostat[, "START_LONGITUDE"], 
  "t_i" = data_geostat[, "YEAR"], 
  "c_i" = rep(0, nrow(data_geostat)), 
  "b_i" = data_geostat[, "WEIGHT"], 
  "a_i" = data_geostat[, "EFFORT"], 
  "v_i" = rep("missing", nrow(data_geostat)),
  "create_strata_per_region" = TRUE,
  "getJointPrecision" = TRUE, 
  "getReportCovariance" = TRUE,
  "X1_formula"= ~ AREA_LTE2_KM2,
  "X2_formula"= ~ AREA_LTE2_KM2,
  "X1config_cp" <- as.matrix(2),
  "X2config_cp" <- as.matrix(2) ,
  "covariate_data" = covariate_data,
  "Npool" = 100,
  "test_fit" = TRUE,
  "working_dir" = paste0(getwd(),"/species_specific_code/BS/",
                         species_name, "/hindcast/results/"), 
  getsd = ifelse(test = finalanalysis, yes = TRUE, no = FALSE),
  newtonsteps = ifelse(test = finalanalysis, yes = 1, no = 0))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save results locally ----
## VAST fit: fit --> VASTfit.RDS
## session info: sessionInfo() --> session_info.txt
## diagnostic plots: written to diagnostics_plots/
## center of gravity: results$Range$COG_Table --> COG.csv
## effective area ln_km2 --> occupied: ln_effective_area.csv
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

## Save VAST fit object
saveRDS(object = fit, 
        file = paste0("species_specific_code/BS/",
                      species_name, "/hindcast/results/VASTfit.RDS"))

## Save package versions
sink(paste0("species_specific_code/BS/",
            species_name, "/hindcast/results/session_info.txt"), 
     type = "output")
sessionInfo()
sink()

# Save diagnostics and other outputs
results <- FishStatsUtils::plot_results( 
  fit = fit, 
  working_dir = paste0(getwd(),  "/species_specific_code/BS/",
                       species_name, "/hindcast/results/output_plots/"),
  zrange = c(-3,3), 
  n_cells = 600, 
  strata_names = c("EBS"), 
  check_residuals=TRUE,
  n_samples = 0 )

saveRDS(object = results, 
        file = paste0("species_specific_code/BS/", species_name, 
                      "/hindcast/results/output_plots/diagnostics.RDS"))

## ESP products
write.csv( x = results$Range$COG_Table, 
           file = paste0("species_specific_code/BS/", species_name, 
                         "/hindcast/results/output_plots/COG.csv"), 
           row.names=FALSE )

## Effective area occupied for ESP request
report = TMB::summary.sdreport(fit$parameter_estimates$SD)
ln_km2 = report[which(rownames(report)=="log_effective_area_ctl"),
                c('Estimate','Std. Error')]
Year <- sort(unique(fit$year_labels))
Area <- rep(x = "EBS", each = length(Year))
ln_km2 <- as.data.frame(cbind(Area, Year, ln_km2))
ln_km2 <- ln_km2[which(ln_km2$Year %in% unique(fit$data_frame$t_i)),]
write.csv(x = ln_km2, 
          file=paste0("species_specific_code/BS/", species_name, 
                      "/hindcast/results/output_plots/ln_effective_area.csv"),
          row.names = FALSE )
