##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Model-based estimate of abundance index for EBS-NBS YFS
## Authors:       Zack Oyafuso (zack.oyafuso#noaa.gov)
##                Emily Markowitz (emily.markowitz#noaa.gov)
##                Jason Conner (jason.conner@noaa.gov)
##                James Thorson (james.thorson@noaa.gov)
## Output POC:    Ingrid Spies (ingrid.spies@noaa.gov)
## Description:   2022 hindcast index standardization for yellowfin sole
##                     (Limanda aspera) in the Eastern and Northern Bering Sea
##  NOTES: 
##                Make sure R and package versions are consistent with those 
##                versions speficied in the 2022 Terms of Reference (TOR)
##                https://docs.google.com/document/d/1t-pIruLZ-F_iNzCysWLH8cdsM0gZpuUb/edit
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
finalanalysis <- TRUE # this will make the model work faster while testing


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages ----
##   Notes: hints for downloading certain packages
##          coldpool: remotes::install_github("afsc-gap-products/coldpool")
##          VAST: see github page: https://github.com/James-Thorson-NOAA/VAST 
##                or remotes::install_github("James-Thorson-NOAA/VAST", 
##                                           ref = "3.9.0") 
##          INLA: source("http://www.math.ntnu.no/inla/givemeINLA.R")  
##          FishStatsUtils: remotes::install_github("James-Thorson-NOAA/FishStatsUtils@2.11.0")
##          Matrix: install.packages(pkgs = "https://cran.r-project.org/bin/windows/contrib/4.2/Matrix_1.4-0.zip", repos = NULL)
##          TMB: install.packages(pkgs = "https://cran.r-project.org/src/contrib/Archive/TMB/TMB_1.7.22.tar.gz", repos = NULL)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(coldpool) 
library(VAST) 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   System preferences ----
##   Updated every year
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
R_version <- "R version 4.0.2 (2020-06-22)"
current_year <- 2022
VAST_cpp_version <- "VAST_v13_1_0"
pck_version <- c("VAST" = "3.9.0", 
                 "FishStatsUtils" = "2.11.0", 
                 "Matrix" = "1.4-0", 
                 "TMB" = "1.7.22", 
                 "DHARMa" = "0.4.5")

{
  if(sessionInfo()$R.version$version.string == R_version) 
    message(paste0(sessionInfo()$R.version$version.string, 
                   " is consistent with the ", current_year, " TOR."))
  
  if(!sessionInfo()$R.version$version.string == R_version) 
    message(paste0("WARNING: ", sessionInfo()$R.version$version.string, 
                   " is NOT consistent with the ", current_year, " TOR. ",
                   "Please update R version to ", R_version))
  
  for (pck in 1:length(pck_version)) {
    temp_version <- packageVersion(pkg = names(pck_version)[pck])
    
    if(temp_version == pck_version[pck])
      message(paste0("The version of the '", names(pck_version)[pck], 
                     "' package (", temp_version, ") is consistent",
                     " with the ", current_year, " TOR."))
    
    if(!temp_version == pck_version[pck])
      message(paste0("WARNING: ", 
                     "The version of the '", names(pck_version)[pck], 
                     "' package (", temp_version, ") is NOT consistent",
                     " with the ", current_year, " TOR. Please update the '", 
                     names(pck_version)[pck], "' package to ", 
                     pck_version[pck]))
  }
  rm(pck, temp_version)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set species ----
##   Record session information
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species <- 10210
species_name <- "yellowfin_sole"

## Save package versions
sink(paste0("species_specific_code/BS/", species_name, 
            "/hindcast/results/session_info.txt"), 
     type = "output")
sessionInfo()
sink()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create results folder ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
folder <- paste0("species_specific_code/BS/", species_name, "/hindcast/")
if(!dir.exists(folder)) dir.create(folder)
folder <- paste0("species_specific_code/BS/", species_name,"/hindcast/results/")
if(!dir.exists(folder)) dir.create(folder)
rm(folder)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load VAST data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
Data_Geostat <- readRDS(file = paste0("species_specific_code/BS/", species_name, 
                                      "/hindcast/data/data_geostat_index.rds"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load coldpool covariate data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cpi <- scale(coldpool:::cold_pool_index$AREA_LTE2_KM2)
covariate_data <- data.frame(Year = c(coldpool:::cold_pool_index$YEAR, 2020),
                             Lat = mean(Data_Geostat$Lat),
                             Lon = mean(Data_Geostat$Lon), 
                             AREA_LTE2_KM2 = c(cpi, 0))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set VAST settings ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
settings <- FishStatsUtils::make_settings( 
  n_x = 750,
  Region = c("Eastern_Bering_Sea", "Northern_Bering_Sea"),
  purpose = "index2",
  fine_scale = TRUE,
  strata.limits = data.frame(STRATA = as.factor('All_areas')),
  ObsModel = c(2,1),
  FieldConfig = c("Omega1" = "IID", "Epsilon1" = "IID", 
                  "Omega2" = "IID", "Epsilon2" = "IID"),
  RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 4, "Epsilon2" = 4),
  OverdispersionConfig = c("Eta1" = 0, "Eta2" = 0),
  Options = c("Calculate_Range" = TRUE, 
              "Calculate_effective_area" = TRUE, 
              "treat_nonencounter_as_zero" = FALSE),
  use_anisotropy = TRUE,
  Version = VAST_cpp_version,
  max_cells = 2000,
  knot_method = "grid",
  bias.correct = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Run VAST model ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit <- FishStatsUtils::fit_model( 
  
  ## Input settings
  "working_dir" = paste0(getwd(),"/species_specific_code/BS/",
                         species_name, "/hindcast/results/"), 
  "settings" = settings, 
  "create_strata_per_region" = TRUE,
  
  ## Input Data
  "Lat_i" = Data_Geostat[,'Lat'], 
  "Lon_i" = Data_Geostat[,'Lon'], 
  "t_i" = Data_Geostat[,'Year'], 
  "c_i" = rep(0, nrow(Data_Geostat)), 
  "b_i" = Data_Geostat[,'Catch_KG'], 
  "a_i" = Data_Geostat[,'AreaSwept_km2'], 
  "v_i" = Data_Geostat[,'Vessel'],
  
  ## Output settings
  "getJointPrecision" = TRUE, 
  "getReportCovariance" = TRUE,
  "getsd" = ifelse(finalanalysis, TRUE, FALSE),
  
  ## Model tuning
  "newtonsteps" = ifelse(finalanalysis, 1, 0),
  "Npool" = 100,
  
  ## Covariate data
  "X1_formula"= ~ AREA_LTE2_KM2,
  "X2_formula"= ~ AREA_LTE2_KM2,
  "X1config_cp" <- as.matrix(2),
  "X2config_cp" <- as.matrix(2) ,
  "covariate_data" = covariate_data)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save results locally ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Save VAST fit object
saveRDS(object = fit, 
        file = paste0("species_specific_code/BS/",
                      species_name, "/hindcast/results/VASTfit.RDS"))

## Save diagnostics and other outputs
results <- FishStatsUtils::plot_results( 
  fit = fit, 
  working_dir = paste0(getwd(),  "/species_specific_code/BS/",
                       species_name, "/hindcast/results/output_plots/"),
  zrange = c(-3, 3), 
  n_cells = 600, 
  strata_names = c("Both", "EBS", "NBS"), 
  check_residuals = TRUE,
  n_samples = 0)

saveRDS(object = results, 
        file = paste0("species_specific_code/BS/", species_name, 
                      "/hindcast/results/output_plots/diagnostics.RDS"))

## ESP products
write.csv(x = results$Range$COG_Table, 
          file = paste0("species_specific_code/BS/", species_name, 
                        "/hindcast/results/output_plots/COG.csv"), 
          row.names = FALSE)

## Effective area occupied for ESP request
report <- TMB::summary.sdreport(fit$parameter_estimates$SD)
ln_km2 <- report[which(rownames(report)=="log_effective_area_ctl"),
                 c('Estimate','Std. Error')]
Year <- sort(unique(fit$year_labels))
Area <- rep(x = c("Both", "EBS", "NBS"), each = length(Year))
ln_km2 <- as.data.frame(cbind(Area, Year, ln_km2))
ln_km2 <- ln_km2[which(ln_km2$Year %in% unique(fit$data_frame$t_i)),]
write.csv(x = ln_km2, 
          file = paste0("species_specific_code/BS/", species_name, 
                        "/hindcast/results/output_plots/ln_effective_area.csv"),
          row.names = FALSE )
