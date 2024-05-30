##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Model-based estimate of abundance index for EBS Kamchatka
## Survey POC:    Lewis Barnett (lewis.barnett@noaa.gov)
## Output POC:    Meaghan Bryan (meaghan.bryan@noaa.gov)
## Description:   2024 hindcast index standardization for Kamchatka flounder
##                     in the Eastern Bering Sea
##  NOTES: 
##                Make sure R and package versions are consistent with those 
##                versions specified in the 2024 Terms of Reference (TOR)
##                https://drive.google.com/file/d/1TmW9gySrWVGxv5OXpTbkMml6vAbi14Xl/view?usp=drive_link
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
which_model <- c("hindcast", "production")[1] # specify by changing index
species_name <- "kamchatka_flounder"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
remotes::install_github("afsc-gap-products/coldpool")
library(coldpool) 
library(VAST) 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   System preferences ----
##   Updated every year
#    if you have newer versions than those below that is completely fine!
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
current_year <- 2024
VAST_cpp_version <- "VAST_v14_0_1" # or VAST_v13_1_0 or later
# if you have newer versions than those below that is completely fine!
pck_version <- c("VAST" = "3.10.0",
                 "FishStatsUtils" = "2.12.0",
                 "Matrix" = "1.5-3",
                 "TMB" = "1.9.2",
                 "DHARMa" = "0.4.6")

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
  
  rm(pck, temp_version)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create results folder ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
folder <- paste0(getwd(), "/species_specific_code/BS/", 
                 species_name, "/", which_model, "/")
if(!dir.exists(paste0(folder, "results/"))) 
  dir.create(path = paste0(folder, "results/"), recursive = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Record sessionInfo ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink(file = paste0(folder, "results/session_info.txt"), 
     type = "output")
sessionInfo()
sink()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load VAST data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
Data_Geostat <- readRDS(file = paste0(folder, "data/data_geostat_index.RDS"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load coldpool covariate data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CPE <- scale(coldpool:::cold_pool_index$AREA_LTE2_KM2)
covariate_data <- data.frame(Year = c(coldpool:::cold_pool_index$YEAR, 2020),
                             Lat = mean(Data_Geostat$Lat),
                             Lon = mean(Data_Geostat$Lon), 
                             CPE = c(CPE, 0))
covariate_data <- covariate_data[covariate_data$Year > 1990,]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set VAST settings ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
settings <- FishStatsUtils::make_settings( 
  n_x = 750,
  Region = "Eastern_Bering_Sea",
  purpose = "index2",
  fine_scale = TRUE,
  strata.limits = data.frame(STRATA = as.factor('All_areas')),
  ObsModel = c(2, 1),
  FieldConfig = c("Omega1" = "IID", "Epsilon1" = "IID", 
                  "Omega2" = "IID", "Epsilon2" = "IID"),
  RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 4, "Epsilon2" = 4),
  OverdispersionConfig = c("Eta1" = 0, "Eta2" = 0),
  # Options = c("Calculate_Range" = TRUE, 
  #             "Calculate_effective_area" = TRUE, 
  #             "treat_nonencounter_as_zero" = FALSE),
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
  "working_dir" = paste0(folder, "results/"), 
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
  "getsd" = TRUE,
  
  ## Model tuning (increase number of newtonsteps if final gradient is too large)
  "newtonsteps" = 1,
  
  ## Covariate data
  "X1_formula"= ~ CPE,
  "X2_formula"= ~ CPE,
  "X1config_cp" <- as.matrix(2),
  "X2config_cp" <- as.matrix(2) ,
  "covariate_data" = covariate_data)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save results locally ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Save VAST fit object
saveRDS(object = fit, 
        file = paste0(folder, "results/", species_name, "_VASTfit.RDS"))

## Save diagnostics and other outputs
results <- FishStatsUtils::plot_results( 
  fit = fit, 
  working_dir = paste0(folder, "results/output_plots/"),
  zrange = c(-3, 3), 
  n_cells = 600, 
  strata_names = c("Both", "EBS", "NBS"), 
  check_residuals = TRUE,
  n_samples = 0)

saveRDS(object = results, 
        file = paste0(folder, "results/output_plots/diagnostics.RDS"))

