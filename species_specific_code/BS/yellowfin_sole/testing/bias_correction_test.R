##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Model-based estimate of age composition for BS YFS
## Authors:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Test sensitivity of age comps to bias correction option
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(VAST) 

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
##   Set species and output directory ----
##   Record session information
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
speciesName <- "yellowfin_sole"
workDir <- paste0("species_specific_code/BS/", speciesName, "/")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import data ----
##   Reduce data to just ages 1 - 4 and yars 2006 - 2009
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
Data_Geostat <- readRDS(paste0(workDir, 
                               "hindcast/data/data_geostat_agecomps.RDS"))
Data_Geostat <- subset(x = Data_Geostat, 
                       subset = Year %in% 2006:2009)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set Initial VAST settings ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
settings <- FishStatsUtils::make_settings(
  n_x = 50,
  Region = c("eastern_bering_sea"),
  purpose = "index2",
  fine_scale = TRUE,
  strata.limits = data.frame(STRATA = as.factor('All_areas')),
  ObsModel = c(2, 1),
  FieldConfig = c("Omega1" = "IID", "Epsilon1" = "IID", 
                  "Omega2" = "IID", "Epsilon2" = "IID"),
  RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0),
  OverdispersionConfig = c("Eta1" = 0, "Eta2" = 0),
  use_anisotropy = FALSE,
  Version = VAST_cpp_version,
  max_cells = 2000,
  knot_method = "grid",
  bias.correct = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Loop over bias.correct options (either T or F), create new result directory
##   change settings, run VAST, save outputs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
for (ibias in c(T, F)) {
  
  ## Set result directory
  res_dir <- 
    paste0(workDir, "/hindcast/test_agecomp_bias_correct_", ibias, "/")
  if(!dir.exists(res_dir)) dir.create(res_dir)
  
  ## Record package versions
  sink(paste0(res_dir, "/session_info.txt"), type = "output")
  sessionInfo()
  sink()
  
  ## Set bias correction option
  settings$bias.correct <- ibias
  
  ##   Run VAST model
  fit <- FishStatsUtils::fit_model(
    
    ## Input settings
    "working_dir" = paste0(getwd(), "/", res_dir),
    "settings" = settings, 
    "create_strata_per_region" = TRUE,
    
    ## Input data
    "Lat_i" = Data_Geostat[, "Lat"], 
    "Lon_i" = Data_Geostat[, "Lon"], 
    "t_i" = Data_Geostat[, "Year"],  
    "c_iz" = Data_Geostat[, "Age"] - 1, 
    "b_i" = Data_Geostat[, "Catch_KG"], 
    "a_i" = Data_Geostat[, "AreaSwept_km2"], 
    "v_i" = Data_Geostat[, "Vessel"],
    
    ## Model tuning
    "refine" = TRUE,
    "newtonsteps" = 2, 
    "test_fit" = FALSE)
  
  ## Save results locally
  saveRDS(object = fit, 
          file = paste0(getwd(), "/", res_dir, "VASTfit.RDS"))
  
  ## General output plots, DHARMa residuals
  results <- FishStatsUtils::plot_results( 
    fit = fit, 
    working_dir = paste0(getwd(), "/", res_dir, "output_plots/"),
    plot_set = NULL,
    check_residuals = TRUE)
  saveRDS(object = results, 
          file = paste0(res_dir, "output_plots/diagnostics.RDS"))
  
  ## Proportions
  Year_Set <- fit$year_labels
  Years2Include <- which(fit$year_labels != 2020)
  proportions <- FishStatsUtils::calculate_proportion( 
    TmbData = fit$data_list, 
    Index = results$Index, 
    year_labels = Year_Set, 
    years_to_plot = Years2Include,
    DirName = paste0(getwd(), "/", res_dir,  "output_plots/"))
  
  prop <- t(data.frame(proportions$Prop_ctl))
  colnames(prop) <- c(paste0("age_", seq(from = 1, length = ncol(prop) - 1 )),
                      paste0("age_", ncol(prop), "+"))
  rownames(prop) <- Year_Set
  
  if(!dir.exists(paste0(res_dir, "proportions/"))){
    dir.create(paste0(res_dir, "proportions/"))
  }
  saveRDS(object = proportions, 
          file = paste0(res_dir, "proportions/VAST_proportions.RDS"))
  write.csv(x = prop, 
            file = paste0(res_dir, 
                          "proportions/clean_proportions.csv"))
  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Compare comps ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
prop_bias_F <- read.csv(paste0(workDir, "hindcast/", 
                               "test_agecomp_bias_correct_FALSE/",
                               "proportions/clean_proportions.csv"))

prop_bias_T <- read.csv(paste0(workDir, "hindcast/", 
                               "test_agecomp_bias_correct_TRUE/",
                               "proportions/clean_proportions.csv"))

round(prop_bias_F[, -1], 3)
round(prop_bias_T[, -1], 3)
