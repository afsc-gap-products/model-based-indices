##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Model-based estimate of age composition for BS YFS
## Authors:       Zack Oyafuso (zack.oyafuso#noaa.gov)
##                Emily Markowitz (emily.markowitz#noaa.gov)
##                Jason Conner (jason.conner@noaa.gov)
##                James Thorson (james.thorson@noaa.gov)
## Output POC:    Ingrid Spies (ingrid.spies@noaa.gov)
## Description:   VAST estimates of age composition for EBS+NBS yellowfin sole
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
which_model <- c("hindcast", "production")[1]
species_name <- "yellowfin_sole"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(VAST) 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   System preferences ----
##   Updated every year
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
current_year <- 2024
VAST_cpp_version <- "VAST_v14_0_1"
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
if(!dir.exists(paste0(folder, "results_age/"))) 
  dir.create(path = paste0(folder, "results_age/"), recursive = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set species and output directory ----
##   Record session information
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
sink(file = paste0(folder, "results_age/session_info.txt"), 
     type = "output")
sessionInfo()
sink()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
Data_Geostat <- readRDS(file = paste0(folder, 
                                      "data/data_geostat_agecomps.RDS"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set VAST settings ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
settings <- FishStatsUtils::make_settings( 
  n_x = 50,
  Region = c("northern_bering_sea", "eastern_bering_sea"),
  purpose = "index2",
  fine_scale = TRUE,
  strata.limits = data.frame(STRATA = as.factor('All_areas')),
  ObsModel = c(2, 1),
  FieldConfig = c("Omega1" = "IID", "Epsilon1" = "IID", 
                  "Omega2" = "IID", "Epsilon2" = "IID"),
  RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 4, "Epsilon2" = 4),
  OverdispersionConfig = c("Eta1" = 0, "Eta2" = 0),
  use_anisotropy = FALSE,
  Version = VAST_cpp_version,
  max_cells = 2000,
  knot_method = "grid",
  bias.correct = TRUE)

strata_names <- c("Both", "EBS", "NBS")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Run VAST model ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
fit <- FishStatsUtils::fit_model(
  
  ## Input settings
  "working_dir" = paste0(folder, "results_age/"),
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
  "Npool" = 100,
  "refine" = TRUE,
  "newtonsteps" = 1, 
  "test_fit" = FALSE#, 
  # "startpar" = parameter_estimates$par,
  # "use_new_epsilon"  = FALSE,
  # "bias.correct.control" = list(nsplit = 20)
  )

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save results locally ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## VAST fit object
saveRDS(object = fit, 
        file = paste0(folder, "results_age/VASTfit.RDS"))

## Parameter estimates
saveRDS(object = fit$ParHat, 
        file = paste0(folder, "results_age/starting_parameters.RDS"))

## General output plots, DHARMa residuals
results <- FishStatsUtils::plot_results( 
  fit = fit, 
  working_dir = paste0(folder, "results_age/"),
  plot_set = NULL,
  strata_names = strata_names, 
  check_residuals = TRUE)

saveRDS(object = results, 
        file = paste0(folder, "results_age/VASTresults.RDS"))

## Mapping information
map_list = FishStatsUtils::make_map_info( 
  "Region" = settings$Region, 
  "spatial_list" = fit$spatial_list, 
  "Extrapolation_List" = fit$extrapolation_list)

## Predicted density maps for each age group across years
n_ages <- length(unique(Data_Geostat$Age))
Year_Set <- fit$year_labels

if(!dir.exists(paste0(folder, "results_age/predicted_density/"))){
  dir.create(paste0(folder, "results_age/predicted_density/"))
}

FishStatsUtils::plot_maps( 
  fit = fit, 
  plot_set = 3,
  n_cells = 200,
  category_names = c(paste0("Age ", 1:(n_ages-1)), 
                     paste0("Age ", n_ages, "+")),
  Obj = fit$tmb_list$Obj, 
  year_labels = Year_Set,
  PlotDF = map_list[["PlotDF"]],
  working_dir = paste0(folder, "results_age/predicted_density/")) 

## Predicted Spatial Fields
if(!dir.exists(paste0(folder, "results_age/spatial_effects/"))){
  dir.create(paste0(folder, "results_age/spatial_effects/"))
}
FishStatsUtils::plot_maps( 
  fit = fit, 
  plot_set = 16:17,
  n_cells = 200,
  category_names = c(paste0("Age ", 1:(n_ages-1)), 
                     paste0("Age ", n_ages, "+")),
  year_labels = Year_Set,
  Obj = fit$tmb_list$Obj, 
  PlotDF = map_list[["PlotDF"]],
  working_dir = paste0(folder, "results_age/spatial_effects/")) 

## Predicted Spatiotemporal Fields
if(!dir.exists(paste0(folder, "results_age/spatiotemporal_effects/"))){
  dir.create(paste0(folder, "results_age/spatiotemporal_effects/"))
}
FishStatsUtils::plot_maps( 
  fit = fit, 
  plot_set = 6:7,
  n_cells = 200,
  category_names = c(paste0("Age ", 1:(n_ages-1)), 
                     paste0("Age ", n_ages, "+")),
  year_labels = Year_Set,
  Obj = fit$tmb_list$Obj, 
  PlotDF = map_list[["PlotDF"]],
  working_dir = paste0(folder, "results_age/spatiotemporal_effects/")) 

## Predicted Proportions
proportions <- FishStatsUtils::calculate_proportion( 
  TmbData = fit$data_list, 
  Index = results$Index, 
  year_labels = Year_Set, 
  years_to_plot = which(fit$year_labels != 2020),
  strata_names = strata_names, 
  DirName = paste0(folder, "results_age/output_plots/"))

prop <- t(data.frame(proportions$Prop_ctl))
colnames(prop) <- c(paste0("age_", seq(from = 1, length = ncol(prop) - 1 )),
                    paste0("age_", ncol(prop), "+"))
rownames(prop) <- as.vector(sapply(X = strata_names, 
                                   FUN = function(x) paste0(x, "_", Year_Set)))

if(!dir.exists(paste0(folder, "results_age/proportions/"))){
  dir.create(paste0(folder, "results_age/proportions/"))
}
saveRDS(object = proportions, 
        file = paste0(folder, "results_age/proportionsVAST_proportions.RDS"))
write.csv(x = prop, 
          file = paste0(folder, "results_age/proportions/clean_proportions.csv"))
