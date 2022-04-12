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
finalanalysis <- F # this will make the model work faster while troubleshooting

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(coldpool)
library(googledrive)
library(tidyverse)
library(tibble)
library(VAST) # VAST 3.6.1, # devtools::install_github('james-thorson/VAST@3.8.2', INSTALL_opts='--no-staged-install')
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
##   Set species and output directory ----
##   Record session information
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
species <- 10210
speciesName <- "yellowfin_sole"
workDir <- paste0("species_specific_code/BS/", speciesName, "/")

## Record package versions
sink(paste0(workDir, "/hindcast/results_age/session_info.txt"), 
     type = "output")
sessionInfo()
sink()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
Data_Geostat <- readRDS(paste0(workDir, 
                               "hindcast/data/data_geostat_agecomps.RDS"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set VAST settings ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
settings <- FishStatsUtils::make_settings( 
  n_x = 50,
  Region = c("northern_bering_sea", "eastern_bering_sea"),
  purpose = "index2",
  fine_scale = ifelse(test = finalanalysis, yes = TRUE, no = FALSE),
  strata.limits = data.frame(STRATA = as.factor('All_areas')),
  ObsModel = c(2, 1),
  FieldConfig = c("Omega1" = "IID", "Epsilon1" = "IID", 
                  "Omega2" = "IID", "Epsilon2" = "IID"),
  RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 4, "Epsilon2" = 4),
  OverdispersionConfig = c("Eta1" = 0, "Eta2" = 0),
  use_anisotropy = TRUE,
  Version = VAST_cpp_version,
  max_cells = 2000,
  knot_method = "grid",
  bias.correct = ifelse(test = finalanalysis, yes = TRUE, no = FALSE))

strata_names <- c("Both", "EBS", "NBS")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Turn off those year-age categories that have 0% encounters
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
Data_Geostat$Age <- Data_Geostat$Age - 1
mean_cpue <- aggregate(formula = Catch_KG ~ Age + Year, 
                       data = Data_Geostat, 
                       FUN = mean)
zero_enc <- subset(x = mean_cpue, subset = Catch_KG == 0)
if (nrow(zero_enc) > 0) {
  remove_idx <- c()
  for (irow in 1:nrow(zero_enc)) {
    remove_idx <- c(remove_idx, with(Data_Geostat, 
                                     which(Year == zero_enc$Year[irow] &
                                             Age == zero_enc$Age[irow])))
  }
  # Data_Geostat[remove_idx, c("Catch_KG")] <- NA
  Data_Geostat <- Data_Geostat[-remove_idx, ]
}

rm(mean_cpue, zero_enc, irow, remove_idx)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Run VAST model ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# Data_Geostat <- subset(x = Data_Geostat)
fit <- FishStatsUtils::fit_model(
  "settings" = settings, 
  "Lat_i" = Data_Geostat[, "Lat"], 
  "Lon_i" = Data_Geostat[, "Lon"], 
  "t_i" = Data_Geostat[, "Year"],  
  "c_iz" = Data_Geostat[, "Age"], 
  "b_i" = Data_Geostat[, "Catch_KG"], 
  "a_i" = Data_Geostat[, "AreaSwept_km2"], 
  "v_i" = Data_Geostat[, "Vessel"],
  "Npool" = 100,
  "refine" = T,
  "test_fit" = F, 
  "working_dir" = paste0(getwd(), "/", workDir, "/hindcast/results_age/"),
  create_strata_per_region = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save fit ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
saveRDS(object = fit, 
        file = paste0(workDir, "/hindcast/results_age/VASTfit.RDS"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Output plots ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
## General output plots, DHARMa residuals
results <- FishStatsUtils::plot_results( 
  fit = fit, 
  working_dir = paste0(getwd(),  "/species_specific_code/BS/",
                       speciesName, "/hindcast/results_age/output_plots/"),
  plot_set = NULL,
  strata_names = strata_names, 
  check_residuals = TRUE)
saveRDS(object = results, 
        file = paste0("species_specific_code/BS/", speciesName, 
                      "/hindcast/results_age/output_plots/diagnostics.RDS"))

## Mapping information
map_list = FishStatsUtils::make_map_info( 
  "Region" = settings$Region, 
  "spatial_list" = fit$spatial_list, 
  "Extrapolation_List" = fit$extrapolation_list)

##   Predicted density maps for each age group across years
n_ages <- length(unique(Data_Geostat$Age)) 
if(!dir.exists(paste0(workDir, "/hindcast/results_age/predicted_density/"))){
  dir.create(paste0(workDir, "/hindcast/results_age/predicted_density/"))
}
FishStatsUtils::plot_maps( 
  fit = fit, 
  plot_set = 3,
  category_names = c(paste0("Age ", 1:(n_ages-1)), 
                     paste0("Age ", n_ages, "+")),
  Obj = fit$tmb_list$Obj, 
  PlotDF = map_list[["PlotDF"]],
  working_dir = paste0(workDir, "/hindcast/results_age/predicted_density/")) 

## Predicted Spatial Fields
if(!dir.exists(paste0(workDir, "/hindcast/results_age/spatial_effects/"))){
  dir.create(paste0(workDir, "/hindcast/results_age/spatial_effects/"))
}
FishStatsUtils::plot_maps( 
  fit = fit, 
  plot_set = 16:17,
  category_names = c(paste0("Age ", 1:(n_ages-1)), 
                     paste0("Age ", n_ages, "+")),
  Obj = fit$tmb_list$Obj, 
  PlotDF = map_list[["PlotDF"]],
  working_dir = paste0(workDir, "/hindcast/results_age/spatial_effects/")) 

## Predicted Spatiotemporal Fields
if(!dir.exists(paste0(workDir, 
                      "/hindcast/results_age/spatiotemporal_effects/"))){
  dir.create(paste0(workDir, "/hindcast/results_age/spatiotemporal_effects/"))
}
FishStatsUtils::plot_maps( 
  fit = fit, 
  plot_set = 6:7,
  category_names = c(paste0("Age ", 1:(n_ages-1)), 
                     paste0("Age ", n_ages, "+")),
  Obj = fit$tmb_list$Obj, 
  PlotDF = map_list[["PlotDF"]],
  working_dir = paste0(workDir, "/hindcast/results_age/",
                       "spatiotemporal_effects/")) 

## Proportions
Year_Set <- fit$year_labels
Years2Include <- which(fit$year_labels != 2020)
proportions <- FishStatsUtils::calculate_proportion( 
  TmbData = fit$data_list, 
  Index = results$Index, 
  year_labels = Year_Set, 
  years_to_plot = Years2Include,
  strata_names = strata_names, 
  DirName = paste0(getwd(),  "/species_specific_code/BS/",
                   speciesName, "/hindcast/results_age/output_plots/"))

prop <- t(data.frame(proportions$Prop_ctl))
colnames(prop) <- c(paste0("age_", seq(from = 1, length = ncol(prop) - 1 )),
                    paste0("age_", ncol(prop), "+"))
rownames(prop) <- as.vector(sapply(X = strata_names, 
                                   FUN = function(x) paste0(x, "_", Year_Set)))

if(!dir.exists(paste0(workDir, "/hindcast/results_age/proportions/"))){
  dir.create(paste0(workDir, "/hindcast/results_age/proportions/"))
}
saveRDS(object = proportions, 
        file = paste0(workDir, 
                      "/hindcast/results_age/proportions/",
                      "VAST_proportions.RDS"))
write.csv(x = prop, 
          file = paste0(workDir, 
                        "/hindcast/results_age/proportions/",
                        "clean_proportions.csv"))
