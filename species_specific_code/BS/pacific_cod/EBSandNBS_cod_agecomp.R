library(dplyr)
library(VAST)
library(tictoc)

# Set species, model -------------------------------------------------------

which_model <- c("hindcast", "production")[2]
compare <- FALSE # If compare = TRUE, using prior year's alk
species <- 21720
species_name <- "pacific_cod"

workDir <- paste0(getwd(),"/species_specific_code/BS/", 
                  species_name, "/", which_model, "/")
if(!dir.exists(workDir))
  dir.create(path = workDir, recursive = TRUE)
if(!dir.exists(paste0(workDir, "results_age/")))
  dir.create(path = paste0(workDir, "results_age/"), recursive = TRUE)

# Record sessionInfo -------------------------------------------------------
sink(file = paste0(workDir, "results_age/session_info.txt"), 
     type = "output")
sessionInfo()
sink()

# Make sure package versions are correct for current year ------------------
current_year <- 2024
prev_year <- current_year-1
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

# Import data -------------------------------------------------------------
Data_Geostat <- readRDS(file = paste0(workDir, 
                                      "data/data_geostat_agecomps.RDS"))

# Settings ----------------------------------------------------------------
Region <- c("eastern_bering_sea","northern_bering_sea")
Method <- "Mesh"
knot_method <- "grid"
grid_size_km <- 25
n_x <- 50 # Specify number of "knots"
FieldConfig <- c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig <- c("Beta1"=0, "Beta2"=0, "Epsilon1"=4, "Epsilon2"=4)
OverdispersionConfig <- c("Eta1"=0, "Eta2"=0)
ObsModel <- c(2,4)
Options <- c("Calculate_Range" = FALSE, "Calculate_effective_area" = FALSE, "treat_nonencounter_as_zero" = TRUE )
Aniso <- FALSE
Npool <- 100 #20
fine_scale <- TRUE
BiasCorr <- TRUE
max_cells <- 2000

strata.limits <- data.frame('STRATA'="All_areas")

# Make settings 
settings <- make_settings( 
                      n_x = n_x,
                      Region = Region,
                      purpose = "index2",
                      fine_scale = fine_scale,
                      strata.limits = strata.limits,
                      ObsModel = ObsModel,
                      FieldConfig = FieldConfig,
                      RhoConfig = RhoConfig,
                      OverdispersionConfig = OverdispersionConfig,
                      Options = Options,
                      use_anisotropy = Aniso,
                      Version = VAST_cpp_version,
                      max_cells = max_cells,
                      knot_method = knot_method,
                      bias.correct = BiasCorr
)

strata_names = c("Both","EBS","NBS")

# Fit model ------------------------------------------------------------

tic("running model")
fit = fit_model( "settings"=settings, 
                 "Lat_i"=Data_Geostat[,'Lat'], 
                 "Lon_i"=Data_Geostat[,'Lon'], 
                 "t_i"=Data_Geostat[,'Year'],  # "t_i"=rep(2019,nrow(Data_Geostat)),
                 "c_i"=Data_Geostat[,'Age'], 
                 "b_i"=Data_Geostat[,'Catch_KG'], 
                 "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                 "v_i"=Data_Geostat[,'Vessel'],
                 Npool = Npool, 
                 test_fit=FALSE, # set to FALSE if want to avoid interruption due to marginal fit issues
                 create_strata_per_region=TRUE,
                 "working_dir" = paste0(workDir,"/results_age/"),
                 "CompileDir" = paste0(workDir,"/results_age/")
)
toc()  

# Save results
saveRDS(fit, file = paste0(workDir,"/results_age/",species_name,"_VASTfit.RDS"))

## Parameter estimates
saveRDS(object = fit$ParHat, 
        file = paste0(workDir, "results_age/starting_parameters.RDS"))

## General output plots, DHARMa residuals
results <- FishStatsUtils::plot_results( 
  fit = fit, 
  working_dir = paste0(workDir, "results_age/"),
  plot_set = NULL,
  strata_names = strata_names, 
  check_residuals = TRUE)

saveRDS(object = results, 
        file = paste0(workDir, "results_age/VASTresults.RDS"))

## Mapping information
map_list = FishStatsUtils::make_map_info( 
  "Region" = settings$Region, 
  "spatial_list" = fit$spatial_list, 
  "Extrapolation_List" = fit$extrapolation_list)

## Predicted density maps for each age group across years
n_ages <- length(unique(Data_Geostat$Age))
Year_Set <- fit$year_labels

if(!dir.exists(paste0(workDir, "results_age/predicted_density/"))){
  dir.create(paste0(workDir, "results_age/predicted_density/"))
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
  working_dir = paste0(workDir, "results_age/predicted_density/")) 

## Predicted Spatial Fields
if(!dir.exists(paste0(workDir, "results_age/spatial_effects/"))){
  dir.create(paste0(workDir, "results_age/spatial_effects/"))
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
  working_dir = paste0(workDir, "results_age/spatial_effects/")) 

## Predicted Spatiotemporal Fields
if(!dir.exists(paste0(workDir, "results_age/spatiotemporal_effects/"))){
  dir.create(paste0(workDir, "results_age/spatiotemporal_effects/"))
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
  working_dir = paste0(workDir, "results_age/spatiotemporal_effects/")) 

## Predicted Proportions
if(!dir.exists(paste0(workDir, "results_age/proportions/"))){
  dir.create(paste0(workDir, "results_age/proportions/"))
}

proportions <- FishStatsUtils::calculate_proportion( 
  TmbData = fit$data_list, 
  Index = results$Index, 
  year_labels = Year_Set, 
  years_to_plot = which(fit$year_labels != 2020),
  strata_names = strata_names, 
  DirName = paste0(workDir, "results_age/proportions/"))

prop <- t(data.frame(proportions$Prop_ctl))
colnames(prop) <- c(paste0("age_", seq(from = 1, length = ncol(prop) - 1 )),
                    paste0("age_", ncol(prop), "+"))
rownames(prop) <- as.vector(sapply(X = strata_names, 
                                   FUN = function(x) paste0(x, "_", Year_Set)))

saveRDS(object = proportions, 
        file = paste0(workDir, "results_age/proportionsVAST_proportions.RDS"))
write.csv(x = prop, 
          file = paste0(workDir, "results_age/proportions/clean_proportions.csv"))