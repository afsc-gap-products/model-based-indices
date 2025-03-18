##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       VAST Run
## Authors:       Zack Oyafuso (zack.oyafuso@noaa.gov)
##                Megsie siple (margaret.siple@noaa.gov)
## Assessor:      Pete Hulson (pete.hulson@noaa.gov)
##
## Notes:         R Version should be at least 4.0.2
##                2023 TOR: 
##     
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Package version checks ----
##   Updated every year
##   2023 TOR is in this google doc: 
##   https://docs.google.com/document/d/18CeXcHhHK48hrtkiC6zygXlHj6YVrWEd/edit
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
current_year <- 2023
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
##   Load VAST
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(VAST)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load Catch and Effort Data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_name <- "Gadus_macrocephalus"
Data_Geostat <- readRDS(file = paste0("species_specific_code/GOA/",
                                      species_name, "/data/Data_Geostat_", 
                                      species_name,".rds"))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import interpolation grid from repo
##   grid area is in m^2 and is converted to km^2
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GOAgrid <- read.csv(file= "extrapolation_grids/GOAThorsonGrid_Less700m.csv")
input_grid <- with(GOAgrid,
                   data.frame(Lat = ShapeLatitude,
                              Lon = ShapeLongitude,
                              Area_km2 = Shape_Area/1000000))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Configure VAST Settings
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
settings <- FishStatsUtils::make_settings(
  Version = VAST_cpp_version,
  Region = "User",
  purpose = "index2", 
  n_x = 750,
  ObsModel= c(2,1), 
  knot_method = "grid", 
  strata.limits = data.frame(STRATA = as.factor('All_areas')),
  fine_scale = TRUE,
  max_cells = 2000,
  bias.correct = TRUE,
  use_anisotropy = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Run model
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!dir.exists(paste0("species_specific_code/GOA/",
                       species_name, "/results"))) {
  dir.create(path = paste0("species_specific_code/GOA/",
                           species_name, "/results"), 
             recursive = TRUE)
}

result_dir <- paste0("species_specific_code/GOA/", species_name, "/results")

fit <- FishStatsUtils::fit_model( 
  "settings" = settings, 
  "Lat_i" = Data_Geostat[, "Lat"], 
  "Lon_i" = Data_Geostat[, "Lon"], 
  "t_i" = Data_Geostat[, "Year"], 
  "b_i" = Data_Geostat[, "Catch_KG"], 
  "a_i" = Data_Geostat[, "AreaSwept_km2"], 
  "v_i" = Data_Geostat[, "Vessel"], ## was ok to leave in because it"s all 
                                    ## "missing" or zero, so no vessel effects
  "refine" = TRUE,
  "Npool" = 100,
  "input_grid" = input_grid, 
  "optimize_args" = list("lower" = -Inf, "upper" = Inf),
  "working_dir" = paste0(getwd(), "/", result_dir, "/"))#,
  # "use_new_epsilon"  = F) #fixes error

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Get model outputs and diagnostics
##   Save VAST model object as RDS file
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(fit, 
     working_dir = paste0(getwd(), "/", result_dir),
     check_residuals = TRUE)

saveRDS(object = fit,
        file = paste0(getwd(), "/", result_dir, "/", 
                      species_name, "_VASTfit.RDS"))

