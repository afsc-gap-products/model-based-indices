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
remotes::install_github("afsc-gap-products/coldpool")
library(tidyverse)
library(VAST) 
library(sf)
library(scales)
library(coldpool)
library(akgfmaps)

# Set species, model -------------------------------------------------------
which_model <- c("hindcast", "production")[1] # specify by changing index
species <- 21720
species_name <- "kamchatka_flounder"

workDir <- paste0(getwd(),"/species_specific_code/BS/", 
                  species_name, "/", which_model, "/")
if(!dir.exists(workDir))
  dir.create(path = workDir, recursive = TRUE)
if(!dir.exists(paste0(workDir, "results/")))
  dir.create(path = paste0(workDir, "results/"), recursive = TRUE)

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
##   Set VAST settings ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Region <- "Eastern_Bering_Sea"
strata_names = "EBS"
knot_method <- "grid"
grid_size_km <- 25
n_x <- 350 # Reduced from standard of 750 knots used in EBS + NBS indices
FieldConfig <- c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig <- c("Beta1"=0, "Beta2"=0, "Epsilon1"=4, "Epsilon2"=4)
OverdispersionConfig <- c("Eta1"=0, "Eta2"=0)
ObsModel <- c(2,1)
Options <-  c("Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE, "treat_nonencounter_as_zero"=FALSE )
Aniso <- TRUE
BiasCorr <- TRUE
getJointPrecision <- TRUE
getReportCovariance <- TRUE
fine_scale <- TRUE
max_cells <- 2000
strata.limits <- data.frame(STRATA = as.factor('All_areas'))


settings <- make_settings( 
  n_x=n_x,
  Region=Region,
  purpose="index2",
  fine_scale = fine_scale,
  strata.limits=strata.limits,
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

# Read catch data ---------------------------------------------------------

dat <- read_rds(paste0(workDir, "data/data_geostat_biomass_index.RDS"))

# Cold Pool covariate -----------------------------------------------------
cpe <- scale(coldpool:::cold_pool_index$AREA_LTE2_KM2)
covariate_data <- data.frame(Year = c(coldpool:::cold_pool_index$YEAR, 2020),
                             Lat = mean(dat$Lat),
                             Lon = mean(dat$Lon), 
                             cpe = c(cpe, 0))

# Load covariates
formula <- ~ cpe
Xconfig_zcp <- array(2, dim=c(2,1,1) )
X1config_cp <- as.matrix(2)
X2config_cp <- as.matrix(2)


# Build the model ---------------------------------------------------------
fit <- fit_model( "settings"=settings, 
                  "Lat_i"=dat[,'Lat'], 
                  "Lon_i"=dat[,'Lon'], 
                  "t_i"=dat[,'Year'], 
                  "c_i"=rep(0,nrow(dat)), 
                  "b_i"=dat[,'Catch_KG'], 
                  "a_i"=dat[,'AreaSwept_km2'], 
                  "v_i"=dat[,'Vessel'],
                  "create_strata_per_region"=TRUE,
                  "getJointPrecision"=getJointPrecision, 
                  "getReportCovariance"=getReportCovariance,
                  "X1_formula"=formula,
                  "X2_formula"=formula,
                  "X1config_cp" <- X1config_cp,
                  "X2config_cp" <- X2config_cp,
                  "covariate_data"=covariate_data,
                  "test_fit" = TRUE,
                  "working_dir" = workDir,
                  "newtonsteps" = 1
)


# Save results

saveRDS(fit, file = paste0(workDir, "results/VASTfit.RDS"))


# Plots -------------------------------------------------------------------
# If you need to load a fit in a new session:
#dyn.load(dynlib(VAST_cpp_version))

# Plot results
results <- plot_results( fit, 
                         zrange = c(-3,3), 
                         n_cells = 600, 
                         strata_names = strata_names, 
                         check_residuals=TRUE,
                         n_samples=0,
                         working_dir = paste0(workDir, "results/")
)

saveRDS(results, file = paste0(workDir, "results/VASTresults.RDS"))

# If residual plots don't... uh... plot...
# plot_quantile_residuals( fit=fit )
# 
# map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
# plot_maps(fit=fit, n_cells = 600, Obj=fit$tmb_list$Obj, PlotDF=map_list[["PlotDF"]] )

# Plot DB vs VAST Comparison ----------------------------------------------
plot(dat %>% group_by(Year) %>% summarise(mean(Catch_KG)), type ="l") # mean cpue for quick look