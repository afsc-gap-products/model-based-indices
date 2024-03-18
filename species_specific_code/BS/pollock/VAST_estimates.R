# This code produces VAST index estimates for the Bering Sea Pollock assessment
# By: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.04.23
# Date updated: 2022.11.17


# # package installs --------------------------------------------------------
# 
# ### 2021 settings (updated packages): Rv4.0.2 VAST v3.6.1, FishStatsUtils v2.8.0, cpp VAST_v12_0_0
# ### 2022 settings: Rv4.0.2 or later VAST v3.8.2, FishStatsUtils v2.10.0, cpp VAST_v13_1_0, TMB v1.7.22, Matrix v1.4-0, DHARMa 0.4.5
# 
# 2022: Rv4.0.2 or later VAST v3.8.2, FishStatsUtils v2.10.0, cpp VAST_v13_1_0, TMB v1.7.22, Matrix v1.4-0, DHARMa 0.4.5
# ## check to ensure you're running the correct package versions
packageVersion('FishStatsUtils')
packageVersion('VAST')
packageVersion('INLA')
packageVersion('TMB')
packageVersion('TMBhelper')

# ## how to check your session info (see R version and all package versions)
sessionInfo()
# # R version 4.0.2 (2020-06-22)

# install vast ------------------------------------------------------------

# install.packages("TMB")
# library(TMB)

# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# library(INLA)
# devtools::install_github("james-thorson/FishStatsUtils", INSTALL_opts="--no-staged-install")
# library(FishStatsUtils)
# devtools::install_github("james-thorson/VAST@3.8.2", INSTALL_opts="--no-staged-install")
# library(VAST)

# libraries ---------------------------------------------------------------

library(here) 
library(tidyverse)
library(sf)
library(scales)
library(renv)
library(lubridate)
library(janitor)
library(VAST)

# Set R environment for consistency ---------------------------------------

# renv::init()
# renv::snapshot()
# renv::restore()

# Set species -------------------------------------------------------------

data_choice <- (readline(prompt = "Which data set do you want to use (Enter: a = EBS+NBS, e = EBS, n = NBS): "))
a

if(data_choice  == 'a') {
  region_select <- "EBS_NBS"
  use_region <- c("Eastern_Bering_Sea", "Northern_Bering_Sea")
  use_strata_names <- c("Both", "EBS", "NBS")
  folder <- "EBS-NBS"
} else if (data_choice  == 'e') {
  region_select <- "EBS"
  
  use_region <- c("Eastern_Bering_Sea")
  use_strata_names <- c("EBS")
  folder <- "EBS_only"
} else if (data_choice  == 'n') {
  region_select <- "NBS"
  
  use_region <- c("Northern_Bering_Sea")
  use_strata_names <- c("NBS")
  folder <- "NBS_only"
} else(print("Invalid selection: please select a for all EBS+NBS, e for EBS, or n for NBS"))

species <- 21740
this_year <- lubridate::year(today())
# this_year <- 2022  # set a different year for debugging
speciesName <- paste0("Walleye_Pollock_index_",this_year,"_",region_select)
workDir <- here::here("species_specific_code","BS","pollock")
dir.create(workDir, showWarnings = FALSE)

# Read in data ------------------------------------------------------------
index_data_all <- read.csv(here("species_specific_code","BS","pollock","data",paste0("VAST_ddc_all_", this_year, ".csv")))

# index_data_EBS <- read.csv(here("output", paste0("VAST_ddc_EBSonly_", this_year,".csv")))
# 
# index_data_NBS <- read.csv(here("output", paste0("VAST_ddc_NBSonly_", this_year,".csv")))

# Format catch data -------------------------------------------------------

if(data_choice  == 'a') {Data <- index_data_all
}else if(data_choice  == 'e') {Data <- index_data_EBS
}else if(data_choice  == 'n') {Data <- index_data_NBS
}else(print("Invalid selection: please select a for all EBS+NBS, e for EBS, or n for NBS"))

Data

# Format the data for VAST
Data_Geostat <-  dplyr::transmute(Data,
                                  Catch_KG_km2 = ddc_cpue_kg_ha*100, #  CPUE in kg/ha this converts it to kg/km^2
                                  Year = year,
                                  Vessel = "missing",
                                  AreaSwept_km2 = 1, # area swept is 1 when using CPUE instead of observed weight
                                  Lat = start_latitude,
                                  Lon = start_longitude,
                                  Pass = 0 )

unique(Data_Geostat$Vessel)
unique(Data_Geostat$AreaSwept_km2)
unique(Data_Geostat$Pass)

write_csv(Data_Geostat, file = here("output", paste0("VAST_datageostat_all_", this_year,".csv")))

Data_Geostat <- Data_Geostat %>% #dplyr::filter(Year <2020) %>%
  as.data.frame()

if(data_choice  == 'a') {Data_Geostat <- read.csv(here("output", paste0("VAST_datageostat_all_", this_year,".csv")))
}else if(data_choice  == 'e') {Data_Geostat <- read.csv()
}else if(data_choice  == 'n') {Data_Geostat <- read.csv()
}else(print("Invalid selection: please select a for all EBS+NBS, e for EBS, or n for NBS"))

head(Data_Geostat)

table(Data_Geostat$Year)

# Cold pool covariate -----------------------------------------------------

# devtools::install_github("afsc-gap-products/coldpool",  force = TRUE, R_REMOTES_NO_ERRORS_FROM_WARNINGS= TRUE)#, lib = .libPaths()[2])
coldpool:::cold_pool_index

cold_pool <- coldpool:::cold_pool_index %>%
  as_tibble() %>%
  clean_names #%>%

# cold_pool <- read_csv(here("output", "cold_pool.csv"))
head(cold_pool)

cp <- cold_pool %>%
  dplyr::select(-year, -last_update) %>%
  scale(center = TRUE, scale = TRUE) %>%
  as_tibble()
missing_yr <- bind_cols(year = 2020,
                area_lte2_km2 = 0,
                area_lte1_km2 = 0,
                area_lte0_km2 = 0,
                area_lteminus1_km2 = 0,
                mean_gear_temperature = 0,
                mean_bt_lt100m = 0,
                mean_surface_temperature = 0)

cp_data <- cold_pool %>%
  dplyr::select(year) %>%
  bind_cols(cp) %>%
  dplyr::filter(year %in% unique(Data_Geostat[,'Year'])) %>%
  full_join(missing_yr) %>%
  arrange(year)

covariate_data <- data.frame( "Year"=cp_data[,"year"],
                              "Lat"=mean(Data_Geostat[,'Lat']),
                              "Lon"=mean(Data_Geostat[,'Lon']),
                              "area_lte2_km2" = cp_data[,"area_lte2_km2"]) %>%
    rename(Year = year) %>%
    data.frame()
write.csv(covariate_data, here(workDir, "data", "cold_pool_scaled_formatted.csv"))

covariate_data <- read.csv(here(workDir, "data", "cold_pool_scaled_formatted.csv"))

# VAST Settings -----------------------------------------------------------
Version <- get_latest_version( package="VAST" )
# Version <- "VAST_v14_0_1"  #2023 not recognized?
# Version <- "VAST_v13_1_0" #2022
# Version <- "VAST_v12_0_0" #2021
Region <- use_region
strata_names <- use_strata_names
Method <- "Mesh"
knot_method <- "grid"
grid_size_km <- 25
# n_x <- 100   # Specify number of stations (a.k.a. "knots") ##CIA: test on 100, final run on 750
n_x <- 750
FieldConfig <- c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig <- c("Beta1"=0, "Beta2"=0, "Epsilon1"=4, "Epsilon2"=4)
OverdispersionConfig <- c("Eta1"=0, "Eta2"=0)
ObsModel <- c(2,1) #Need (2,4) if there are some years with 100% encounter rate
Options <-  c("Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE, "treat_nonencounter_as_zero"= FALSE) #TRUE )
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
  Version = Version,
  max_cells = max_cells,
  # bias.correct = FALSE, #turn off for full run
  knot_method = knot_method
  
)

# covariate formula
formula <- ~ area_lte2_km2
Xconfig_zcp <- array(2, dim=c(2,1,1) )
X1config_cp <- as.matrix(2)
X2config_cp <- as.matrix(2)

# read saved image (for troubleshooting) ----------------------------------

# load('pollock_index.RData')

# Build the model ---------------------------------------------------------
# quick fit:

# settings_quick <- settings
# settings_quick$n_x <- 100  # 100 knots for testing
# settings_quick$bias.correct <- TRUE

options(max.print = .Machine$integer.max)

# fit <- fit_model( "settings"=settings_quick, 
#                   "Lat_i"=Data_Geostat[,'Lat'], 
#                   "Lon_i"=Data_Geostat[,'Lon'], 
#                   "t_i"=Data_Geostat[,'Year'], 
#                   "c_i"=rep(0,nrow(Data_Geostat)),
#                   "b_i"=Data_Geostat[,'Catch_KG_km2'],
#                   "a_i"=Data_Geostat[,'AreaSwept_km2'], 
#                   "v_i"=Data_Geostat[,'Vessel'],
#                   "create_strata_per_region"=TRUE,
#                   #"getJointPrecision"=getJointPrecision, # turn on for full run
#                   #"getReportCovariance"=getReportCovariance,  # turn on for full run
#                   "X1_formula"=formula,
#                   "X2_formula"=formula,
#                   "X1config_cp" = X1config_cp,
#                   "X2config_cp" = X2config_cp ,
#                   "covariate_data"= covariate_data,
#                   # getsd=TRUE,
#                   getsd=FALSE,            #for testing
#                   # "run_model" = FALSE,    #for testing
#                   # "build_model" = FALSE,  #for testing
#                   # test_fit = FALSE,       #for testing
#                   newtonsteps=0,          #for testing
#                   # CheckForBugs = FALSE,   #for testing
#                   "working_dir" = workDir
# )
# 
# saveRDS(fit, file = here(workDir, "VAST Index" ,"VASTfit_quick.RDS"))
# fit_check <- readRDS(file = here(workDir, "VAST Index", "VASTfit_quick.RDS"))
# fit_check$ParHat

# full model fit:
start.time <- Sys.time() #"2024-03-14 16:49 PDT"
full_fit <- fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,'Lat'], 
                  "Lon_i"=Data_Geostat[,'Lon'], 
                  "t_i"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)),
                  "b_i"=Data_Geostat[,'Catch_KG_km2'],
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'],
                  #parameters=fit_check$ParHat,   #use params from quick run as starting point
                  "create_strata_per_region"=TRUE,
                  "getJointPrecision"=getJointPrecision, # turn on for full run
                  "getReportCovariance"=getReportCovariance,  # turn on for full run
                  "X1_formula"=formula,
                  "X2_formula"=formula,
                  "X1config_cp" = X1config_cp,
                  "X2config_cp" = X2config_cp ,
                  "covariate_data"= covariate_data,
                  getsd=TRUE,
                  # getsd=FALSE,            #for testing
                  # "run_model" = FALSE,    #for testing -- if an issue, try uncommenting this one
                  # "build_model" = FALSE,  #for testing
                  # test_fit = FALSE,       #for testing
                  newtonsteps=1,          #for testing
                  # CheckForBugs = FALSE,   #for testing
                  "working_dir" = workDir
)

# 2022: setting questions- "input_grid", "refine", optimize_args
stop.time <- Sys.time()

# Save results

saveRDS(full_fit, file = here(workDir, "results", "VAST Index", "VASTfit_full.RDS"))
# full_fit <- readRDS(file = paste0(workDir,"/VASTfit_full.RDS"))



# Plots -------------------------------------------------------------------
# 
# plot( full_fit)
# 
# plot( full_fit, zrange = c(-3,3), n_cells = 600, strata_names = strata_names )
# 
# 
# ## Save COG (center of gravity) for ESP request
# results = plot( full_fit, n_cells=200^2 )
# write.csv( results$Range$COG_Table, file=paste0(workDir, "/COG.csv"), row.names=FALSE )
# 
# ##save effective area occupied for ESP request
# report = TMB::summary.sdreport(full_fit$parameter_estimates$SD)
# ln_km2 = report[which(rownames(report)=="log_effective_area_ctl"),c('Estimate','Std. Error')]
# Year <- sort(unique(full_fit$year_labels))
# ln_km2 <- as.data.frame(cbind(ln_km2, Year))
# ln_km2 <- ln_km2[which(ln_km2$Year %in% unique(fit$data_frame$t_i)),]
# write.csv( ln_km2, file=paste0(workDir,"/ln_effective_area.csv"), row.names=FALSE )


# Plots 2022 ---------------------------------------------------------------

# If you need to load a fit in a new session:
dyn.load(dynlib("VAST_v12_0_0"))

# Record package versions
sink("session_info.txt", type = "output")
sessionInfo()
sink()

# Plot results
results <- plot_results( full_fit, 
                         zrange = c(-3,3),
                         n_cells = 600, 
                         strata_names = strata_names, 
                         check_residuals=TRUE,
                         n_samples=0)

# plot_results( full_fit, 
#               # zrange = c(-3,3), 
#               n_cells = 50^2, 
#               strata_names = strata_names)#, 
# check_residuals=TRUE,
# n_samples=0)

saveRDS(results, file = here(workDir, "results", "VAST Index", "VASTresults.RDS"))

map_list = make_map_info( "Region"=settings$Region, "spatial_list"=full_fit$spatial_list, "Extrapolation_List"=full_fit$extrapolation_list )
plot_maps(fit = full_fit, Obj=full_fit$tmb_list$Obj, PlotDF=map_list[["PlotDF"]] )


# ESP products
write.csv( results$Range$COG_Table, file="COG.csv", row.names=FALSE )

##save effective area occupied for ESP request
report = TMB::summary.sdreport(full_fit$parameter_estimates$SD)
ln_km2 = report[which(rownames(report)=="log_effective_area_ctl"),c('Estimate','Std. Error')]
Year <- sort(unique(full_fit$year_labels))
ln_km2 <- as.data.frame(cbind(ln_km2, Year))
ln_km2 <- ln_km2[which(ln_km2$Year %in% unique(full_fit$data_frame$t_i)),]
write.csv( ln_km2, file="ln_effective_area.csv", row.names=FALSE )
