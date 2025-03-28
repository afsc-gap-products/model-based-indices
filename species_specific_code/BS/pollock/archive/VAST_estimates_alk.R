# This code produces VAST agecomp estimates for the Bering Sea Pollock assessment
# By: Caitlin I. Allen Akselrud; modified from O'Leary base code 2022 and Thorson 2019 code
# Maintained by: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2021.04.23
# Date updated: 2024.09.06

# Notes -------------------------------------------------------------------

### 2022 settings: Rv4.0.2 or later VAST v3.8.2, FishStatsUtils v2.10.0, cpp VAST_v13_1_0, TMB v1.7.22, Matrix v1.4-0, DHARMa 0.4.5

# libraries ---------------------------------------------------------------

library(tidyverse)
library(rgdal)
library(VAST)
library(lubridate)
library(here)
library(FishStatsUtils)
library(INLA)
library(TMB)

# Set R environment for consistency ---------------------------------------

# renv::init()
# renv::snapshot()
# renv::restore()

# Set species -------------------------------------------------------------

species <- 21740
this_year <- lubridate::year(today())
# this_year <- 2022  # different year for debugging
Species = "pollock"
speciesName <- paste0("Walleye_Pollock_age_",lubridate::year(today()),"_EBS-NBS")
workDir <- here::here("species_specific_code","BS", Species)
# dir.create(workDir, showWarnings = FALSE)

# Settings ----------------------------------------------------------------

Version <- get_latest_version( package="VAST" )
# Version <- "VAST_v14_0_1"
Region <- c("Eastern_Bering_Sea","Northern_Bering_Sea")
knot_method <- "grid"
grid_size_km = 25
n_x = 50   # Specify number of stations (a.k.a. "knots")
FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=4, "Epsilon2"=4)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
ObsModel = c(2,1)
Options =  c("Calculate_Range"=FALSE, "Calculate_effective_area"=FALSE, "treat_nonencounter_as_zero"=TRUE )
Aniso = TRUE
Npool = 100
fine_scale <- TRUE
BiasCorr = TRUE
max_cells <- 2000

strata.limits <- data.frame('STRATA'="All_areas")


# Make settings 
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
  knot_method = knot_method,
  bias.correct = BiasCorr
)

strata_names = c("Both","EBS","NBS")

# Age Composition ---------------------------------------------------------

# pollock data ------------------------------------------------------------

Data_Geostat <- read.csv(here("species_specific_code","BS",Species,"data",paste0("VAST_ddc_alk_", this_year, ".csv")))
Data_Geostat$AreaSwept_km2 <- 1  # set to 1 b/c we're using CPUE, not catch

# check for sample size
table(Data_Geostat$Age)


# Run Analysis ------------------------------------------------------------

#### Explore the data ####

Date = Sys.Date()
RunDir = paste0(workDir,"/results","/Comps/")
dir.create(RunDir, recursive = TRUE)
setwd(RunDir)

# Run model
start.time <- Sys.time() 
fit = fit_model( "settings"=settings, 
                 "Lat_i"=Data_Geostat[,'Lat'], 
                 "Lon_i"=Data_Geostat[,'Lon'], 
                 "t_i"=Data_Geostat[,'Year'],  # "t_i"=rep(2019,nrow(Data_Geostat)),
                 "c_iz"=Data_Geostat[,'Age'] - 1,  # need to change this so index starts at 0
                 "b_i"=Data_Geostat[,'CPUE_num'],  # SNW: updated column name in 2024 - actually numbers, not kg
                 # "b_i"=as_units(Data_Geostat[,'Catch_KG'], "count"), #new for 2023 changes
                 "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                 # "v_i"=Data_Geostat[,'Vessel'],  # not using vessel data
                 Npool = Npool, 
                 test_fit=F,
                 newtonsteps = 0,       # for testing
                 # getsd = FALSE,         # for testing
                 # test_fit = FALSE,      # for testing
                 # "run_model" = FALSE,   # for testing
                 # "build_model" = FALSE, # for testing
                 # CheckForBugs = FALSE,  # for testing
                 create_strata_per_region=TRUE)
stop.time <- Sys.time()

# Check gradient
fit$parameter_estimates$max_gradient

# Save results
dir.create(here(workDir, "results", "Comps"))
saveRDS(fit, file = here(workDir, "results", "Comps", "VASTfit_age.RDS"))

#Load results if using a previous model run
#fit <- readRDS(file = paste0(workDir,"/VASTfit_age.RDS"))

# Plots -------------------------------------------------------------------
# If you need to load a fit in a new session:
# dyn.load(dynlib("VAST_v12_0_0"))
# load(here("VAST_results", "Walleye_Pollock_age_2022_EBS-NBS", "Comps_2022-06-28_Walleye Pollock Agecomp_npool=100_BiasCorr=TRUE","VASTresults.RDS"))
# fit2 <- readRDS(file = paste0(workDir,"/VASTfit_age.RDS"))

# Record package versions
sink("session_info.txt", type = "output")
sessionInfo()
sink()

# Plot results
results <- plot_results(fit, #zrange = c(-3,3), n_cells = 600, 
                        strata_names = strata_names, 
                        check_residuals = TRUE)
saveRDS(results, file = here(workDir, "results", "Comps", "VASTresults_age.RDS"))

#Load results if taking a previous model run
# readRDS(file = paste0(getwd(), "/VASTresults_age.RDS"))

# Expand to proportional population numbers -------------------------------
VASTfit <- fit

## Variable names have changed in FishStatsUtils check for consistency ?calculate_proportion
Year_Set =VASTfit$year_labels
Years2Include = which(VASTfit$year_labels != 2020)
proportions = calculate_proportion( TmbData=VASTfit$data_list, 
                                    Index=results$Index, 
                                    year_labels = Year_Set, 
                                    years_to_plot = Years2Include,
                                    strata_names=strata_names)

prop <- data.frame(t(data.frame(proportions$Prop_ctl))) %>%
    # drop_na() %>%
    dplyr::rename("age_1"=1,"age_2"=2,"age_3"=3,"age_4"=4,"age_5"=5,"age_6"=6,"age_7"=7,"age_8"=8,"age_9"=9,"age_10"=10,
           "age_11"=11,"age_12"=12,"age_13"=13,"age_14"=14,"age_15+"=15) %>%
    dplyr::bind_cols(data.frame(Year = rep(Year_Set,3) ), 
              data.frame(Region = c(rep(strata_names[1],length(Year_Set)),
                                    rep(strata_names[2],length(Year_Set)),
                                    rep(strata_names[3],length(Year_Set)) ) 
              )
    )

write.csv(prop, file = here(workDir, "results", "Comps", "proportions.csv"), row.names = FALSE)
