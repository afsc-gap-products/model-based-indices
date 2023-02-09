# This code produces VAST agecomp estimates for the Bering Sea Pollock assessment
# This code is converted from example code by Jason Connor/Jim Thorson
# Modifided by: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.04.23
# Date updated: 2022.04.19

# Notes -------------------------------------------------------------------

### 2022 settings: Rv4.0.2 or later VAST v3.8.2, FishStatsUtils v2.10.0, cpp VAST_v13_1_0, TMB v1.7.22, Matrix v1.4-0, DHARMa 0.4.5

# libraries ---------------------------------------------------------------

library(googledrive)
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
Species = "Walleye Pollock Agecomp"
speciesName <- paste0("Walleye_Pollock_age_",lubridate::year(today()),"_EBS-NBS")
workDir <- here::here("VAST_results", speciesName)
# RootDir = workDir
dir.create(workDir, showWarnings = FALSE)
# dir.create(workDir)
# dir.create(here::here("data","shapefiles"))
# dir.create(here::here("output"))
# setwd(workDir)


# Get data from Google drive ----------------------------------------------------------------

# drive_auth()
# 1
# 
# # Summary catch data for EBS and NBS with NBS 2018 included
# googledrive::drive_download(file=as_id("1TctmzLjuFUopvdD9jqBnhNxw1NuAi87c"), 
#                             path=here::here("data","EBS_NBS_Index.RDS"), 
#                             overwrite = TRUE)
# 
# # Summary size compositions
# googledrive::drive_download(file=as_id("1EWD-0HM_WOyfbVa66o5KUVbIwDJczHKT"), 
#                             path=here::here("data","EBS_NBS_Sizecomp.RDS"), 
#                             overwrite = TRUE)
# 
# # Pacific cod unstratified age-length key
# googledrive::drive_download(file=as_id("1UEXMTmDZbUAVS6aKG8nNHfkuQVfSPpek"), 
#                             path=here::here("data","pcod_unstratified_alk_2019.RDS"), 
#                             overwrite = TRUE)
# 
# # EBS Strata
# googledrive::drive_download(file=as_id("1UZ4OBGuwqSpPhTD3idDUNqO57Fm1GYnL"), 
#                             path=here::here("data", "shapefiles","EBS_NBS_2019_1983.cpg"), 
#                             overwrite = TRUE)
# googledrive::drive_download(file=as_id("1GhH47aoQ42kx3TYqMo_w0zTV4UeXYWsN"), 
#                             path=here::here("data", "shapefiles","EBS_NBS_2019_1983.dbf"), 
#                             overwrite = TRUE)
# googledrive::drive_download(file=as_id("1SLOH6Ggp8PufL0XZPXLa8ZzPrwIgz68S"), 
#                             path=here::here("data", "shapefiles","EBS_NBS_2019_1983.sbn"), 
#                             overwrite = TRUE)
# googledrive::drive_download(file=as_id("1Hk_t3RvwqwHp6ypL4yfUsACXfxq9IynZ"), 
#                             path=here::here("data", "shapefiles","EBS_NBS_2019_1983.prj"), 
#                             overwrite = TRUE)
# googledrive::drive_download(file=as_id("16GPjJfiWL5ZNCHdimKvoQVp63PGvrVmn"), 
#                             path=here::here("data", "shapefiles","EBS_NBS_2019_1983.sbx"), 
#                             overwrite = TRUE)
# googledrive::drive_download(file=as_id("1Pfg-HxarbSU2M0u-HzSHAIj7E8v-MeO4"), 
#                             path=here::here("data", "shapefiles","EBS_NBS_2019_1983.shp"), 
#                             overwrite = TRUE)
# googledrive::drive_download(file=as_id("1Ke_9cy5wwXolzx34TBM1gbRPGVaS2g7b"), 
#                             path=here::here("data", "shapefiles","EBS_NBS_2019_1983.shp.xml"), 
#                             overwrite = TRUE)
# googledrive::drive_download(file=as_id("1wqdoTKjVSziRdQnUQa7COuYU_2H_TyAt"), 
#                             path=here::here("data", "shapefiles","EBS_NBS_2019_1983.shx"), 
#                             overwrite = TRUE)
# 
# 


# Settings ----------------------------------------------------------------

Version <- "VAST_v13_1_0"
Region <- c("Eastern_Bering_Sea","Northern_Bering_Sea")
Method = "Mesh"
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

#### Explore the data ####

Date = Sys.Date()
RunDir = paste0(workDir,"/Comps_",Date,"_",Species,"_npool=",Npool,"_BiasCorr=",BiasCorr,"/")
dir.create(RunDir)
setwd(RunDir)


# Age Composition ---------------------------------------------------------
# plus_group <- 12
# min_year <- 1994
# 
# # Load age-length keys produced by sumfish
# alk_all <- readRDS(here::here("data","pcod_unstratified_alk_2019.RDS") )
# alk_ebs <- alk_all$EBS %>%
#   dplyr::filter(SPECIES_CODE == species) %>%
#   mutate(REGION = "EBS")
# alk_nbs <- alk_all$NBS %>%
#   bind_rows( filter(alk_ebs, YEAR == 2018) ) %>%   # Use EBS ALK for 2018 ad hoc sampling in NBS
#   dplyr::filter(SPECIES_CODE == species) %>%
#   mutate(REGION = "NBS")
# 
# alk <- bind_rows(alk_ebs, alk_nbs)
# alk <- alk_all
# 
# sizeComp <- readRDS(here::here("data","EBS_NBS_SizeComp.RDS") ) %>%
#   dplyr::filter(YEAR >= min_year,
#                 SPECIES_CODE == species,
#                 !is.na(EFFORT)
#   )
# 
# haulData <- readRDS(here::here("data","EBS_NBS_Index.RDS") ) %>%
#   dplyr::filter(YEAR >= min_year,
#                 SPECIES_CODE == species,
#                 !is.na(EFFORT)
#   )
# 
# # Get summary calculations - this also fills zeroes within year - add sex if generating sex-specific agecomps
# allCats <- expand.grid(HAULJOIN=unique(haulData$HAULJOIN), AGE = unique(alk$AGE[alk$AGE<=plus_group]), noAge = 0) %>%
#   inner_join(haulData, by = c("HAULJOIN")) 
# 
# 
# # Aggregate by Age key
# Data <- sizeComp %>%
#   # left_join(alk, by = c("YEAR", "REGION", "LENGTH","SEX","SPECIES_CODE")) %>%
#   left_join(alk, by = c("YEAR", "LENGTH","SEX")) %>%
#   # mutate(ageCPUE = nSizeCPUE * probability,
#     mutate(ageCPUE = nSizeCPUE * proportion,
#          AGE = ifelse(AGE > plus_group,plus_group, AGE)) %>% 
#   group_by(YEAR,HAULJOIN,STRATUM,START_LONGITUDE, START_LATITUDE,nCPUE, AGE) %>%
#   summarize(ageCPUE = sum(ageCPUE),
#             count=n()) %>%
#   ungroup() %>%
#   select(HAULJOIN,AGE, ageCPUE, count) %>%
#   right_join(allCats, by= c("HAULJOIN","AGE")) %>%
#   mutate(ageCPUE = ifelse(is.na(ageCPUE), noAge, ageCPUE)
#   ) 
# 
# # Test CPUE
# # checkData <- Data %>%
# #   group_by(YEAR,HAULJOIN, nCPUE) %>%
# #   summarize(sum_age_cpue = sum(ageCPUE)) %>%
# #   mutate(diff = nCPUE - sum_age_cpue) 
# 
# # dbSummary <- Data %>%
# #   group_by(YEAR, STRATUM, AGE) %>%
# #   summarize(meanAgeCPUE = mean(ageCPUE)) %>%
# #   inner_join(bind_rows(NBS$stratum,EBS$stratum), by="STRATUM") %>%
# #   mutate(agePopStratum = meanAgeCPUE * STRATUM_AREA) %>%
# #   group_by(YEAR, AGE) %>%
# #   summarize(agePopTotal = sum(agePopStratum)) %>%
# #   ungroup()
# 
# # write.csv(dbSummary, "design-estimate.csv")
# 
# # Format the data for VAST
# Data_Geostat <-  transmute(Data,
#                            Catch_KG = ageCPUE,
#                            Year = YEAR,
#                            Vessel = "missing",
#                            Age = AGE,
#                            AreaSwept_km2 = .01, # Converts CPUE to km^2
#                            Lat = START_LATITUDE,
#                            Lon = START_LONGITUDE,
#                            Pass = 0
# ) %>%
#   data.frame()
# 
# write.csv(Data_Geostat, "Data_Geostat.csv")


# pollock data ------------------------------------------------------------

alk_data <- read_csv(here::here("output", "VAST_ddc_alk_2022.csv"))
table(alk_data$Age)

Data_Geostat <- alk_data %>%
  mutate(Vessel = "missing",
         AreaSwept_km2 = .01, # Converts CPUE to km^2
         Pass = 0) %>% 
  dplyr::filter(Age != 0) %>% 
  mutate(Age = Age -1) %>% 
  data.frame()

# check for sample size
table(Data_Geostat$Age)


# Run Analysis ------------------------------------------------------------

# Run model
fit = fit_model( "settings"=settings, 
                 "Lat_i"=Data_Geostat[,'Lat'], 
                 "Lon_i"=Data_Geostat[,'Lon'], 
                 "t_i"=Data_Geostat[,'Year'],  # "t_i"=rep(2019,nrow(Data_Geostat)),
                 "c_i"=Data_Geostat[,'Age'], 
                 "b_i"=Data_Geostat[,'Catch_KG'], 
                 "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                 "v_i"=Data_Geostat[,'Vessel'],
                 Npool = Npool, 
                 test_fit=T,
                 # newtonsteps = 0,       # for testing
                 # getsd = FALSE,         # for testing
                 # test_fit = FALSE,      # for testing
                 # "run_model" = FALSE,   # for testing
                 # "build_model" = FALSE, # for testing
                 # CheckForBugs = FALSE,  # for testing
                 create_strata_per_region=TRUE)

# Save results

saveRDS(fit, file = paste0(workDir,"VASTfit.RDS"))


# Plots -------------------------------------------------------------------
# If you need to load a fit in a new session:
# dyn.load(dynlib("VAST_v12_0_0"))
# load(here("VAST_results", "Walleye_Pollock_age_2022_EBS-NBS", "Comps_2022-06-28_Walleye Pollock Agecomp_npool=100_BiasCorr=TRUE","VASTresults.RDS"))
fit <- readRDS(file = paste0(workDir,"/VASTfit.RDS"))

# Record package versions
sink("session_info.txt", type = "output")
sessionInfo()
sink()

# Plot results
results <- plot_results( fit, #zrange = c(-3,3), n_cells = 600, 
                         strata_names = strata_names, 
                         check_residuals = TRUE )
saveRDS(results, file = "VASTresults.RDS")


# If residual plots don't... uh... plot...
plot_quantile_residuals( fit=fit ) 

map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
plot_maps( Obj=fit$tmb_list$Obj, PlotDF=map_list[["PlotDF"]] )



# Expand to proportional population numbers -------------------------------
VASTfit <- fit
# results = plot_results( settings=settings, fit=VASTfit, check_residuals=FALSE )

## Variable names have changed in FishStatsUtils check for consistency ?calculate_proportion
Year_Set =VASTfit$year_labels
Years2Include = which(VASTfit$year_labels != 2020)
proportions = calculate_proportion( TmbData=VASTfit$data_list, 
                                    Index=results$Index, 
                                    Year_Set= Year_Set, 
                                    Years2Include = Years2Include,
                                    strata_names=strata_names)

prop <- data.frame(t(data.frame(proportions$Prop_ctl))) %>%
    # drop_na() %>%
    rename("age_1"=1,"age_2"=2,"age_3"=3,"age_4"=4,"age_5"=5,"age_6"=6,"age_7"=7,"age_8"=8,"age_9"=9,"age_10"=10,
           "age_11"=11,"age_12"=12,"age_13"=13,"age_14"=14,"age_15+"=15) %>%
    bind_cols(data.frame(Year = rep(Year_Set,3) ), 
              data.frame(Region = c(rep(strata_names[1],length(Year_Set)),
                                    rep(strata_names[2],length(Year_Set)),
                                    rep(strata_names[3],length(Year_Set)) ) 
              )
    )

write.csv(prop,"proportions.csv")



# check results -----------------------------------------------------------

load("parameter_estimates.RData")
cbind(parameter_estimates$SD$value, parameter_estimates$SD$sd)

