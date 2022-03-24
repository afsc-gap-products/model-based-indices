#' ---
#' Project:       Hindcast Index Calculation, EBS-NBS YFS
#' Authors:       Zack Oyafuso (zack.oyafuso AT noaa.gov)
#'                Emily Markowitz (emily.markowitz AT noaa.gov)
#'                Jason Conner (jason.conner AT noaa.gov)
#' Output POC:    Ingrid Spies (ingrid.spies AT noaa.gov)
#' Description:   2022 hindcast index standardization for yellowfin sole
#'                     (Limanda aspera) in the Eastern and Northern Bering Sea
#' NOTES: 
#'                Make sure R and package versions are consistent with those versions
#'                speficied in the 2022 Terms of Reference (TOR)
#'                https://docs.google.com/document/d/1t-pIruLZ-F_iNzCysWLH8cdsM0gZpuUb/edit
#' ---
 
rm(list = ls())
finalanalysis <- TRUE # this will make the model work faster while troubleshooting
  
# Import packages ---------------------------------------------------------

library(googledrive)
library(tidyverse)
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


# System preferences ---------------------------------------------------------

R_version <- "R version 4.0.2 (2020-06-22)"
VAST_cpp_version <- "VAST_v13_1_0"
pck_version <- c("VAST" = "3.8.2", 
                 "FishStatsUtils" = "2.10.2", 
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

# Authorize googledrive to view and manage your Drive files. ------------------

googledrive::drive_deauth()
googledrive::drive_auth() 
1

google_drive_dir <- "1JCUP8VL-22BVXT8owDqWcG6WH9Qn5CZm"

# Set species ------------------------------------------------------------------

species <- 10210
species_name <- "yellowfin_sole"

# Set up folder to store species specific results ----------------------------

folder <- paste0("species_specific_code/BS/", species_name, "/hindcast/")
if(!dir.exists(folder)) dir.create(folder)
folder <- paste0("species_specific_code/BS/", species_name,"/hindcast/results/")
if(!dir.exists(folder)) dir.create(folder)
rm(folder)

# Load the data for VAST -------------------------------------------------------

Data_Geostat <- readRDS(file = paste0("species_specific_code/BS/", 
                                      species_name, 
                                      "/hindcast/data/Data_Geostat.rds"))
Data_Geostat$Catch_KG[which(is.na(Data_Geostat$Catch_KG))] <- 0

# Load Cold Pool Covariate Data ------------------------------------------------

covariate_data <- readRDS(file = paste0("species_specific_code/BS/",
                                        species_name, 
                                        "/hindcast/data/Data_ColdPool.rds"))

# Set VAST settings ------------------------------------------------------------

settings <- FishStatsUtils::make_settings( 
  n_x = 750,
  Region = c("Eastern_Bering_Sea", "Northern_Bering_Sea"),
  purpose = "index2",
  fine_scale = TRUE,
  strata.limits = data.frame(STRATA = as.factor('All_areas')),
  ObsModel = c(2,1),
  FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", 
                  "Omega2"="IID", "Epsilon2"="IID"),
  RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 4, "Epsilon2" = 4),
  OverdispersionConfig = c("Eta1"=0, "Eta2"=0),
  Options = c("Calculate_Range" = TRUE, 
              "Calculate_effective_area" = TRUE, 
              "treat_nonencounter_as_zero" = FALSE ),
  use_anisotropy = TRUE,
  Version = VAST_cpp_version,
  max_cells = 2000,
  knot_method = "grid",
  bias.correct = TRUE
)

# Run VAST model ---------------------------------------------------------------

fit <- fit_model( "settings" = settings, 
                  "Lat_i" = Data_Geostat[,'Lat'], 
                  "Lon_i" = Data_Geostat[,'Lon'], 
                  "t_i" = Data_Geostat[,'Year'], 
                  "c_i" = rep(0, nrow(Data_Geostat)), 
                  "b_i" = Data_Geostat[,'Catch_KG'], 
                  "a_i" = Data_Geostat[,'AreaSwept_km2'], 
                  "v_i" = Data_Geostat[,'Vessel'],
                  "create_strata_per_region" = TRUE,
                  "getJointPrecision" = TRUE, 
                  "getReportCovariance" = TRUE,
                  "X1_formula"= ~ AREA_LTE2_KM2,
                  "X2_formula"= ~ AREA_LTE2_KM2,
                  "X1config_cp" <- as.matrix(2),
                  "X2config_cp" <- as.matrix(2) ,
                  "covariate_data" = covariate_data,
                  "Npool" = 100,
                  "test_fit" = TRUE,
                  "working_dir" = paste0(getwd(),"/species_specific_code/BS/",
                                         species_name, "/hindcast/results/"), 
                  getsd = ifelse(finalanalysis, TRUE, FALSE), # default
                  newtonsteps = ifelse(finalanalysis, 1, 0)#,  # default
                  # build_model = ifelse(finalanalysis, TRUE, FALSE) # default
)


# Save results locally: --------------------------------------------------------

# VAST fit: fit --> VASTfit.RDS
# session info: sessionInfo() --> session_info.txt
# diagnostic plots: written to diagnostics_plots/
# center of gravity: results$Range$COG_Table --> COG.csv
# effective area ln_km2 --> occupied: ln_effective_area.csv

## Save VAST fit object
saveRDS(object = fit, 
        file = paste0("species_specific_code/BS/",
                      species_name, "/hindcast/results/VASTfit.RDS"))

## Save package versions
sink(paste0("species_specific_code/BS/",
            species_name, "/hindcast/results/session_info.txt"), 
     type = "output")
sessionInfo()
sink()

# Save diagnostics and other outputs
results <- FishStatsUtils::plot_results( 
  fit = fit, 
  working_dir = paste0(getwd(),  "/species_specific_code/BS/",
                       species_name, "/hindcast/results/output_plots/"),
  zrange = c(-3,3), 
  n_cells = 600, 
  strata_names = c("Both", "EBS", "NBS"), 
  check_residuals=TRUE,
  n_samples = 0 )

saveRDS(object = results, 
        file = paste0("species_specific_code/BS/", species_name, 
                      "/hindcast/results/output_plots/diagnostics.RDS"))

## ESP products
write.csv( x = results$Range$COG_Table, 
           file = paste0("species_specific_code/BS/", species_name, 
                         "/hindcast/results/output_plots/COG.csv"), 
           row.names=FALSE )

## Effective area occupied for ESP request
report = TMB::summary.sdreport(fit$parameter_estimates$SD)
ln_km2 = report[which(rownames(report)=="log_effective_area_ctl"),
                c('Estimate','Std. Error')]
Year <- sort(unique(fit$year_labels))
Area <- rep(x = c("Both", "EBS", "NBS"), each = length(Year))
ln_km2 <- as.data.frame(cbind(Area, Year, ln_km2))
ln_km2 <- ln_km2[which(ln_km2$Year %in% unique(fit$data_frame$t_i)),]
write.csv(x = ln_km2, 
          file=paste0("species_specific_code/BS/", species_name, 
                      "/hindcast/results/output_plots/ln_effective_area.csv"),
          row.names=FALSE )

# Save local results to google drive -------------------------------------------

for (ifile in dir(paste0("species_specific_code/BS/", species_name, 
                         "/hindcast/results/"), 
                  full.names = T, 
                  recursive = T) ) {
  googledrive::drive_upload(media = ifile,
                            path = as_id(google_drive_dir), 
                            overwrite = TRUE)
}

# map_list = FishStatsUtils::make_map_info(
#   "Region" = settings$Region, 
#   "spatial_list" = fit$spatial_list, 
#   "Extrapolation_List" = fit$extrapolation_list)
# # Plot DB vs VAST Comparison ----------------------------------------------
# 
# # Load shapefile with strata boundaries
# EBS_strata <- st_read( here("data","shapefiles"), layer = "EBS_NBS_2019_1983")
# total_area <- sum(EBS_strata$AREA_KM2)
# EBS_strata <- mutate(EBS_strata, area_ratio = AREA_KM2/total_area)
# 
# # Convert Data_Geostat to sf and transform projection
# sf_geostat <- st_as_sf(Data_Geostat, coords = c("Lon", "Lat"), remove = FALSE,
#                        crs = "+proj=longlat +datum=NAD83" ) %>%
#   st_transform(crs = st_crs(EBS_strata))
# 
# # Spatially intersect the two dataframes
# sf_geostat_strata <- st_join(sf_geostat, EBS_strata)
# na_check <- sf_geostat_strata[is.na(sf_geostat_strata$Stratum),]
# sf_geostat_strata[is.na(sf_geostat_strata$Catch_KG), "Catch_KG"] <- 0
# 
# # Calculate the DBEs
# Strata_Geostat <- sf_geostat_strata %>%
#   filter(!is.na(STRATUM) ) %>%
#   group_by(Year, STRATUM) %>%
#   summarize( mean_cpue = mean(Catch_KG),
#              area = first(AREA_KM2),
#              var_cpue = ( var(Catch_KG) / n() ),
#              Stratum_biomass = mean_cpue * area,
#              Stratum_biomass_var = ifelse(is.na(var_cpue),0, var_cpue * area^2)
#   ) %>%
#   ungroup()
# 
# DB_Geostat <- Strata_Geostat %>%
#   group_by(Year) %>%
#   summarize( total_biomass_mt = sum(Stratum_biomass)/1000,
#              total_biomass_mt_se = sqrt(sum(Stratum_biomass_var))/1000
#   ) %>%
#   mutate(Estimator = "Design-based") %>%
#   ungroup()
# 
# VASTindex <- read_csv( paste0(workDir,"/Table_for_SS3.csv") ) %>%
#   mutate(Estimator = "VAST") %>%
#   filter(Fleet == "Both")
# 
# index_compare <- bind_rows(
#   select(DB_Geostat, Year, total_biomass_mt, total_biomass_mt_se,Estimator),
#   select(VASTindex, Year, total_biomass_mt = Estimate_metric_tons, total_biomass_mt_se = SD_mt,Estimator)
# )
# 
# 
# pd <- position_dodge(.9)
# plotLimit <- 1.1*max(index_compare$total_biomass_mt+2*index_compare$total_biomass_mt_se)
# 
# spPlot <- ggplot(index_compare, aes(x=Year, y=total_biomass_mt, color=Estimator, group=Estimator)) +
#   geom_errorbar(aes(ymin=ifelse(total_biomass_mt-(2*total_biomass_mt_se) < 0, 0,total_biomass_mt-(2*total_biomass_mt_se)),ymax=total_biomass_mt+(2*total_biomass_mt_se)),
#                 width=.9,
#                 show.legend=T,
#                 position=pd) +
#   geom_point(aes(shape=Estimator), size=1.6, position=pd) +
#   theme_bw() +
#   ggtitle(paste0("Comparison of design-based and VAST estimators for ",speciesName)) +
#   ylab("Index (1,000 fish) with 2 Standard Errors") +
#   scale_y_continuous(expand = c(0, 0), 
#                      breaks=pretty_breaks(n=5), 
#                      labels=comma,
#                      limits = c(0, plotLimit)
#   ) + 
#   scale_x_continuous(breaks=seq(1982,2020,2)) +
#   scale_color_manual(breaks= c("Design-based","VAST","Database"), values= c("blue","red","purple")) +
#   scale_shape_manual(breaks= c("Design-based","VAST","Database"), values= c(15,16,17)) +
#   theme(axis.text.x=element_text(angle=90,vjust=.5))
# print(spPlot)
# 
# ggsave(paste0(speciesName,"_Index_compare.png"), width=8, height=6, units="in")
