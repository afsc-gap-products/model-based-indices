##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Model-based estimate of abundance index for EBS-NBS YFS
## Survey POC:    Zack Oyafuso (zack.oyafuso@noaa.gov)
##                Emily Markowitz (emily.markowitz@noaa.gov)
## Output POC:    Ingrid Spies (ingrid.spies@noaa.gov)
## Description:   2023 hindcast index standardization for yellowfin sole
##                     (Limanda aspera) in the Eastern and Northern Bering Sea
##  NOTES: 
##                Make sure R and package versions are consistent with those 
##                versions speficied in the 2023 Terms of Reference (TOR)
##                https://docs.google.com/document/d/18CeXcHhHK48hrtkiC6zygXlHj6YVrWEd/edit
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
which_model <- c("hindcast", "production")[1]
species_name <- "yellowfin_sole"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(coldpool) 
library(VAST) 
library(gapindex)

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
if(!dir.exists(paste0(folder, "results_index/"))) 
  dir.create(path = paste0(folder, "results_index/"), recursive = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Record sessionInfo ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink(file = paste0(folder, "results_index/session_info.txt"), 
     type = "output")
sessionInfo()
sink()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load VAST data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
Data_Geostat <- 
  readRDS(file = paste0(folder, "data/data_geostat_biomass_index.rds"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load coldpool covariate data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mean_bt_lt100m <- scale(coldpool:::cold_pool_index$MEAN_BT_LT100M)
# covariate_data <- data.frame(Year = c(coldpool:::cold_pool_index$YEAR, 2020),
#                              Lat = mean(Data_Geostat$Lat),
#                              Lon = mean(Data_Geostat$Lon),
#                              mean_bt_lt100m = c(mean_bt_lt100m, 0))

# coldpool package wouldnt load on this version of R
# write.csv(file = paste0(folder, "data/mean_bt_lt100m.csv"),
#           x = data.frame(year = coldpool:::cold_pool_index$YEAR,
#                          mean_bt_lt100m = scale(coldpool:::cold_pool_index$MEAN_BT_LT100M)))
mean_bt_lt100m <- read.csv(file = paste0(folder, "/data/mean_bt_lt100m.csv"))
covariate_data <- data.frame(Year = c(mean_bt_lt100m$year, 2020),
                             Lat = mean(Data_Geostat$Lat),
                             Lon = mean(Data_Geostat$Lon), 
                             mean_bt_lt100m = c(mean_bt_lt100m$mean_bt_lt100m, 0))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set VAST settings ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
settings <- FishStatsUtils::make_settings( 
  n_x = 750,
  Region = c("Eastern_Bering_Sea", "Northern_Bering_Sea"),
  purpose = "index2",
  fine_scale = TRUE,
  strata.limits = data.frame(STRATA = as.factor('All_areas')),
  ObsModel = c(2, 1),
  FieldConfig = c("Omega1" = "IID", "Epsilon1" = "IID", 
                  "Omega2" = "IID", "Epsilon2" = "IID"),
  RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 4, "Epsilon2" = 4),
  OverdispersionConfig = c("Eta1" = 0, "Eta2" = 0),
  Options = c("Calculate_Range" = TRUE, 
              "Calculate_effective_area" = TRUE, 
              "treat_nonencounter_as_zero" = FALSE),
  use_anisotropy = TRUE,
  Version = VAST_cpp_version,
  max_cells = 2000,
  knot_method = "grid",
  bias.correct = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Run VAST model ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit <- FishStatsUtils::fit_model( 
  
  ## Input settings
  "working_dir" = paste0(folder, "results_index/"), 
  "settings" = settings, 
  "create_strata_per_region" = TRUE,
  
  ## Input Data
  "Lat_i" = Data_Geostat[,'Lat'], 
  "Lon_i" = Data_Geostat[,'Lon'], 
  "t_i" = Data_Geostat[,'Year'], 
  "c_i" = rep(0, nrow(Data_Geostat)), 
  "b_i" = Data_Geostat[,'Catch_KG'], 
  "a_i" = Data_Geostat[,'AreaSwept_km2'], 
  "v_i" = Data_Geostat[,'Vessel'],
  
  ## Output settings
  "getJointPrecision" = TRUE, 
  "getReportCovariance" = TRUE,
  "getsd" = TRUE,
  
  ## Model tuning
  "newtonsteps" = 2,
  "Npool" = 100,
  ## Covariate data
  "X1_formula"= ~ mean_bt_lt100m,
  "X2_formula"= ~ mean_bt_lt100m,
  "X1config_cp" <- as.matrix(2),
  "X2config_cp" <- as.matrix(2) ,
  "covariate_data" = covariate_data)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save results locally ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Save VAST fit object
saveRDS(object = fit, 
        file = paste0(folder, "results_index/", species_name, "_VASTfit.RDS"))

## Save diagnostics and other outputs
results <- FishStatsUtils::plot_results( 
  fit = fit, 
  working_dir = paste0(folder, "results_index/output_plots/"),
  zrange = c(-3, 3), 
  n_cells = 600, 
  strata_names = c("Both", "EBS", "NBS"), 
  check_residuals = TRUE,
  n_samples = 0)

saveRDS(object = results, 
        file = paste0(folder, "results_index/output_plots/diagnostics.RDS"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Compare to Design-Based Estiamtes ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
chl <- gapindex::get_connected()

## Extract DBEs from GAP_PRODUCTS
yfs_dbe <- RODBC::sqlQuery(channel = chl,
                           query = "SELECT 
                                    YEAR, 
                                    CASE
                                     WHEN AREA_ID = 99900 THEN 'EBS'
                                     WHEN AREA_ID = 99902 THEN 'NBS'
                                    END AS REGION, 
                                    BIOMASS_MT,
                                    SQRT(BIOMASS_VAR)/BIOMASS_MT AS BIOMASS_CV
                                    FROM GAP_PRODUCTS.BIOMASS
                                    WHERE SPECIES_CODE = 10210
                                    AND AREA_ID in (99900, 99902)")

## Extract VAST Index Estimates
mbe <- 
  subset(x = read.csv(file = paste0(folder, 
                                    "results_index/output_plots/Index.csv")),
         subset = Stratum %in% c("EBS", "NBS"),
         select = c(Time, Stratum, Estimate, Std..Error.for.ln.Estimate.))
mbe$Estimate <- mbe$Estimate / 1000
names(x = mbe) <- c("YEAR", "REGION", "BIOMASS_MT", "BIOMASS_CV")

## Merge the two sources using YEAR and REGION as a composite key.
dbe_vast_comp <- merge(x = mbe, all.x = TRUE,
                       y = yfs_dbe, 
                       by = c("REGION", "YEAR"),
                       suffixes = c("_VAST", "_DBE"))

## Save comparison
write.csv(x = dbe_vast_comp,
          file = paste0(folder, 
                        "results_index/output_plots/dbe_vs_vast_comp.csv"),
          row.names = F)
