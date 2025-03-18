library(tidyverse)
library(VAST) 
library(sf)
library(scales)
#remotes::install_github("afsc-gap-products/coldpool")
library(coldpool)
library(akgfmaps)

# Set species, model -------------------------------------------------------

which_model <- c("hindcast", "production")[1]
species <- 21720
species_name <- "pacific_cod"

workDir <- paste0(getwd(),"/species_specific_code/BS/", 
                      species_name, "/", which_model, "/")
if(!dir.exists(workDir))
  dir.create(path = workDir, recursive = TRUE)
if(!dir.exists(paste0(workDir, "results_b/")))
  dir.create(path = paste0(workDir, "results_b/"), recursive = TRUE)


# Record sessionInfo -------------------------------------------------------
sink(file = paste0(workDir, "results_b/session_info.txt"), 
     type = "output")
sessionInfo()
sink()

# Make sure package versions are correct for current year ------------------
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

# VAST Settings -----------------------------------------------------------
Region <- c("Eastern_Bering_Sea","Northern_Bering_Sea")
strata_names = c("Both","EBS","NBS")
knot_method <- "grid"
grid_size_km <- 25
n_x <- 750
FieldConfig <- c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig <- c("Beta1"=0, "Beta2"=0, "Epsilon1"=4, "Epsilon2"=4)
OverdispersionConfig <- c("Eta1"=0, "Eta2"=0)
ObsModel <- c(2,4)
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
saveRDS(fit, file = paste0(workDir, "results_b/VASTfit.RDS"))

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
                         working_dir = paste0(workDir, "results_b/" )
                        )

saveRDS(results, file = paste0(workDir, "results_b/VASTresults.RDS"))

# residuals
plot_quantile_residuals( fit=fit )
map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
plot_maps(fit=fit, n_cells = 600, Obj=fit$tmb_list$Obj, PlotDF=map_list[["PlotDF"]] )

# ESP products
write.csv( results$Range$COG_Table, file=paste0(workDir,"results_b/COG.csv"), row.names=FALSE )

##save effective area occupied for ESP request
report = TMB::summary.sdreport(fit$parameter_estimates$SD)
ln_km2 = report[which(rownames(report)=="log_effective_area_ctl"),c('Estimate','Std. Error')]
Year <- sort(unique(fit$year_labels))
ln_km2 <- as.data.frame(cbind(ln_km2, Year))
ln_km2 <- ln_km2[which(ln_km2$Year %in% unique(fit$data_frame$t_i)),]
write.csv( ln_km2, file=paste0(workDir,"results_b/ln_effective_area.csv"), row.names=FALSE )


# Plot design-based estimate vs VAST estimate ------------------------------

# query db index
sql_channel <- gapindex::get_connected() # enter credentials in pop-out window

db <- RODBC::sqlQuery(channel = sql_channel, 
                      query = 
                        "
WITH FILTERED_STRATA AS (
SELECT AREA_ID, DESCRIPTION FROM GAP_PRODUCTS.AKFIN_AREA
WHERE AREA_TYPE = 'REGION'
AND SURVEY_DEFINITION_ID IN (143, 98))

-- Select columns for output data
SELECT 
BIOMASS_MT,
BIOMASS_VAR,
YEAR, 
DESCRIPTION

-- Identify what tables to pull data from
FROM GAP_PRODUCTS.AKFIN_BIOMASS BIOMASS
JOIN FILTERED_STRATA STRATA 
ON STRATA.AREA_ID = BIOMASS.AREA_ID

-- Filter data results
WHERE BIOMASS.SPECIES_CODE = 21720")
    
db <- db %>% 
  mutate(stratum = recode(DESCRIPTION, 
                          `EBS Standard Plus NW Region: All Strata` = "EBS", 
                          `NBS Region: All Strata` = "NBS"),
         se_estimate = sqrt(BIOMASS_VAR) * 1000,
         estimate_kg = BIOMASS_MT * 1000)  %>% 
  filter(stratum != "EBS Standard Region: All Strata") %>%
  select(stratum,
         year = YEAR, 
         estimate_kg,
         se_estimate) %>%
  mutate(estimator = "db")

# load model-based index
mb <- read.csv(paste0(workDir, "results_b/Index.csv")) 
mb <- mb %>% filter(Stratum != "Both") %>%
  select(stratum = Stratum,
         year = Time,
         estimate_kg = Estimate,
         se_estimate = Std..Error.for.Estimate) %>%
  mutate(estimator = "mb")

d <- bind_rows(db, mb) %>% filter(year != 2020)
saveRDS(d, file=paste0(workDir,"results_b/db_mb_index.RDS"))

library(ggplot2)
ggplot(d, aes(year, estimate_kg, group = interaction(stratum, estimator), 
              color = estimator, shape = stratum)) +
  geom_pointrange(aes(ymin = estimate_kg - se_estimate, 
                      ymax = estimate_kg + se_estimate),
                  size = 0.3,
                  position = position_dodge(width = 0.5)) + 
  theme_bw()
ggsave(file=paste0(workDir,"results_b/db_mb_index_comparison.png"), height = 4, width = 6, units = c("in"))
