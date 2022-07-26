# This code produces VAST index estimates for the Bering Sea Pollock assessment
# This code is converted from example code by Jason Connor
# Modifided by: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.04.23
# Date updated: 2021.09.24


# package installs --------------------------------------------------------

### 2021 settings (updated packages): Rv4.0.2 VAST v3.6.1, FishStatsUtils v2.8.0, cpp VAST_v12_0_0

## for more VAST install info, go to https://github.com/James-Thorson-NOAA/VAST
# Install TMB from CRAN
# install.packages("TMB")
# Install INLA using currently recommended method
# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# # Install FishStatsUtils from CRAN
# devtools::install_github("james-thorson/FishStatsUtils@2.8.0") #, INSTALL_opts="--no-staged-install")
# # Install VAST at specific version
# devtools::install_github("james-thorson/VAST@3.6.1") #for a specific VAST package version

## Load packages and set working directory
library(TMB)               
library(VAST)

## check to ensure you're runing the correct package versions
packageVersion('FishStatsUtils')
packageVersion('VAST')

## how to check your session info (see R version and all package versions)
sessionInfo()
# R version 4.0.2 (2020-06-22)

# libraries ---------------------------------------------------------------

library(here) 
library(googledrive)
library(tidyverse)
library(VAST)
library(sf)
library(scales)
# library(renv)
library(lubridate)
library(janitor)

# Set R environment for consistency ---------------------------------------

# renv::init()
# renv::snapshot()
# renv::restore()


# Set species -------------------------------------------------------------

data_choice <- (readline(prompt = "Which data set do you want to use (Enter: a = EBS+NBS, e = EBS, n = NBS): "))
a

if(data_choice  == 'a') {region_select <- "EBS_NBS" 
use_region <- c("Eastern_Bering_Sea","Northern_Bering_Sea")
use_strata_names <- c("Both","EBS","NBS")
folder <- "EBS-NBS"
} else if(data_choice  == 'e') {region_select <- "EBS"; 
use_region <- c("Eastern_Bering_Sea")
use_strata_names <- c("EBS")
folder <- "EBS_only"
}else if(data_choice  == 'n') {region_select <- "NBS"; 
use_region <- c("Northern_Bering_Sea")
use_strata_names <- c("NBS")
folder <- "NBS_only"
}else(print("Invalid selection: please select a for all EBS+NBS, e for EBS, or n for NBS"))



species <- 21740
speciesName <- paste0("Walleye_Pollock_",lubridate::year(today()),"_",region_select)
workDir <- here::here("results", speciesName)
dir.create(workDir, showWarnings = FALSE)
# workDir2 <- here::here("results", speciesName)
# dir.create(workDir2, showWarnings = FALSE)
setwd(workDir)


# Get data from Google drive ----------------------------------------------------------------
# Summary catch data for EBS and NBS with NBS 2018 included
# googledrive::drive_download(file=as_id("1DUEPWgAtCJRJxo_hg2fYV8qTKJkKOwwD"), 
#                             path=here::here("data","EBS_NBS_Index.RDS"), 
#                             overwrite = TRUE)
# # EBS Grid
# googledrive::drive_download(file=as_id("1ll3uRvZDoQcDHMZFmX1TPK92xBDLYKXm"), 
#                             path=here::here("data","EBSThorsonGrid.csv"), 
#                             overwrite = TRUE)
# 
# # NBS grid 
# googledrive::drive_download(file=as_id("1tX-THrTDBBACraP1sL-Chcq74k8YGxqd"), 
#                             path=here::here("data","NBSThorsonGrid.csv"), 
#                             overwrite = TRUE)
# 
# # Cold Pool Area covariate
# googledrive::drive_download(file=as_id("119pehim03WxFjGH9tzjUF7gZDcHJeUe8"), 
#                             path=here::here("data","cpa_areas2019.csv"), 
#                             overwrite = TRUE)
# 
# googledrive::drive_download(file=as_id("1G0VnlpPf0tS3-7GPULKbWFfPfkcsxYgS"), 
#                             path=here::here("data","cold_pool_2021.csv"), 
#                             overwrite = TRUE)

# 2021 cold pool: 2021-09-29 version
# googledrive::drive_download(file=as_id("14nJ4HFkyazif8dkEU-q_Oeihz2zJlEZL"),
#                             path=here::here("data","cold_pool_2021.csv"),
#                             overwrite = TRUE)
# 
# # EBS Strata
# # CIA: NOTE-- this is 2019 strata ONLY; not 2010 for comparison
# googledrive::drive_download(file=as_id("1Pfg-HxarbSU2M0u-HzSHAIj7E8v-MeO4"), 
#                             path=here::here("data", "shapefiles","EBS_NBS_2019_1983.shp"), 
#                             overwrite = TRUE)


# Read in data ------------------------------------------------------------
index_data_all <- read_csv(here("data", "VAST_ddc_all_2021.csv"))

# index_data_EBS <- read_csv(here("output", "VAST_ddc_EBSonly_2021.csv"))
# 
# index_data_NBS <- read_csv(here("output", "VAST_ddc_NBSonly_2021.csv"))

# VAST Settings -----------------------------------------------------------
# Version <- get_latest_version( package="VAST" )
Version <- "VAST_v13_1_0"
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
Options <-  c("Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE, "treat_nonencounter_as_zero"= TRUE) #TRUE )
Aniso <- TRUE
Npool <- 100
BiasCorr <- TRUE
getJointPrecision <- TRUE
getReportCovariance <- TRUE
fine_scale <- TRUE
max_cells <- 4000
strata.limits <- data.frame(STRATA = as.factor('All_areas'))
# strata.limits <- list('All_areas' = 1:1e5)


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
    knot_method = knot_method
)

# Format catch data -------------------------------------------------------
# * CIA: modify this section for density dependent pollock data ####
# sumAll <- read_rds(here::here("data","EBS_NBS_Index.RDS"))
# Data <- sumAll %>%
# dplyr::filter(SPECIES_CODE==species)

if(data_choice  == 'a') {Data <- index_data_all
} else if(data_choice  == 'e') {Data <- index_data_EBS
}else if(data_choice  == 'n') {Data <- index_data_NBS
}else(print("Invalid selection: please select a for all EBS+NBS, e for EBS, or n for NBS"))

Data

# check data
# # check for years with 100% encounter rate



# Format the data for VAST
Data_Geostat <-  dplyr::transmute(Data,
                                  Catch_KG = ddc_cpue_kg_ha*100, #  CPUE in kg/ha this converts it to kg/km^2
                                  Year = year,
                                  Vessel = "missing",
                                  AreaSwept_km2 = 1, # area swept is 1 when using CPUE instead of observed weight
                                  Lat = start_latitude,
                                  Lon = start_longitude,
                                  Pass = 0
)

unique(Data_Geostat$Vessel)
unique(Data_Geostat$AreaSwept_km2)
unique(Data_Geostat$Pass)

Data_Geostat <- Data_Geostat %>% #dplyr::filter(Year <2020) %>% 
    as.data.frame()

# VAST data format:
# # COlumns: "Catch_KG" "Year" "Vessel" "AreaSwept_km2" "Lat" "Lon" "Pass" 
# #     Vessel = "missing
# #     AreaSwept_km2 = 1
# #     Pass = 0


# Cold Pool covariate -----------------------------------------------------

# CIA: there is something WRONG in this section with the apply function; getting different values
# 
#     CP <- read.csv(here::here("data","cpa_areas2019.csv"))
#     covariate_data <- CP[ which(as.integer(CP[,'YEAR']) %in% unique(Data_Geostat[,'Year'])), ]
#     covariate_data <- data.frame( "Year"=covariate_data[,"YEAR"], "Lat"=mean(Data_Geostat[,'Lat']), "Lon"=mean(Data_Geostat[,'Lon']), apply(covariate_data[-1],MARGIN=2,FUN=function(vec){(vec-mean(vec))/sd(vec)}) )
#     
#     
#     covariate_data_setup <- CP %>% 
#         dplyr::filter(YEAR %in% unique(Data_Geostat$Year)) %>% 
#         rename(Year = YEAR) %>%    
#         mutate(Lat = mean(Data_Geostat$Lat),
#                Lon = mean(Data_Geostat$Lat)) 
#         
#         
#     cov_data_setup2 <- covariate_data_setup %>% 
#         dplyr::select(-Year, -Lat, -Lon) %>% 
#         map(function(.x) .x-mean(.x)/sd(.x)) %>% 
#         as_tibble()
#         
#         
#     covariate_data <- dplyr::select(covariate_data_setup, Year, Lat, Lon) %>%    
#         bind_cols(cov_data_setup2) %>% 
#         dplyr::select(Year, Lat, Lon, AREA_SUM_KM2_LTE2) %>% 
#         mutate(AREA_SUM_KM2_LTE2 = AREA_SUM_KM2_LTE2/1000000) %>%  #SCALE (more reasonable)
#         as.data.frame()



# 2021 cold pool --------------------------------------------------------------------
cpa_raw <- read_csv("F:/R/VAST2021/data/cpa_out_ste_simplified.csv")
CP <- cpa_raw %>% #read.csv(here::here("data","cold_pool_2021.csv")) %>% 
    clean_names() %>% 
    dplyr::select(-last_update)

covariate_data_setup <- CP %>% 
    dplyr::filter(year %in% unique(Data_Geostat$Year)) %>% 
    rename(Year = year) %>%    
    mutate(Lat = mean(Data_Geostat$Lat),
           Lon = mean(Data_Geostat$Lat)) 


cov_data_setup2 <- covariate_data_setup %>% 
    dplyr::select(-Year, -Lat, -Lon) %>% 
    map(function(.x) .x-mean(.x)/sd(.x)) %>% 
    as_tibble()

cov_2020 <- bind_cols(Year = 2020,
                      Lat = covariate_data_setup$Lat[1],
                      Lon = covariate_data_setup$Lon[1],
                      ste_lte2 = 0)
# ste_lte2 = 999999)

covariate_data <- dplyr::select(covariate_data_setup, Year, Lat, Lon) %>%    
    bind_cols(cov_data_setup2) %>% 
    dplyr::select(Year, Lat, Lon, area_lte2_km2) %>% 
    mutate(ste_lte2 = area_lte2_km2/1000000) %>%  #SCALE (more reasonable)
    bind_rows(cov_2020) %>%
    arrange(Year) %>% 
    as.data.frame()

# Load covariates
formula <- ~ ste_lte2
Xconfig_zcp <- array(2, dim=c(2,1,1) )
X1config_cp <- as.matrix(2)
X2config_cp <- as.matrix(2)


# read saved image (for troubleshooting) ----------------------------------

# load('pollock_index.RData')

# Build the model ---------------------------------------------------------
fit <- fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,'Lat'], 
                  "Lon_i"=Data_Geostat[,'Lon'], 
                  "t_i"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)),
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  # "v_i"=Data_Geostat[,'Vessel'],
                  "create_strata_per_region"=TRUE,
                  "getJointPrecision"=getJointPrecision, 
                  "getReportCovariance"=getReportCovariance,
                  "X1_formula"=formula,
                  "X2_formula"=formula,
                  "X1config_cp" = X1config_cp,
                  "X2config_cp" = X2config_cp ,
                  "covariate_data"=covariate_data,
                  #"Npool" = Npool,
                  # "run_model" = FALSE,
                  # "build_model" = FALSE,
                  # test_fit = FALSE,
                  getsd=TRUE,
                  # newtonsteps=0,
                  # CheckForBugs = FALSE,
                  "working_dir" = workDir,
                  # "knot_method" = "grid",
                  # "create_strata_per_region"=TRUE,
                  "getJointPrecision"=getJointPrecision,
                  "getReportCovariance"=getReportCovariance
                  # "Xconfig_zcp"=Xconfig_zcp
                  
                  # debug: run_model = FALSE; build_model = FALSE
                  # CheckForBugs = FALSE
)


# Save results

saveRDS(fit, file = "VASTfit.RDS")

# Plots -------------------------------------------------------------------
# If you need to load a fit in a new session:
dyn.load(dynlib("VAST_v12_0_0"))

# Record package versions
sink("session_info.txt", type = "output")
sessionInfo()
sink()

# Plot results
results <- plot_results( fit, 
                         zrange = c(-3,3), 
                         n_cells = 600, 
                         strata_names = strata_names, 
                         check_residuals=TRUE,
                         n_samples=0 )

saveRDS(results, file = "VASTresults.RDS")


# If residual plots don't... uh... plot...
plot_quantile_residuals( fit=fit ) 

map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
plot_maps( Obj=fit$tmb_list$Obj, PlotDF=map_list[["PlotDF"]] )


# ESP products
write.csv( results$Range$COG_Table, file="COG.csv", row.names=FALSE )

##save effective area occupied for ESP request
report = TMB::summary.sdreport(fit$parameter_estimates$SD)
ln_km2 = report[which(rownames(report)=="log_effective_area_ctl"),c('Estimate','Std. Error')]
Year <- sort(unique(fit$year_labels))
ln_km2 <- as.data.frame(cbind(ln_km2, Year))
ln_km2 <- ln_km2[which(ln_km2$Year %in% unique(fit$data_frame$t_i)),]
write.csv( ln_km2, file="ln_effective_area.csv", row.names=FALSE )



# Plot DB vs VAST Comparison ----------------------------------------------

# Load shapefile with strata boundaries
EBS_strata <- st_read( here("data","shapefiles"), layer = "EBS_NBS_2019_1983") 
total_area <- sum(EBS_strata$AREA_KM2)
EBS_strata <- mutate(EBS_strata, area_ratio = AREA_KM2/total_area)

# Convert Data_Geostat to sf and transform projection
sf_geostat <- st_as_sf(Data_Geostat, coords = c("Lon", "Lat"), remove = FALSE, 
                       crs = "+proj=longlat +datum=NAD83" ) %>%
    st_transform(crs = st_crs(EBS_strata))

# Spatially intersect the two dataframes
sf_geostat_strata <- st_join(sf_geostat, EBS_strata)
na_check <- sf_geostat_strata[is.na(sf_geostat_strata$Stratum),]
sf_geostat_strata[is.na(sf_geostat_strata$Catch_KG), "Catch_KG"] <- 0

# Calculate the DBEs
Strata_Geostat <- sf_geostat_strata %>%
    filter(!is.na(STRATUM) ) %>%
    group_by(Year, STRATUM) %>%
    summarize( mean_cpue = mean(Catch_KG),
               area = first(AREA_KM2),
               var_cpue = ( var(Catch_KG) / n() ),
               Stratum_biomass = mean_cpue * area,
               Stratum_biomass_var = ifelse(is.na(var_cpue),0, var_cpue * area^2)
    ) %>%
    ungroup()

DB_Geostat <- Strata_Geostat %>%
    group_by(Year) %>%
    summarize( total_biomass_mt = sum(Stratum_biomass)/1000,
               total_biomass_mt_se = sqrt(sum(Stratum_biomass_var))/1000
    ) %>%
    mutate(Estimator = "Design-based") %>%
    ungroup()

VASTindex <- read_csv( paste0(workDir,"/Table_for_SS3.csv") ) %>%
    mutate(Estimator = "VAST") %>%
    filter(Fleet == "Both")

index_compare <- bind_rows(
    select(DB_Geostat, Year, total_biomass_mt, total_biomass_mt_se,Estimator),
    select(VASTindex, Year, total_biomass_mt = Estimate_metric_tons, total_biomass_mt_se = SD_mt,Estimator)
)


pd <- position_dodge(.9)
plotLimit <- 1.1*max(index_compare$total_biomass_mt+2*index_compare$total_biomass_mt_se)

spPlot <- ggplot(index_compare, aes(x=Year, y=total_biomass_mt, color=Estimator, group=Estimator)) +
    geom_errorbar(aes(ymin=ifelse(total_biomass_mt-(2*total_biomass_mt_se) < 0, 0,total_biomass_mt-(2*total_biomass_mt_se)),ymax=total_biomass_mt+(2*total_biomass_mt_se)),
                  width=.9,
                  show.legend=T,
                  position=pd) +
    geom_point(aes(shape=Estimator), size=1.6, position=pd) +
    theme_bw() +
    ggtitle(paste0("Comparison of design-based and VAST estimators for ",speciesName)) +
    ylab("Index (1,000 fish) with 2 Standard Errors") +
    scale_y_continuous(expand = c(0, 0), 
                       breaks=pretty_breaks(n=5), 
                       labels=comma,
                       limits = c(0, plotLimit)
    ) + 
    scale_x_continuous(breaks=seq(1982,2020,2)) +
    scale_color_manual(breaks= c("Design-based","VAST","Database"), values= c("blue","red","purple")) +
    scale_shape_manual(breaks= c("Design-based","VAST","Database"), values= c(15,16,17)) +
    theme(axis.text.x=element_text(angle=90,vjust=.5))
print(spPlot)

ggsave(paste0(speciesName,"_Index_compare.png"), width=8, height=6, units="in")


# 
# # Plots -------------------------------------------------------------------
# 
# plot( fit, zrange = c(-3,3), n_cells = 600, strata_names = strata_names )
# 
# 
# 
# 
# # Plot DB vs VAST Comparison ----------------------------------------------
# 
# # Load shapefile with strata boundaries
# EBS_strata <- st_read( here("data", "shapefiles"), layer = "EBS_NBS_2019_1983") 
# 
# # Convert Data_Geostat to sf and transform projection
# sf_geostat <- st_as_sf(Data_Geostat, coords = c("Lon", "Lat"), remove = FALSE, 
#                        crs = "+proj=longlat +datum=NAD83" ) %>%
#     st_transform(crs = st_crs(EBS_strata))
# 
# # Spatially intersect the two dataframes
# sf_geostat_strata <- st_join(sf_geostat, EBS_strata)
# na_check <- sf_geostat_strata[is.na(sf_geostat_strata$Stratum),]
# 
# # Calculate the DBEs
# Strata_Geostat <- sf_geostat_strata %>%
#     filter(!is.na(STRATUM) & !STRATUM %in% c(70,71,81)) %>%
#     group_by(Year, STRATUM) %>%
#     summarize( mean_cpue = mean(Catch_KG),
#                area = first(AREA_KM2),
#                var_cpue = ( var(Catch_KG) / n() ),
#                Stratum_biomass = mean_cpue * area,
#                Stratum_biomass_se = ifelse(is.na(var_cpue),0, sqrt(var_cpue) * area)
#     ) %>%
#     ungroup()
# 
# DB_Geostat <- Strata_Geostat %>%
#     group_by(Year) %>%
#     summarize( total_biomass_mt = sum(Stratum_biomass)/1000,
#                total_biomass_mt_se = sum(Stratum_biomass_se)/1000
#     ) %>%
#     mutate(Estimator = "Design-based") %>%
#     ungroup()
# 
# VASTindex <- read_csv( paste0(workDir,"/Table_for_SS3.csv") ) %>%
#     mutate(Estimator = "VAST") %>%
#     filter(Fleet == "EBS")
# 
# index_compare <- bind_rows(
#     select(DB_Geostat, Year, total_biomass_mt, total_biomass_mt_se,Estimator),
#     select(VASTindex, Year, total_biomass_mt = Estimate_metric_tons, total_biomass_mt_se = SD_mt,Estimator)
# )
# 
# 
# pd <- position_dodge(.9)
# plotLimit <- 1.1*max(index_compare$total_biomass_mt+1*index_compare$total_biomass_mt_se)
# 
# spPlot <- ggplot(index_compare, aes(x=Year, y=total_biomass_mt, color=Estimator, group=Estimator)) +
#     geom_errorbar(aes(ymin=ifelse(total_biomass_mt-(1*total_biomass_mt_se) < 0, 0,total_biomass_mt-(1*total_biomass_mt_se)),ymax=total_biomass_mt+(1*total_biomass_mt_se)),
#                   width=.9,
#                   show.legend=T,
#                   position=pd) +
#     geom_point(aes(shape=Estimator), size=1.6, position=pd) +
#     theme_bw() +
#     ggtitle(paste0("EBS Comparison of design-based and VAST estimators of biomass for ",speciesName)) +
#     ylab("Biomass (mt) with Standard Error") +
#     scale_y_continuous(expand = c(0, 0), 
#                        breaks=pretty_breaks(n=5), 
#                        labels=comma,
#                        limits = c(0, plotLimit)
#     ) + 
#     scale_x_continuous(breaks=seq(1984,2019,2)) +
#     scale_color_manual(breaks= c("Design-based","VAST","Database"), values= c("blue","red","purple")) +
#     scale_shape_manual(breaks= c("Design-based","VAST","Database"), values= c(15,16,17)) +
#     theme(axis.text.x=element_text(angle=90,vjust=.5))
# print(spPlot)
# 
# ggsave(paste0(speciesName,"_Index_compare.png"), width=8, height=6, units="in")
