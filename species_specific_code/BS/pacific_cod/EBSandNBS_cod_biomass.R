## TODO: set up results directories properly, some of the output is going to root dir

library(tidyverse)
library(VAST) 
library(sf)
library(scales)
library(coldpool)
library(akgfmaps)

# Set species, model -------------------------------------------------------

which_model <- c("hindcast", "production")[2]
species <- 21720
species_name <- "pacific_cod"

workDir <- paste0(getwd(),"/species_specific_code/BS/", 
                      species_name, "/", which_model, "/")
if(!dir.exists(workDir))
  dir.create(path = workDir, recursive = TRUE)
if(!dir.exists(paste0(workDir, "results/")))
  dir.create(path = paste0(workDir, "results/"), recursive = TRUE)


# Record sessionInfo -------------------------------------------------------
sink(file = paste0(workDir, "results/session_info.txt"), 
     type = "output")
sessionInfo()
sink()

# Make sure package versions are correct for current year ------------------
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

# Format catch data -------------------------------------------------------

dat <- read_rds(paste0(workDir, "data/EBS_NBS_Index.RDS"))

# Format the data for VAST
Data_Geostat <-  dplyr::transmute(dat,
                               Catch_KG = wCPUE*100, # sumfish calculates CPUE in kg/ha this converts it to kg/km^2
                               Year = YEAR,
                               Vessel = "missing",
                               AreaSwept_km2 = 1, # area swept is 1 when using CPUE instead of observed weight
                               Lat = START_LATITUDE,
                               Lon = START_LONGITUDE,
                               Pass = 0
) %>%
    data.frame()

saveRDS(Data_Geostat, file = paste0(workDir,"data/Data_Geostat_",species_name,".RDS"))


# Cold Pool covariate -----------------------------------------------------
cpe <- scale(coldpool:::cold_pool_index$AREA_LTE2_KM2)
covariate_data <- data.frame(Year = c(coldpool:::cold_pool_index$YEAR, 2020),
                             Lat = mean(Data_Geostat$Lat),
                             Lon = mean(Data_Geostat$Lon), 
                             cpe = c(cpe, 0))

    # Load covariates
    formula <- ~ cpe
    Xconfig_zcp <- array(2, dim=c(2,1,1) )
    X1config_cp <- as.matrix(2)
    X2config_cp <- as.matrix(2)

    
# Build the model ---------------------------------------------------------
fit <- fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,'Lat'], 
                  "Lon_i"=Data_Geostat[,'Lon'], 
                  "t_i"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'],
                  "create_strata_per_region"=TRUE,
                  "getJointPrecision"=getJointPrecision, 
                  "getReportCovariance"=getReportCovariance,
                  "X1_formula"=formula,
                  "X2_formula"=formula,
                  "X1config_cp" <- X1config_cp,
                  "X2config_cp" <- X2config_cp,
                  "covariate_data"=covariate_data,
                  "test_fit" = TRUE,
                  "working_dir" = workDir
)


# Save results

saveRDS(fit, file = "VASTfit.RDS")


# Plots -------------------------------------------------------------------
    # If you need to load a fit in a new session:
    #dyn.load(dynlib(VAST_cpp_version))
    
    # Plot results
    results <- plot_results( fit, 
                             zrange = c(-3,3), 
                             n_cells = 600, 
                             strata_names = strata_names, 
                             check_residuals=TRUE,
                             n_samples=0 )
    
    saveRDS(results, file = paste0(workDir, "results/VASTresults.RDS"))
    
    
    # If residual plots don't... uh... plot...
    # plot_quantile_residuals( fit=fit )
    # 
    # map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
    # plot_maps(fit=fit, n_cells = 600, Obj=fit$tmb_list$Obj, PlotDF=map_list[["PlotDF"]] )
    
    
    # ESP products
    write.csv( results$Range$COG_Table, file=paste0(workDir,"results/COG.csv"), row.names=FALSE )
    
    ##save effective area occupied for ESP request
    report = TMB::summary.sdreport(fit$parameter_estimates$SD)
    ln_km2 = report[which(rownames(report)=="log_effective_area_ctl"),c('Estimate','Std. Error')]
    Year <- sort(unique(fit$year_labels))
    ln_km2 <- as.data.frame(cbind(ln_km2, Year))
    ln_km2 <- ln_km2[which(ln_km2$Year %in% unique(fit$data_frame$t_i)),]
    write.csv( ln_km2, file=paste0(workDir,"results/ln_effective_area.csv"), row.names=FALSE )
    


# Plot DB vs VAST Comparison ----------------------------------------------
    
    # Load shapefile with strata boundaries
    EBS_strata <- akgfmaps::get_base_layers("ebs")
    total_area <- sum(EBS_strata$survey.strata$F_AREA)
    EBS_strata <- mutate(EBS_strata$survey.strata, area_ratio = F_AREA/total_area)
    
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
        filter(!is.na(Stratum) ) %>%
        group_by(Year, Stratum) %>%
        summarize( mean_cpue = mean(Catch_KG),
                   area = first(F_AREA),
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
    
    VASTindex <- read_csv( paste0(workDir,"results/Index.csv") ) %>%
        mutate(Estimator = "VAST") %>%
        filter(Stratum == "Both")
    
    index_compare <- bind_rows(
        dplyr::select(as.data.frame(DB_Geostat), Year, total_biomass_mt, total_biomass_mt_se, Estimator),
        dplyr::mutate(VASTindex, Year=Time, total_biomass_mt = Estimate/1000, total_biomass_mt_se = `Std. Error for Estimate`/1000,Estimator)
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
        ggtitle(paste0("Comparison of design-based and VAST estimators for ",species_name)) +
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

    ggsave(paste0(workDir,"results/Index_compare.png"), width=8, height=6, units="in")
   

