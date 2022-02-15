#' VAST 2021 Example Script
#' 
#' Use this script to standardize VAST runs for Mod Squad products in 2021. 

library(here) 
library(googledrive)
library(tidyverse)
library(VAST) 
library(sf)
library(scales)

# Set species -------------------------------------------------------------

species <- 10210
speciesName <- "YellowfinSole_EBS_NBS_index_weight_gamma"
species_name <- "yellowfin_sole"

workDir <- here::here("results",speciesName)
#dir.create(workDir)
#setwd(workDir)
# Set up folder to store species specific results
folder <- paste0(getwd(),"/species_specific_code/Bering/",species_name)
dir.create(folder)
folder <- paste0(getwd(),"/species_specific_code/Bering/",species_name,"/data")
dir.create(folder)
folder <- paste0(getwd(),"/species_specific_code/Bering/",species_name,"/results")
dir.create(folder)

# Load the data for VAST
Data_Geostat <- readRDS(file = paste0(getwd(),"/species_specific_code/Bering/",species_name,"/data/Data_Geostat.rds"))
Data_Geostat$Catch_KG[which(is.na(Data_Geostat$Catch_KG))] <- 0

# VAST Settings -----------------------------------------------------------
Version <- "VAST_v13_1_0"
Region <- c("Eastern_Bering_Sea","Northern_Bering_Sea")
strata_names = c("Both","EBS","NBS")
Method <- "Mesh"
knot_method <- "grid"
grid_size_km <- 25
n_x <- 750   # Specify number of stations (a.k.a. "knots")
FieldConfig <- c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig <- c("Beta1"=0, "Beta2"=0, "Epsilon1"=4, "Epsilon2"=4)
OverdispersionConfig <- c("Eta1"=0, "Eta2"=0)
ObsModel <- c(2,1)
Options <-  c("Calculate_Range"=TRUE, "Calculate_effective_area"=TRUE, "treat_nonencounter_as_zero"=FALSE )
Aniso <- TRUE
Npool <- 100
BiasCorr <- TRUE
getJointPrecision <- TRUE
getReportCovariance <- TRUE
fine_scale <- TRUE
max_cells <- 4000
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
    knot_method = knot_method,
    bias.correct = BiasCorr
    )



# Cold Pool covariate -----------------------------------------------------

covariate_data <- readRDS(file = paste0(getwd(),"/species_specific_code/Bering/",species_name,"/data/Data_ColdPool.rds"))

    
    # Load covariates
    formula <- ~ AREA_LTE2_KM2
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
                  "X2config_cp" <- X2config_cp ,
                  "covariate_data"=covariate_data,
                  "Npool" = Npool,
                  "test_fit" = TRUE,
                  #"working_dir" = workDir
                  "working_dir" = paste0(getwd(),"/species_specific_code/Bering/",species_name,"/results"),
                  
                  
                  
)


# Save results

saveRDS(fit, file = paste0(getwd(),"/species_specific_code/Bering/",species_name,"/results/VASTfit.RDS"))


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
   

