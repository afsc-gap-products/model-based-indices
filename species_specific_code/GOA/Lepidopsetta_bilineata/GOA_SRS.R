### Cecilia O'Leary, cecilia.oleary@noaa.gov
### code to produce model-based indices using VAST software implemented in R
### created March 2020
### 2021 settings: Rv4.0.2 VAST v3.6.1, FishStatsUtils v2.8.0, cpp VAST_v12_0_0

#devtools::install_local("C:/Users/cecilia.oleary/Downloads/FishStatsUtils-2.8.0/FishStatsUtils-2.8.0") 
#devtools::install_local("C:/Users/cecilia.oleary/Downloads/VAST-3.6.1/VAST-3.6.1") 

# Load packages and set working directory
library(TMB)               
library(VAST)

packageVersion('VAST')
packageVersion('FishStatsUtils')
packageVersion('Matrix')
packageVersion('TMB')

# Set species 
Species <- c("Gadus_macrocephalus","Sebastes_variabilis","Sebastes_polyspinis", "Sebastes_alutus","Gadus_chalcogrammus",
             "Lepidopsetta_polyxystra","Lepidopsetta_bilineata","Hippoglossoides_elassodon","Atheresthes_stomias")[7] ##change number 1 to select a difference species from the vector
## c(P. cod, dusky rockfish, northern rockfish, POP, pollock, 
## northern rock sole, southern rock sole, flathead sole, arrowtooth)
#c(P. cod, dusky rockfish, northern rockfish, POP,pollock, arrowtooth)
#dplyr::filter(COMMON_NAME %in% c("dusky and dark rockfishes unid.", "dusky rockfish"))

# Load the data for VAST
Data_Geostat <- readRDS(file = paste0(getwd(),"/species_specific_code/GOA/",species_name,"/data/Data_Geostat_",species_name,".rds"))
Data_Geostat$Catch_KG[which(is.na(Data_Geostat$Catch_KG))] <- 0
Data_Geostat <- Data_Geostat[which(Data_Geostat$Year >= 1996),]


# Define strata
strata.limits <- data.frame(STRATA = as.factor('All_areas'))
#strata.limits <- data.frame('STRATA' = "west_of_140W", 'west_border' = -Inf, 'east_border' = -140 )
FieldConfig = matrix( c("IID","IID","IID","IID","IID","IID"), ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
RhoConfig  = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0)

# Make settings 
settings = make_settings( Version = "VAST_v13_1_0",
                          n_x = 750,#1000, 
                          Region = "User", #"gulf_of_alaska",
                          purpose = "index2", 
                          fine_scale = TRUE, 
                          ObsModel= c(2,1), #c(1,1) #c(10,2)
                          strata.limits=strata.limits, 
                          knot_method = "grid", 
                          bias.correct = TRUE,
                          use_anisotropy = TRUE)

# Import extrapolation grid, these will be available on Jason's Google drive: VASTGAP\Extrapolation Grids
GOAgrid <- read.csv(file= paste0(getwd(),"/extrapolation_grids/GOAThorsonGrid_Less700m.csv"))

input_grid=cbind(Lat=GOAgrid$Lat,Lon=GOAgrid$Lon,Area_km2=GOAgrid$Shape_Area/1000000)  # Extrapolation grid area is in m^2 and is converted to km^2
gc()

# Run model
fit = fit_model( "settings"=settings, 
                 "Lat_i"=Data_Geostat[,'Lat'], 
                 "Lon_i"=Data_Geostat[,'Lon'], 
                 "t_i"=Data_Geostat[,'Year'], 
                 "b_i"=Data_Geostat[,'Catch_KG'], 
                 "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                 "v_i"=Data_Geostat[,'Vessel'], #### ##was ok to leave in because it's all "missing" or zero, so no vessel effects
                 "refine" = TRUE,
                 "input_grid"=input_grid, 
                 optimize_args=list("lower"=-Inf,"upper"=Inf),
                 "working_dir" = paste0(paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results"))

# Plot results
plot( fit )

# save the VAST model
saveRDS(fit,file = paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results/",species_name,"VASTfit.RDS"))


## Save COG (center of gravity) for ESP request
results = plot( fit, n_cells=200^2 )
write.csv( results$Range$COG_Table, file=paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results/",species_name,"COG.csv"), row.names=FALSE )

##save effective area occupied for ESP request
report = TMB::summary.sdreport(fit$parameter_estimates$SD)
ln_km2 = report[which(rownames(report)=="log_effective_area_ctl"),c('Estimate','Std. Error')]
Year <- sort(unique(fit$year_labels))
ln_km2 <- as.data.frame(cbind(ln_km2, Year))
ln_km2 <- ln_km2[which(ln_km2$Year %in% unique(fit$data_frame$t_i)),]
write.csv( ln_km2, file=paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results/",species_name,"ln_effective_area.csv"), row.names=FALSE )