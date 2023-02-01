### Cecilia O'Leary, cecilia.oleary@noaa.gov
### code to produce model-based indices using VAST software implemented in R
### created March 2020
### 2020 settings: VAST v3.3.0, FishStatsUtils v2.5.0, cpp VAST_v8_2_0
### 2021 settings: Rv4.0.2 VAST v3.6.1, FishStatsUtils v2.8.0, cpp VAST_v12_0_0
### 2023 settings: Rv4.0.2 or later, VAST v3.10.0, FishStatsUtils v2.12.0, cpp VAST_v14_0_1, Matrix 1.5-3; TMB 1.9.2; DHARMa 0.4.6

# Load packages and set working directory
library(TMB)               
library(VAST)

packageVersion('VAST')
packageVersion('FishStatsUtils')
packageVersion('Matrix')
packageVersion('TMB')
packageVersion('DHARMa')

# Set species 
species_name <- "Sebastes_variabilis"
## dusky rockfish
#dplyr::filter(COMMON_NAME %in% c("dusky and dark rockfishes unid.", "dusky rockfish"))
knots <- c(500,750,100)[2]
model <- c("poisson_delta","delta")[1]
obs <- c("gamma","lognormal")[2]

# Load the data for VAST
Data_Geostat <- readRDS(file = paste0(getwd(),"/species_specific_code/GOA/",species_name,"/Data_Geostat_",species_name,".rds"))
#Data_Geostat <- readRDS(file = paste0(getwd(),"/data/Data_Geostat_",species_name,".rds"))
Data_Geostat$Catch_KG[which(is.na(Data_Geostat$Catch_KG))] <- 0

# Define strata
#strata.limits <- data.frame(STRATA = as.factor('All_areas'))
#strata c(3, 4, 1, 2)
strata.limits <- data.frame('STRATA' = as.factor(c("Total","Western","Central","Eastern")), 
                  'west_border' = c(-Inf,-Inf, -159,-147), 
                  'east_border' = c(Inf,-159,-147,Inf ))
FieldConfig = matrix( c("IID","IID","IID","IID","IID","IID"), ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) )
RhoConfig  = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0)

# Make settings 
settings = make_settings( Version = "VAST_v14_0_1",
                          n_x = knots,#1000, 
                          Region = "User", #"gulf_of_alaska",
                          purpose = "index2", 
                          fine_scale = TRUE, 
                          ObsModel= c(1,1), #c(2,1) #c(10,2)
                          strata.limits=strata.limits, 
                          knot_method = "grid", 
                          bias.correct = TRUE,
                          use_anisotropy = TRUE#,
                          #max_cells = 2000
                          )

# Import extrapolation grid, these will be available on Jason's Google drive: VASTGAP\Extrapolation Grids
GOAgrid <- read.csv(file= paste0(getwd(),"/GOA_extrapolation_grids/GOAThorsonGrid_Less700m.csv"))

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
                 "input_grid"=input_grid, 
                 "refine" = TRUE,
                 optimize_args=list("lower"=-Inf,"upper"=Inf),
                 "working_dir" = paste0(paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results")) )


# Plot results
plot( fit, "working_dir" = paste0(paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results") ) )

# save the VAST model
saveRDS(fit,file = paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results/",species_name,"VASTfit.RDS"))


## Save COG (center of gravity) for ESP request
results = plot( fit, n_cells=200^2 , "working_dir" = paste0(paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results") ))
write.csv( results$Range$COG_Table, file=paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results/",species_name,"COG.csv"), row.names=FALSE )

##save effective area occupied for ESP request
report = TMB::summary.sdreport(fit$parameter_estimates$SD)
ln_km2 = report[which(rownames(report)=="log_effective_area_ctl"),c('Estimate','Std. Error')]
Year <- sort(unique(fit$year_labels))
ln_km2 <- as.data.frame(cbind(ln_km2, Year))
ln_km2 <- ln_km2[which(ln_km2$Year %in% unique(fit$data_frame$t_i)),]
write.csv( ln_km2, file=paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results/",species_name,"ln_effective_area.csv"), row.names=FALSE )
