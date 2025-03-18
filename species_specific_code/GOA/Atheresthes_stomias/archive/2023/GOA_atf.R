library(FishStatsUtils)
library(TMB)
library(purrr)
library(RODBC)
#library(sumfish)
library(VAST)

packageVersion('FishStatsUtils')
packageVersion('VAST')

species_name <- c("Gadus_chalcogrammus","Atheresthes_stomias")[2]
Species <- c("Pollock","ATF")[2]

#load data
Data_ATF <- readRDS(file = paste0(getwd(),"/species_specific_code/GOA/",species_name,"/data/Data_",Species,".rds"))

#GOAGRID changed from default "Gulf of Alaska" to GOAThorsonGrid_Less700m on August 31,2020
GOAgrid<-read.csv(paste0(getwd(),"/extrapolation_grids/GOAThorsonGrid_Less700m.csv"))
input_grid<-cbind(Lat=GOAgrid$Latitude, Lon=GOAgrid$Longitude,
                  Area_km2=GOAgrid$Shape_Area/1000000)

#settings<-make_settings(n_x=100, Region=example$Region, purpose = "index",strata.limits = example$strata.limits, bias.correct = FALSE)
settings <- make_settings( Version = "VAST_v12_0_0", 
                         n_x = 750,
                         Region = "User", 
                         purpose = "index2", 
                         bias.correct = TRUE, 
                         fine_scale = TRUE, 
                         ObsModel = c(2,1),
                         use_anisotropy = TRUE)

fit_ATF <- fit_model("settings" = settings, 
                     "Lat_i" = Data_ATF$lat,
                     "Lon_i" = Data_ATF$lon, 
                     "t_i" = Data_ATF$year,
                     "c_i" = rep(0,nrow(Data_ATF)), 
                     "b_i" = Data_ATF$catch_kg,
                     "a_i" = Data_ATF$areaswept_km2, 
                     "input_grid" = input_grid, 
                     "knot_method"='grid',
                     "working_dir" = paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results"))  


plot(fit_ATF)

#Save the model fit to the ATF folder inside the working directory
#saveRDS(fit_ATF, file = paste0(getwd(),"/",Species,"/","VASTfit.RDS"))
saveRDS(fit_ATF, file = paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results/",species_name,"VASTfit.RDS"))


###Save COG (center of gravity) for ESP request
dyn.load(dynlib("ATF/VAST_v12_0_0"))
results = plot(fit_ATF, n_cells=200^2)
write.csv( results$Range$COG_Table, file=paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results/",species_name,"COG.csv"), row.names=FALSE )

#Save effective area occupied for ESP request
report=TMB::summary.sdreport(fit_ATF$parameter_estimates$SD)
In_km2=report[which(rownames(report)=="log_effective_area_ctl"),c('Estimate','Std. Error')]
Year<-sort(unique(fit_ATF$year_labels))
In_km2<-as.data.frame(cbind(In_km2, Year))
In_km2<-In_km2[which(In_km2$Year %in% unique(fit_ATF$data_frame$t_i)),]
write.csv( ln_km2, file=paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results/",species_name,"ln_effective_area.csv"), row.names=FALSE )