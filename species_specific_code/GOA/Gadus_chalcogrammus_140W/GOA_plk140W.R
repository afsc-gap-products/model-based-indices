library(FishStatsUtils)
library(TMB)
library(purrr)
library(RODBC)
#library(sumfish)
library(VAST)

packageVersion('FishStatsUtils')
packageVersion('VAST')

species_name <- c("Gadus_chalcogrammus","Atheresthes_stomias")[1]
Species<-c("PLK","ATF")[1]

# Data_PLK<-read.csv("Data_PLK2021.csv")
Data_PLK <- readRDS(file = paste0(getwd(),"/species_specific_code/GOA/",species_name,"_140W/data/Data_",Species,".rds"))

#GOAGRID changed from default "Gulf of Alaska" to GOAThorsonGrid_Less700m on August 31,2020
GOAgrid <- read.csv(paste0(getwd(),"/extrapolation_grids/GOAThorsonGrid_Less700m.csv"))
input_grid <- cbind(Lat=GOAgrid$Latitude, 
                    Lon=GOAgrid$Longitude,
                    Area_km2=GOAgrid$Shape_Area/1000000)


#Per Martin Dorn's request, customize strata limits to limit extrapolation grid to west of 140
strata.limits <- data.frame(
  'STRATA' = "west_of_140W",
  'west_border' = -Inf,
  'east_border' = -140
)

settings <- make_settings( Version = "VAST_v12_0_0", 
                         n_x = 750,
                         Region = "User",
                         strata.limits = strata.limits, 
                         purpose = "index2", 
                         bias.correct = TRUE, 
                         fine_scale = TRUE, 
                         ObsModel = c(2,1),
                         use_anisotropy = TRUE)

fit_PLK <- fit_model("settings" = settings, 
                   "Lat_i" = Data_PLK$lat,
                   "Lon_i" = Data_PLK$lon, 
                   "t_i" = Data_PLK$year,
                   "c_i" = rep(0,nrow(Data_PLK)), 
                   "b_i" = Data_PLK$catch_kg,
                   "a_i" = Data_PLK$areaswept_km2, 
                   "input_grid" = input_grid, 
                   "knot_method" = 'grid',
                   "working_dir" = paste0(getwd(),"/species_specific_code/GOA/",species_name,"_140W/results"))  


plot(fit_PLK)

#saveRDS(fit_PLK, file = paste0(getwd(),"/",Species,"/","VASTfit.RDS"))
saveRDS(fit_PLK, file = paste0(getwd(),"/species_specific_code/GOA/",species_name,"_140W/results/",species_name,"VASTfit.RDS"))
