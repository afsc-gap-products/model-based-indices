### Cecilia O'Leary, cecilia.oleary@noaa.gov; Madison Hall, madison.hall@noaa.gov
### code to produce model-based indices using VAST software implemented in R
### created March 2020
### 2020 settings: VAST v3.3.0, FishStatsUtils v2.5.0, cpp VAST_v8_2_0
### 2021 settings (updated packages): Rv4.0.2 VAST v3.6.1, FishStatsUtils v2.8.0, cpp VAST_v12_0_0

## Load packages and set working directory
library(TMB)               
library(VAST)

## check to ensure you're runing the correct package versions
packageVersion('FishStatsUtils')
packageVersion('VAST')

## how to check your session info (see R version and all package versions)
sessionInfo()

## Set species 
Species <- c("Gadus_macrocephalus","Sebastes_variabilis","Sebastes_polyspinis", "Sebastes alutus","Gadus chalcogrammus", "Lepidopsetta polyxystra","Lepidopsetta bilineata","Hippoglossoides elassodon","Atheresthes stomias")[3] 
### 3 is for Northern
## c(P. cod, dusky rockfish, northern rockfish, POP, pollock, 
## northern rock sole, southern rock sole, flathead sole, arrowtooth)

# Set up folder to store species specific results
folder <- paste0(getwd(),"/",Species)
dir.create(folder)

## ## 
## Load the data for VAST
Ndata<-readRDS("Data_Geostat_Sebastes_polyspinis.RDS")
Ndata<-Ndata[Ndata$Year>1984, ]## Ben Williams requested 1984 onwards.

## Define settings


settings = make_settings( Version = "VAST_v12_0_0", #.cpp version, not software #e.g., "VAST_v12_0_0"
                          n_x = 750, #knots aka spatial resolution of our estimates
                          #Ben asked for 750 on Oct 4
                          Region = "User", #Region = "gulf_of_alaska" , go to ?make_settings for other built in extrapolation grids
                          purpose = "index2", #changes default settings
                          ObsModel= c(2,1)#, this is the default gamma obs error, other two options: #c(1,1) #c(10,2) 
) 

## ## 
## Import extrapolation grid
GOAgrid <- read.csv(file= paste0(getwd(),"/Extrapolation_Grids/GOAThorsonGrid_Less700m.csv"))

input_grid <- cbind(Lat = GOAgrid$Lat,
                    Lon = GOAgrid$Lon,
                    Area_km2 = GOAgrid$Shape_Area/1000000)  # Extrapolation grid area is in m^2 & is converted to km^2 with this line
gc() #garbage collector

## ## 
## Run model
fit <- fit_model( "settings"= settings, #all of the settings we set up above
                  "Lat_i"= Ndata[,'Lat'], #latitude of observation
                  "Lon_i"= Ndata[,'Lon'],  #longitude of observation
                  "t_i"= Ndata[,'Year'], #time for each observation
                  "b_i"= Ndata[,'Catch_KG'], #in kg, raw catch or in CPUE per tow
                  "a_i"= Ndata[,'AreaSwept_km2'], #sampled area for eahc observation
                  "v_i"= Ndata[,'Vessel'], #ok to leave in because it's all "missing" in data, so NO vessel effects
                  "input_grid"= input_grid, #only needed if you have a user input extrapolation grid, which we do for GOA
                  "optimize_args" =list("lower"=-Inf,"upper"=Inf), #TMB argument (?fit_tmb) tath can be used if you're having optimization issues, shouldn't need to change
                  "working_dir" = paste0(getwd(),"/",Species,"/"))

## ## 
## Plot results
plot( fit )
#if you need to link VAST as a .dll:
dyn.load(dynlib("VAST_v12_0_0"))


## ## 
## save the VAST model
saveRDS(fit,file = paste0(getwd(),"/",Species,"/",Species,"VASTfit.RDS"))

