
#setwd("C:/VAST_GOA_2019/PLK2020_140W")
setwd("C:/Paul/VAST_2021")
library(FishStatsUtils)
library(TMB)
library(purrr)
library(RODBC)
#library(sumfish)
library(VAST)

packageVersion('FishStatsUtils')
packageVersion('VAST')

#Species<-c("Gadus_chalcogrammus","Atheresthes_stomias")
Species<-c("Pollock","ATF")[1]

#Set up folder to store results
folder<-paste0(getwd(),"/",Species)
dir.create(folder)


#SKIP TO: "START HERE AFTER setwd and library functions

#This function (Ned L) extracts relevant data from racebase.haul,goa.cpue, and racebase.species tables
get_data_geostat <- function(channel = NA, region = "GOA", sp.code = 21740){
  
  if(is.na(channel)){
    require(RODBC)
    channel <- odbcConnect(dsn = "AFSC", uid = "vszalay", pwd = "PlatinumSnapper$4", believeNRows = FALSE)
    close.channel = TRUE
  }else{
    close.channel <- FALSE
  }
  
  sqry <- paste0("select a.vessel, start_latitude lat, start_longitude lon, a.cruise||'_'||b.haul||'_'||b.vessel towid, 
                 floor(a.cruise/100) year, weight catch_kg, effort areaswept_km2, c.common_name
                 from racebase.haul a, goa.cpue b, racebase.species c
                 where a.hauljoin = b.hauljoin and a.region = '", region, "' and a.cruise >= 198401 and b.species_code = ", 
                 sp.code, " and b.species_code = c.species_code and a.start_longitude < -140")
  
  dat <- sqlQuery(channel = channel, query = sqry, rows_at_time = 1, errors = TRUE)
  names(dat) <- tolower(names(dat))
  
  if(close.channel)close(channel)
  
  dat
  
}

#Run Ned's function above for indicated species and create a df with required columns
Data<-get_data_geostat(region="GOA",sp.code = 21740)
Data_PLK<-dplyr::select(Data, "catch_kg","year","vessel","areaswept_km2","lat","lon")

#Add a column 'Pass' containing 0's to match Thorson's example
Pass<-rep(0,nrow(Data_PLK))
Data_PLK<-cbind(Data_PLK,Pass)
Data_PLK<-na.omit(Data_PLK)


#Create csv file of Data_PLK and save in working directory
write.csv(Data_PLK,file="C:/Paul/VAST_2021/Data_PLK2021.csv")

#using Jason's sumfish package, extract relevant tables from Racebase
###data<-sumHaul(getRacebase(year=c(1984,2017),'GOA'))
###data.red<-data[,c(1,6,10,14:15,19:22)]
###data.red.pop<-data.red[data.red$SPECIES_CODE==30060,]
#x<-sort(unique(data.red.plk$YEAR))

###Data_POP<-dplyr::select(data.red.pop,WEIGHT,YEAR,VESSEL,EFFORT,START_LATITUDE,START_LONGITUDE)
###names(Data_POP)<-c("Catch_KG","Year","Vessel","AreaSwept_km2","Lat","Lon")

#Add a column 'Pass' containing 0's to match Thorson's example
###Pass<-rep(0,nrow(Data_POP))
###Data_POP<-cbind(Data_POP,Pass)        


#################################################################################################################################
#START HERE AFTER setwd and library functions above
#################################################################################################################################
#In lieu of running Ned's Oracle query above, start here if that's already been done and the output is a file called Data_PLK.csv
Data_PLK<-read.csv("Data_PLK2021.csv")

#GOAGRID changed from default "Gulf of Alaska" to GOAThorsonGrid_Less700m on August 31,2020
GOAgrid<-read.csv("GOAThorsonGrid_Less700m.csv")
input_grid<-cbind(Lat=GOAgrid$Latitude, Lon=GOAgrid$Longitude,
                  Area_km2=GOAgrid$Shape_Area/1000000)


#Per Martin Dorn's request, customize strata limits to limit extrapolation grid to west of 140
strata.limits <- data.frame(
  'STRATA' = "west_of_140W",
  'west_border' = -Inf,
  'east_border' = -140
)

#setwd("C:/Paul/VAST2021_test/PLK/output")


#settings<-make_settings(n_x=100, Region=example$Region, purpose = "index",strata.limits = example$strata.limits, bias.correct = FALSE)
settings<-make_settings( Version="VAST_v12_0_0", n_x=750,Region="User",strata.limits = strata.limits, purpose = "index2", bias.correct = TRUE, fine_scale = TRUE, ObsModel = c(2,1),use_anisotropy = TRUE)
fit_PLK<-fit_model("settings"=settings, "Lat_i"=Data_PLK$lat,
                   "Lon_i"=Data_PLK$lon, "t_i"=Data_PLK$year,
                   "c_i"=rep(0,nrow(Data_PLK)), "b_i"=Data_PLK$catch_kg,
                   "a_i"=Data_PLK$areaswept_km2, "input_grid"=input_grid, knot_method='grid',
                   "working_dir"= paste0(getwd(),"/",Species,"/"))  #"v_i"=Data_PLK$vessel,


plot(fit_PLK)

saveRDS(fit_PLK, file = paste0(getwd(),"/",Species,"/","VASTfit.RDS"))
