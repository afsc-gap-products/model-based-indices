species_name <- c("Gadus_chalcogrammus","Atheresthes_stomias")[1]
Species<-c("Pollock","ATF")[1]

#This function (Ned L) extracts relevant data from racebase.haul,goa.cpue, and racebase.species tables
get_data_geostat <- function(region = "GOA", sp.code = 21740){
  
  PKG <- c("RODBC")
  for (p in PKG) {
    if(!require(p,character.only = TRUE)) {  
      install.packages(p)
      require(p,character.only = TRUE)}
  }
  
  # odbcGetInfo(channel)
  source("R/get_connected.R")
  
  sqry <- paste0("select a.vessel, start_latitude lat, start_longitude lon, a.cruise||'_'||b.haul||'_'||b.vessel towid, 
                 floor(a.cruise/100) year, weight catch_kg, effort areaswept_km2, c.common_name
                 from racebase.haul a, goa.cpue b, racebase.species c
                 where a.hauljoin = b.hauljoin and a.region = '", region, "' and a.cruise >= 198401 and b.species_code = ", 
                 sp.code, " and b.species_code = c.species_code and a.start_longitude < -140")
  
  dat <- sqlQuery(channel = channel, query = sqry, rows_at_time = 1, errors = TRUE)
  names(dat) <- tolower(names(dat))
  
  dat
  
}

#Run Ned's function above for indicated species and create a df with required columns
Data <- get_data_geostat(region="GOA",sp.code = 21740)
Data_PLK <- dplyr::select(Data, "catch_kg","year","vessel","areaswept_km2","lat","lon")

#Add a column 'Pass' containing 0's to match Thorson's example
Pass <- rep(0,nrow(Data_PLK))
Data_PLK <- cbind(Data_PLK,Pass)
Data_PLK <- na.omit(Data_PLK)


#Create csv file of Data_PLK and save in working directory
#write.csv(Data_PLK,file="C:/Paul/VAST_2021/Data_PLK2021.csv")
saveRDS(Data_PLK, paste0(getwd(),"/species_specific_code/GOA/",species_name,"_140W/data/Data_PLK.rds"))
