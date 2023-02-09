##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Catch and Effort Data pull
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Synethesize zero-filled CPUE (kg/km2) dataset used for VAST
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(RODBC)
library(getPass)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set constants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
requested_start_year <- 1990
species_code <- c("Gadus_macrocephalus" = 21720, 
                  "Gadus_chalcogrammus" = 21740)[1]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Connect to Oracle. Make sure to connect to network or VPN
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("R/get_connected.R")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull Haul data with abundance haul == Y, haul_type == 3, Performance >= 0
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
haul <- RODBC::sqlQuery(channel = channel, 
                        query = paste("SELECT * FROM RACEBASE.HAUL WHERE",
                                      "REGION = 'GOA' AND",
                                      "ABUNDANCE_HAUL = 'Y' AND",
                                      "HAUL_TYPE = 3 AND",
                                      "PERFORMANCE >= 0"))
haul$YEAR = as.numeric(x = substring(text = haul$START_TIME, 
                                     first = 1, last = 4))
haul <- subset(x = haul, subset = YEAR >= requested_start_year)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull catch data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
catch <- RODBC::sqlQuery(channel = channel, 
                         query = paste("SELECT * FROM RACEBASE.CATCH WHERE",
                                       "REGION = 'GOA' AND",
                                       "SPECIES_CODE = ", species_code))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge catch data to the haul df using HAULJOIN as the key
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
catch_effort <- merge(x = haul, 
                      y = catch, 
                      by = "HAULJOIN", 
                      all.x = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Convert NA Weights to 0
##   Effort is calculated as the product of distance fished (km) * net_width
##          (units of m, converted to km) 
##   Calculate CPUE
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
catch_effort$WEIGHT[is.na(x = catch_effort$WEIGHT)] <- 0
catch_effort$EFFORT <- with(catch_effort, DISTANCE_FISHED * NET_WIDTH * 0.001)
catch_effort$wCPUE <- with(catch_effort, WEIGHT / EFFORT )

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create VAST input df
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Data_Geostat <-
  with(catch_effort,
       data.frame(Catch_KG = wCPUE, 
                  Year = YEAR,
                  Vessel = "missing",
                  # area swept is 1 when using wCPUE instead of obs weight
                  AreaSwept_km2 = 1, 
                  Lat = START_LATITUDE,
                  Lon = START_LONGITUDE,
                  Pass = 0) )

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save object as both an RDS and RData file.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!dir.exists(paths = paste0("species_specific_code/GOA/",
                               names(species_code), "/data/")))
  dir.create(path = paste0("species_specific_code/GOA/",
                           names(species_code), "/data/"), 
             recursive = TRUE)

saveRDS(object = Data_Geostat,
        file = paste0("species_specific_code/GOA/",
                      names(species_code), "/data/Data_Geostat_", 
                      names(species_code), ".rds"))

save(Data_Geostat, 
     file = paste0("species_specific_code/GOA/",
                   names(species_code), "/data/Data_Geostat_", 
                   names(species_code), ".RData"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   For ESP runs, create a subsetted df for those data west of -140 Lon
##   Save object as both an RDS and RData file.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Data_Geostat_w140 <- subset(x = Data_Geostat,
                            subset = Lon < -140)

if (!dir.exists(paths = paste0("species_specific_code/GOA/Pcod_WPlk_ESP/",
                               names(species_code), "_w140/data/")))
  dir.create(path = paste0("species_specific_code/GOA/Pcod_WPlk_ESP/",
                            names(species_code), "_w140/data/"), 
             recursive = TRUE)

saveRDS(object = Data_Geostat_w140,
        file = paste0("species_specific_code/GOA/Pcod_WPlk_ESP/",
                      names(species_code), "_w140/data/Data_Geostat_", 
                      names(species_code), "_w140.rds"))

save(Data_Geostat_w140, 
     file = paste0("species_specific_code/GOA/Pcod_WPlk_ESP/",
                   names(species_code), "_w140/data/Data_Geostat_", 
                   names(species_code), "_w140.RData"))
