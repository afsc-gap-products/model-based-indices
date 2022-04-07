###############################################################################
## Project:         GOA Groundfish CPUE Data Synthesis 2021
##
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
##                  adapted from Lewis Barnett (lewis.barnett@noaa.gov)
##
## Description:     Create CPUE dataset used for VAST for 
##                  Pacific cod and walleye pollock
##                  Get AK groundfish bottom trawl survey data for 3 primary 
##                  surveys Result is cleaned cpue (kg/km^2) by haul, with 
##                  zeros included
##
## Notes:           sumfish package is used to pull haul data
##                  R version 3.6.1
###############################################################################

##################################################
#### Import Packages
##################################################
library(sumfish)
library(dplyr)

##################################################
#### Log into ROracle and pull GoA haul data from 1984 - 2021
##################################################
sumfish::setUser()
goa_all <- getRacebase(year=c(2021, 1984), survey = "GOA")
cruise <- goa_all$cruise
catch <- goa_all$catch
species <- goa_all$species
haul_with_zeros <- sumHaul(goa_all)
haul <- goa_all$haul

##################################################
#### Subset haul data for pollock and PCod
##################################################
spp_idx <- species$COMMON_NAME %in% c("walleye pollock", "Pacific cod")
spp_code <- species$SPECIES_CODE[spp_idx]
data <- subset(x = haul_with_zeros, subset = SPECIES_CODE %in% spp_code)

##################################################
#### join haul data to get coordinates
##################################################
haul <- cbind(haul, 
              geosphere::midPoint(cbind(haul$START_LONGITUDE, 
                                        haul$START_LATITUDE), 
                                  cbind(haul$END_LONGITUDE, 
                                        haul$END_LATITUDE))) 
haul <- haul %>% select(HAULJOIN, LATITUDE = lat, LONGITUDE = lon)
data <- inner_join(data, haul)

##################################################
#### Convert units of effort from hectares to square kilometers
##################################################
data$EFFORT <- data$EFFORT * 0.01

##################################################
#### Subset data to those west of -140 degrees longitude
##################################################
boundary <- -140
data <- data %>% filter(LONGITUDE <= boundary)

##################################################
#### Write output
##################################################
write.csv(x = data, 
          file = "processed_data/AK_BTS_GOA_PCod_pollock_140.csv",
          row.names = F)

write.csv(x = species,
          file = "processed_data/species.csv",
          row.names = F)
