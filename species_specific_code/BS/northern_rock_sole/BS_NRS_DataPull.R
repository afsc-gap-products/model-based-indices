##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Prepare EBS/NBS catch weight data
## Authors:       Lewis Barnett (lewis.barnett@noaa.gov) 
##                adapted from 
##                Zack Oyafuso (zack.oyafuso@noaa.gov)
##                Emily Markowitz (emily.markowitz@noaa.gov)
##                Jason Conner (jason.conner@noaa.gov)
## Description:   Prepare haul-level CPUE for the 
##                EBS and NBS for northern rock sole
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import sumfish ----
##   Installation instructions: https://github.com/afsc-gap-products/sumfish
##   Setup Oracle Account
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
library(sumfish)
sumfish::getSQL()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Select species and years ----
##   start_yer and current_year are used when querying RACEBASE via sumfish
##
##   For some species, we want to use all the years from 
##   start_year:current_year when calculating the ALKs but only use 
##   data from after a particular year when fitting VAST (e.g., PCod). Thus, 
##   we declare another variable called min_year to filter years at and after 
##   min_year for that purpose. By default, we set min_year <- start_year.
##
##   plus_group is used for the age composition calculations, where ages at or
##   older than the plus group are grouped as one group. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_name <- "northern_rock_sole"
species_code <- 10261
start_year <- 1982
current_year <- 2023
min_year <- start_year

which_model <- c("hindcast", "production")[2] # specify by changing index

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create directory to store data products
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
folder <- paste0(getwd(),"/species_specific_code/BS/",species_name)
dir.create(folder)

res_dir <- paste0(getwd(), "/species_specific_code/BS/", 
                  species_name, "/", which_model, "/", "data/")
if(!dir.exists(paste0(res_dir))) 
  dir.create(path = paste0(res_dir), recursive = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull EBS database ----
##   Query database for the standard EBS survey haul and catch data.
##   Append nonstandard well-performing EBS hauls that occurred in 1994, 2001, 
##   2005, and 2006. These tows are not used in the design-based estimates but
##   becasue they follow standard procedures, can be included in a model-
##   based estiamte.
##
##   Notes: haul_type == 3 = standard bottom sample (preprogrammed station)
##          PERFORMANCE 0 is a good performance code, > 0 are satisfactory and
##          < 0 are bad performance codes for a particular haul.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
EBS <- sumfish::getRacebase(year = c(start_year, current_year), 
                            survey = 'EBS_SHELF', 
                            speciesCode = species_code) 

EBS_nonstandard_hauls <- EBS$haul_other %>%
  filter(CRUISE %in% c(200101, 199401, 200501, 200601) 
         & PERFORMANCE >= 0
         & HAUL_TYPE == 3
         & !is.na(STRATUM))

EBS_nonstandard_catch <- EBS$catch_other %>%
  filter(HAULJOIN %in% EBS_nonstandard_hauls$HAULJOIN)

EBS$catch <- bind_rows(EBS$catch, EBS_nonstandard_catch)
EBS$haul <- bind_rows(EBS$haul, EBS_nonstandard_hauls)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull NBS database ----
##   Query database for the standard NBS survey haul and catch data.
##   Append nonstandard well-performing NBS hauls that occurred in 2018 that
##   are in the EBS dataset. The way that this is queried is by selecting
##   records from the Bering Sea (BS), CRUISE is 201801 (Year 2018),
##   and haul_type == 13, meaning it is an index sample tow.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

## Standard NBS Hauls
NBS <- sumfish::getRacebase(year = c(start_year, current_year), 
                            survey = 'NBS_SHELF', 
                            speciesCode = species_code) 

## 2018 NBS cruise
NBS18.cruise <- dplyr::filter(EBS$cruise, CRUISE == 201801)
NBS18.haul <- sumfish::getSQL("select * from racebase.haul where region = 'BS' and cruise = 201801 and haul_type = 13 and performance>=0")
NBS18.haul$STRATUM <- '99'
NBS18.catch <- sumfish::getSQL("select * from racebase.catch where hauljoin in (select hauljoin from racebase.haul where region = 'BS' and cruise = 201801 and haul_type = 13 and performance>=0)")
NBS18.stratum <- data.frame(SURVEY = 'NBS_18', 
                            STRATUM = '99', 
                            STRATUM_AREA = 158286)
NBS18.length <- sumfish::getSQL("select * from racebase.length where hauljoin in (select hauljoin from racebase.haul where region = 'BS' and cruise = 201801 and haul_type = 13 and performance>=0)")
NBS18.specimen <- sumfish::getSQL("select * from racebase.specimen where hauljoin in (select hauljoin from racebase.haul where region = 'BS' and cruise = 201801 and haul_type = 13 and performance>=0)")

## Combine all the NBS18 data types into a list
NBS18 <- list(cruise = NBS18.cruise,
              haul = NBS18.haul,
              catch = NBS18.catch,
              species = NBS$species,
              stratum = NBS18.stratum,
              length = NBS18.length,
              specimen = NBS18.specimen) 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate haul-level CPUEs ----
##   Fill in zeros for hauls that did not encounter the species, 
##   Bind EBS, NBS and NBS18 haul data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
sumEBS <- sumfish::sumHaul(EBS) %>%
  dplyr::mutate(REGION = "EBS")
sumNBS <- sumfish::sumHaul(NBS) %>%
  dplyr::mutate(REGION = "NBS")
sum18 <- sumfish::sumHaul(NBS18) %>%
  dplyr::mutate(STRATUM = as.character(STRATUM),
                REGION = "NBS")

weightAll <- sumAll <- dplyr::bind_rows(sumEBS, sumNBS,sum18) %>%
  dplyr::filter(SPECIES_CODE %in% species_code, 
                YEAR >= min_year,
                !is.na(EFFORT))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Test for missing values
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
check_hauls <- c(EBS$haul$HAULJOIN, NBS$haul$HAULJOIN, NBS18$haul$HAULJOIN)
ifelse(test = length(sumAll$HAULJOIN[!sumAll$HAULJOIN %in% check_hauls]) == 0,
       yes = "No missing hauls",
       no = paste0(length(sumAll$HAULJOIN[!sumAll$HAULJOIN %in% check_hauls]),
                   " missing hauls. Check code."))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format VAST Data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
data_geostat_index <- 
  with(sumAll,
       data.frame(Region = REGION,
                  Catch_KG = wCPUE, #cpue units: kg per hectare
                  Year = YEAR,
                  Vessel = "missing", 
                  AreaSwept_km2 = 0.01, # converts cpue units to: kg per km^2
                  Lat = START_LATITUDE,
                  Lon = START_LONGITUDE,
                  Pass = 0 ))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save output ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

## Index data
write_rds(x = weightAll, 
          file = paste0(res_dir, "EBS_NBS_Index.RDS"))
write_rds(x = data_geostat_index, 
          file = paste0(res_dir, "data_geostat_index.RDS"))

## Strata data
strata <- dplyr::bind_rows(EBS$stratum, NBS$stratum, NBS18$stratum) 
write_rds(x = strata, 
          file = paste0(res_dir, "EBS_NBS_strata.RDS"))

## Raw data
write_rds(x = list(EBS = EBS, NBS = NBS, NBS18 = NBS18), 
          file = paste0(res_dir, "raw_data.RDS"))
