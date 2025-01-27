##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Prepare EBS catch weight data for VAST
## Authors:       Lewis Barnett (lewis.barnett@noaa.gov) 
##                Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Prepare haul-level CPUE for  
##                Kamchatka flounder in the EBS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import gapindex v3.0.2 and connect to Oracle. Make sure you are connected
##   to the internal network or VPN. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
devtools::install_github(repo = "afsc-gap-products/gapindex@v3.0.2", 
                         dependencies = TRUE)
library(gapindex)

channel <- gapindex::get_connected(check_access = F)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Species-Specific Constants. Toggle species row
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## specify whether hindcast or production phase
phase <- c("hindcast", "production")[1]

start_year <- 1991
current_year <- 2024
species_code <- 10112
species_name <- "kamchatka_flounder"
start_year_age <- 1991
plus_group <- 25

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull EBS catch and effort data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## First, pull data from the standard EBS stations
ebs_standard_data <- gapindex::get_data(year_set = start_year:current_year,
                                        survey_set = "EBS",
                                        spp_codes = species_code,
                                        pull_lengths = FALSE, 
                                        haul_type = 3, 
                                        abundance_haul = "Y",
                                        channel = channel,
                                        remove_na_strata = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate CPUE for EBS and reorder columns
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ebs_cpue <- subset(x = gapindex::calc_cpue(gapdata = ebs_standard_data),
                   select = c("SURVEY", "YEAR", "STRATUM", "HAULJOIN",
                              "LATITUDE_DD_START", "LATITUDE_DD_END",
                              "LONGITUDE_DD_START", "LONGITUDE_DD_END",
                              "SPECIES_CODE", "WEIGHT_KG", "COUNT",
                              "AREA_SWEPT_KM2", "CPUE_KGKM2", "CPUE_NOKM2"))         

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format Data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
data_geostat_index <- with(ebs_cpue,
                           data.frame(region = SURVEY,
                                      hauljoin = HAULJOIN,
                                      year = as.integer(YEAR),
                                      lon = LONGITUDE_DD_START,
                                      lat = LATITUDE_DD_START,
                                      catch_kg = WEIGHT_KG, 
                                      catch_n = COUNT,
                                      effort_km2 = AREA_SWEPT_KM2,
                                      cpue_kg_km2 = CPUE_KGKM2,
                                      cpue_n_km2 = CPUE_NOKM2))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save output. Set dir_out to the appropriate directory. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir_out <- paste0("species_specific_code/BS/", species_name, "/", phase, "/data/")
if (!dir.exists(paths = dir_out)) dir.create(path = dir_out, recursive = T)
for (ifile in c("data_geostat_index", "ebs_standard_data")) 
  saveRDS(object = get(x = ifile), 
          file = paste0(dir_out, ifile, ".RDS"))