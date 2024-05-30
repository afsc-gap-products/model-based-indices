##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Prepare EBS catch weight data for VAST
## Authors:       Lewis Barnett (lewis.barnett@noaa.gov) 
##                Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Prepare haul-level CPUE for  
##                Kamchatka flounder in the EBS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import gapindex package, connect to Oracle. Make sure you are connected
##   to the internal network or VPN. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
remotes::install_github("afsc-gap-products/gapindex")
library(gapindex)
sql_channel <- gapindex::get_connected()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Select constants ----
##   start_year and current_year are used when querying RACEBASE via sumfish
##
##   For some species, we want to use all the years from 
##   start_year:current_year when calculating the ALKs but only use 
##   data from after a particular year when fitting VAST (e.g., PCod). Thus, 
##   we declare another variable called min_year to filter years at and after 
##   min_year for that purpose. By default, we set min_year <- start_year.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
which_model <- c("hindcast", "production")[1] # specify by changing index

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
ebs_data <- gapindex::get_data(year_set = start_year:current_year,
                                        survey_set = "EBS",
                                        spp_codes = species_code,
                                        pull_lengths = TRUE, 
                                        haul_type = 3, 
                                        abundance_haul = "Y",
                                        sql_channel = sql_channel,
                                        na_rm_strata = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate CPUE for EBS and reorder columns
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ebs_cpue <- subset(x = gapindex::calc_cpue(racebase_tables = ebs_data),
                   select = c("SURVEY_DEFINITION_ID", "SURVEY", "YEAR", 
                              "DESIGN_YEAR", "STRATUM", "HAULJOIN",
                              "LATITUDE_DD_START", "LATITUDE_DD_END",
                              "LONGITUDE_DD_START", "LONGITUDE_DD_END",
                              "SPECIES_CODE", "WEIGHT_KG", "COUNT",
                              "AREA_SWEPT_KM2", "CPUE_KGKM2", "CPUE_NOKM2"))          

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate Total Abundance 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ebs_popn_stratum <- gapindex::calc_biomass_stratum(racebase_tables = ebs_data,
                                                   cpue = ebs_cpue)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate numerical CPUE for a given haul/sex/length bin
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sizecomp <- gapindex::calc_sizecomp_stratum(
  racebase_tables = ebs_data,
  racebase_cpue = ebs_cpue,
  racebase_stratum_popn = ebs_popn_stratum,
  spatial_level = "haul",
  fill_NA_method = "BS")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate globally filled Age-Length Key for EBS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (start_year != start_year_age) { 
  ebs_data$size <- merge(x = ebs_data$size,
                         y = ebs_data$cruise[, c("CRUISEJOIN", "CRUISE")],
                         by = "CRUISEJOIN")
  ebs_data$size <- subset(x = ebs_data$size, 
                          CRUISE >= start_year_age * 100)
  
  ebs_data$specimen <- subset(x = ebs_data$specimen, 
                              CRUISE >= start_year_age * 100)
}

ebs_alk <- gapindex::calc_alk(racebase_tables = ebs_data, 
                              unsex = "unsex", 
                              global = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Decompose the numerical CPUE (S_ijklm) for a given haul and sex/length 
##   bin across ages. S_ijklm is the numerical CPUE of the ith station in 
##   stratum j for species k, length bin l, and sex m. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alk <- merge(x = sizecomp,
                     y = ebs_alk,
                     by = c("SURVEY", "YEAR", "SPECIES_CODE", 
                            "SEX", "LENGTH_MM"))
alk$AGE_CPUE_NOKM2 <- 
  alk$S_ijklm_NOKM2 * alk$AGE_FRAC

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Aggregate numerical CPUE across lengths for a given age/sex/haul. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
age_cpue <- rbind(
  ## Aggregate ages younger than the `plus_group`
  stats::aggregate(AGE_CPUE_NOKM2 ~ AGE + HAULJOIN + SPECIES_CODE,
                   data = alk,
                   FUN = sum,
                   drop = F,
                   subset = AGE < plus_group),
  ## Aggregate ages at or older than the `plus_group` as one age
  cbind(AGE = plus_group,
        stats::aggregate(AGE_CPUE_NOKM2 ~ HAULJOIN + SPECIES_CODE,
                         data = alk,
                         FUN = sum,
                         drop = F,
                         subset = AGE >= plus_group))
)

## Zero-fill CPUEs for missing ages. First we create a grid of all possible
## HAULJOINs and ages. But, there are hauls with positive numerical catches
## but no associated size data. These hauls will be removed.

## Hauls with zero catch
unique_hauls_cpue_zeros <- 
  sort(x = unique(x = ebs_cpue$HAULJOIN[ebs_cpue$CPUE_NOKM2 == 0]))
## Hauls with postive count data
unique_hauls_cpue_pos <- 
  sort(x = unique(x = ebs_cpue$HAULJOIN[ebs_cpue$CPUE_NOKM2 > 0]))
## Hauls with size data
unique_hauls_sizecomp <- 
  sort(x = unique(x = sizecomp$HAULJOIN))

every_combo_of_ages <- 
  expand.grid(HAULJOIN = c(unique_hauls_cpue_zeros,
                           unique_hauls_cpue_pos[unique_hauls_cpue_pos %in% 
                                                   unique_hauls_sizecomp]),
              SPECIES_CODE = species_code,
              AGE = min(age_cpue$AGE, na.rm = TRUE):plus_group)

## Merge the age_cpue onto the `every_combo_of_ages` df using "HAULJOIN", 
## "SPECIES_CODE", and "AGE" as a composite key. all.x = TRUE will reveal
## missing ages denoted by an NA AGE_CPUE_NOKM2 value
every_combo_of_ages <- merge(x = every_combo_of_ages,
                             y = age_cpue,
                             by = c("HAULJOIN", "SPECIES_CODE", "AGE"),
                             all.x = TRUE)

## Turn NA AGE_CPUE_NOKM2 values to zero
every_combo_of_ages$AGE_CPUE_NOKM2[
  is.na(x = every_combo_of_ages$AGE_CPUE_NOKM2)
] <- 0

## Append the latitude and longitude information from the `ebs_cpue` df
## using "HAULJOIN" as a key. 
age_cpue <- merge(x = every_combo_of_ages,
                  y = ebs_cpue[, c("HAULJOIN", "SURVEY", "YEAR",
                                       "LATITUDE_DD_START",
                                       "LONGITUDE_DD_START")],
                  by = "HAULJOIN")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format VAST Data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
data_geostat_biomass_index <- with(ebs_cpue,
                                   data.frame(Hauljoin = HAULJOIN,
                                              Region = SURVEY,
                                              Catch_KG = CPUE_KGKM2,
                                              Year = YEAR,
                                              Vessel = "missing", 
                                              AreaSwept_km2 = 1,
                                              Lat = LATITUDE_DD_START,
                                              Lon = LONGITUDE_DD_START,
                                              Pass = 0))

data_geostat_numerical_index <- with(ebs_cpue,
                                     data.frame(Hauljoin = HAULJOIN,
                                                Region = SURVEY,
                                                Catch_KG = CPUE_NOKM2,
                                                Year = YEAR,
                                                Vessel = "missing", 
                                                AreaSwept_km2 = 1,
                                                Lat = LATITUDE_DD_START,
                                                Lon = LONGITUDE_DD_START,
                                                Pass = 0))

data_geostat_agecomps <- with(age_cpue, 
                              data.frame(Hauljoin = HAULJOIN,
                                         Region = SURVEY,
                                         Catch_KG = AGE_CPUE_NOKM2,
                                         Year = YEAR,
                                         Vessel = "missing",
                                         Age = AGE,
                                         AreaSwept_km2 = 1,
                                         Lat = LATITUDE_DD_START,
                                         Lon = LONGITUDE_DD_START,
                                         Pass = 0)) 

if (start_year != start_year_age) {
  data_geostat_agecomps <- subset(x = data_geostat_agecomps,
                                  subset = Year >= start_year_age)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save output. Set dir_out to the appropriate directory. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir_out <- paste0(getwd(),"/species_specific_code/BS/", 
                  species_name, "/", which_model, "/data/")
if (!dir.exists(paths = dir_out)) dir.create(path = dir_out, recursive = T)
for (ifile in c("data_geostat_biomass_index", 
                "data_geostat_numerical_index",
                "data_geostat_agecomps", 
                "sizecomp", "alk", "ebs_data")) 
  saveRDS(object = get(x = ifile), 
          file = paste0(dir_out, ifile, ".RDS"))