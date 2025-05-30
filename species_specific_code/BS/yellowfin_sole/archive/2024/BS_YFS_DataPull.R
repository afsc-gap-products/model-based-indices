##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       EBS/NBS data pull for ModSquad VAST input
## Description:   via gapindex R package
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import gapindex package, connect to Oracle. Make sure you are connected
##   to the internal network or VPN. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
remotes::install_github("afsc-gap-products/gapindex")
library(gapindex)
sql_channel <- gapindex::get_connected()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Species-Specific Constants. Toogle species row
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_info <- data.frame(species_name = c("yellowfin_sole", "Pacific_cod"),
                           species_code = c(10210, 21720),
                           start_year = 1982,
                           current_year = 2024,
                           plus_group = c(20, 12), 
                           start_year_age = c(1982, 1994))[1, ]

## Set constants
which_model <- c("hindcast", "production")[2]
start_year <- species_info$start_year
current_year <- species_info$current_year
species_code <- species_info$species_code
species_name <- species_info$species_name
start_year_age <- species_info$start_year_age
plus_group <- species_info$plus_group

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull EBS catch and effort data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## First, pull data from the standard EBS stations
ebs_standard_data <- gapindex::get_data(year_set = start_year:current_year,
                                        survey_set = "EBS",
                                        spp_codes = species_code,
                                        pull_lengths = TRUE, 
                                        haul_type = 3, 
                                        abundance_haul = "Y",
                                        sql_channel = sql_channel,
                                        na_rm_strata = TRUE)

## Next, pull data from hauls that are not included in the design-based index
## production (abundance_haul == "N") but are included in VAST. By default, 
## the gapindex::get_data() function will filter out hauls with negative 
## performance codes (i.e., poor-performing hauls).
ebs_other_data <- gapindex::get_data(year_set = c(1994, 2001, 2005, 2006),
                                     survey_set = "EBS",
                                     spp_codes = species_code,
                                     pull_lengths = TRUE, 
                                     haul_type = 3, 
                                     abundance_haul = "N",
                                     sql_channel = sql_channel, 
                                     na_rm_strata = TRUE)

## Combine the EBS standard and EBS other data into one list. 
ebs_data <- list(
  ## Some cruises are shared between the standard and other EBS cruises, 
  ## so the unique() wrapper is there to remove duplicate cruise records. 
  cruise = unique(rbind(ebs_standard_data$cruise,
                        ebs_other_data$cruise)),
  haul = rbind(ebs_standard_data$haul,
               ebs_other_data$haul),
  catch = rbind(ebs_standard_data$catch,
                ebs_other_data$catch),
  size = rbind(ebs_standard_data$size,
               ebs_other_data$size),
  species = ebs_standard_data$species,
  specimen = rbind(ebs_standard_data$specimen,
                   ebs_other_data$specimen),
  strata = ebs_standard_data$strata)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull NBS catch and effort data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Pull data from the standard NBS stations
nbs_standard_data <- gapindex::get_data(year_set = start_year:current_year,
                                        survey_set = "NBS",
                                        spp_codes = species_code,
                                        pull_lengths = TRUE, 
                                        haul_type = 3, 
                                        abundance_haul = "Y",
                                        sql_channel = sql_channel)

## Pull data from miscellaneous stations in 1985, 1988, and 1991 that were 
## sampled in the NBS but are not a part of the standard stations used in the
## design-based index production. These cruises have a different survey 
## definition IDs from the EBS (98), NBS (143) and BSS (78), so they will not
## come up in gapindex::get_data(). Thus, the next set of code consists of
## hard-coded SQL pulls from RACEBASE.
nbs_other <- list()

nbs_other$cruise <-
  data.frame(YEAR = c(1988, 1985, 1991),
             SURVEY_DEFINITION_ID = 112,
             SURVEY = "NBS",
             RODBC::sqlQuery(
               channel = sql_channel,
               query = "SELECT CRUISEJOIN, CRUISE
                          FROM RACEBASE.CRUISE 
                          WHERE CRUISE IN (198502, 198808, 199102) 
                          AND REGION = 'BS'"),
             VESSEL_ID = NA,
             VESSEL_NAME = NA,
             DESIGN_YEAR = 2022)

nbs_other$haul <- 
  RODBC::sqlQuery(channel = sql_channel,
                  query = "SELECT *
                             FROM RACEBASE.HAUL 
                             WHERE CRUISE IN (198502, 198808, 199102) 
                             AND HAUL_TYPE = 3 
                             AND PERFORMANCE >= 0 
                             AND ABUNDANCE_HAUL = 'Y'")
nbs_other$catch <- 
  RODBC::sqlQuery(
    channel = sql_channel,
    query = paste("SELECT HAULJOIN, SPECIES_CODE, WEIGHT, NUMBER_FISH
                     FROM RACEBASE.CATCH 
                     WHERE SPECIES_CODE IN", 
                  gapindex::stitch_entries(species_code),
                  "AND HAULJOIN IN",
                  gapindex::stitch_entries(nbs_other$haul$HAULJOIN)))

nbs_other$size <-
  RODBC::sqlQuery(
    channel = sql_channel,
    query = paste("SELECT CRUISEJOIN, HAULJOIN, SPECIES_CODE, LENGTH, 
                     SEX, FREQUENCY 
                     FROM RACEBASE.LENGTH 
                     WHERE SPECIES_CODE IN", 
                  gapindex::stitch_entries(species_code),
                  "AND HAULJOIN IN",
                  gapindex::stitch_entries(nbs_other$haul$HAULJOIN)))

nbs_other$specimen <-
  RODBC::sqlQuery(
    channel = sql_channel,
    query = paste("SELECT SPECIES_CODE, CRUISEJOIN, HAULJOIN, REGION,
                     VESSEL, CRUISE, HAUL, SPECIMENID, LENGTH, SEX, WEIGHT, AGE
                     FROM RACEBASE.SPECIMEN 
                     WHERE SPECIES_CODE IN", 
                  gapindex::stitch_entries(species_code),
                  "AND HAULJOIN IN",
                  gapindex::stitch_entries(nbs_other$haul$HAULJOIN)))

## Pull data from a 2018 NBS survey that is defined as an EBS survey but just
## with a different haul_type value. ABUNDANCE_TYPE for these hauls is 'N' 
## but by default, the gapindex::get_data() function will filter out hauls 
## with negative performance codes (i.e., poor-performing hauls).
nbs18_data <- gapindex::get_data(year_set = 2018,
                                 survey_set = "EBS",
                                 spp_codes = species_code,
                                 pull_lengths = TRUE, 
                                 haul_type = 13, 
                                 abundance_haul = "N",
                                 sql_channel = sql_channel)
nbs18_data$cruise$SURVEY <- "NBS"
nbs18_data$cruise$SURVEY_DEFINITION_ID <- 143

## Combine all the NBS data into one list. 
nbs_data <- 
  list(cruise = rbind(nbs_standard_data$cruise[,names(x = nbs18_data$cruise)],
                      nbs_other$cruise,
                      nbs18_data$cruise),
       catch = rbind(nbs_standard_data$catch,
                     nbs_other$catch,
                     nbs18_data$catch),
       haul = rbind(nbs_standard_data$haul,
                    nbs_other$haul,
                    nbs18_data$haul),
       size = rbind(nbs_standard_data$size,
                    nbs_other$size,
                    nbs18_data$size),
       specimen = rbind(nbs_standard_data$specimen,
                        nbs_other$specimen,
                        nbs18_data$specimen),
       species = nbs_standard_data$species,
       strata = nbs_standard_data$strata)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate CPUE for EBS and NBS, rbind, and reorder columns
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ebs_cpue <- gapindex::calc_cpue(racebase_tables = ebs_data)
nbs_cpue <- gapindex::calc_cpue(racebase_tables = nbs_data)

ebs_nbs_cpue <- 
  subset(x = rbind(ebs_cpue, nbs_cpue),
         select = c("SURVEY", "YEAR", "STRATUM", "HAULJOIN",
                    "LATITUDE_DD_START", "LATITUDE_DD_END",
                    "LONGITUDE_DD_START", "LONGITUDE_DD_END",
                    "SPECIES_CODE", "WEIGHT_KG", "COUNT",
                    "AREA_SWEPT_KM2", "CPUE_KGKM2", "CPUE_NOKM2"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate Total Abundance 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ebs_popn_stratum <- gapindex::calc_biomass_stratum(racebase_tables = ebs_data,
                                                   cpue = ebs_cpue)
nbs_popn_stratum <- gapindex::calc_biomass_stratum(racebase_tables = nbs_data,
                                                   cpue = nbs_cpue)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate numerical CPUE for a given haul/sex/length bin
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sizecomp <- rbind(gapindex::calc_sizecomp_stratum(
  racebase_tables = ebs_data,
  racebase_cpue = ebs_cpue,
  racebase_stratum_popn = ebs_popn_stratum,
  spatial_level = "haul",
  fill_NA_method = "BS"),
  gapindex::calc_sizecomp_stratum(
    racebase_tables = nbs_data,
    racebase_cpue = nbs_cpue,
    racebase_stratum_popn = nbs_popn_stratum,
    spatial_level = "haul",
    fill_NA_method = "BS"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate globally filled Age-Length Key for EBS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (start_year != start_year_age) { ## i.e., Pcod when age data start 1994
  ebs_data$size <- merge(x = ebs_data$size,
                         y = ebs_data$cruise[, c("CRUISEJOIN", "CRUISE")],
                         by = "CRUISEJOIN")
  ebs_data$size <- subset(x = ebs_data$size, 
                          CRUISE >= start_year_age * 100)
  
  ebs_data$specimen <- subset(x = ebs_data$specimen, 
                              CRUISE >= start_year_age * 100)
  
  nbs_data$size <- merge(x = nbs_data$size,
                         y = nbs_data$cruise[, c("CRUISEJOIN", "CRUISE")],
                         by = "CRUISEJOIN")
  nbs_data$size <- subset(x = nbs_data$size, 
                          CRUISE >= start_year_age * 100)
  
  nbs_data$specimen <- subset(x = nbs_data$specimen, 
                              CRUISE >= start_year_age * 100 & CRUISE != 201801)
}

ebs_alk <- gapindex::calc_alk(racebase_tables = ebs_data, 
                              unsex = "unsex", 
                              global = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate Age-Length Key (ALK) for NBS. First calculate the non-global 
##   ALK. Then fill in the missing ALKs with the globally filled EBS ALK.
##   Then fill in any remaining missing ALKs with the globally filled NBS ALK. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nbs_alk_non_global <- gapindex::calc_alk(racebase_tables = nbs_data, 
                                         unsex = "unsex", 
                                         global = FALSE)
nbs_alk_global <-  gapindex::calc_alk(racebase_tables = nbs_data, 
                                      unsex = "unsex", 
                                      global = TRUE)

## By right merging with the globally filled NBS ALK. The missing ALK values 
## are denoted by NA AGE_FRAC values.
nbs_alk <- 
  merge(x = nbs_alk_non_global,
        y = subset(x = nbs_alk_global, select = -AGE_FRAC),
        all.y = TRUE,
        by = c("SURVEY", "YEAR", "SPECIES_CODE", "SEX", "LENGTH_MM", "AGE"))

## Calculate which ALK records are missing, i.e., the combinations of 
## year/sex/length that don't have probabilities across ages. 
missing_alk <-
  do.call(what = rbind,
          args = lapply(X = split(x = nbs_alk,
                                  f = with(nbs_alk,
                                           list(YEAR,
                                                SPECIES_CODE,
                                                SEX,
                                                LENGTH_MM))),
                        FUN = function(df) {
                          if (sum(df$AGE_FRAC, na.rm = TRUE) == 0)
                            return(data.frame(
                              YEAR =  unique(df$YEAR),
                              SPECIES_CODE = unique(df$SPECIES_CODE),
                              SEX = unique(df$SEX),
                              LENGTH_MM = unique(df$LENGTH_MM)))
                        }
          )
  )

## Remove the records in nbs_alk with the missing alk values
for (irow in 1:nrow(x = missing_alk)) 
  nbs_alk <- 
  subset(x = nbs_alk,
         subset = !(YEAR == missing_alk$YEAR[irow] & 
                      SPECIES_CODE == missing_alk$SPECIES_CODE[irow] & 
                      SEX == missing_alk$SEX[irow] & 
                      LENGTH_MM == missing_alk$LENGTH_MM[irow]) ) 

## Fill in missing alk values using the global EBS alk
nbs_alk_fill <-
  merge(x = missing_alk,
        y = subset(x = ebs_alk, select = -SURVEY),
        all.x = TRUE, suffixes = c("_NG", "_EBS_G"),
        by = c("YEAR", "SPECIES_CODE", "SEX", "LENGTH_MM"))

## Append the EBS-filled alks to nbs_alk
nbs_alk <- rbind(nbs_alk,
                 data.frame(SURVEY = "NBS", nbs_alk_fill))

## Repeat the process of right merging with the global NBS ALK. The missing 
## ALK values are denoted by NA AGE_FRAC values.
nbs_alk <-
  merge(x = nbs_alk,
        y = subset(x = nbs_alk_global, select = -AGE_FRAC),
        all = TRUE,
        by = c("SURVEY", "YEAR", "SPECIES_CODE", "SEX", "LENGTH_MM", "AGE"))

## Find missing ALK values again
missing_alk <-
  do.call(what = rbind,
          args = lapply(X = split(x = nbs_alk,
                                  f = with(nbs_alk,
                                           list(YEAR,
                                                SPECIES_CODE,
                                                SEX,
                                                LENGTH_MM))),
                        FUN = function(df) {
                          if (sum(df$AGE_FRAC, na.rm = TRUE) == 0)
                            return(data.frame(
                              YEAR =  unique(df$YEAR),
                              SPECIES_CODE = unique(df$SPECIES_CODE),
                              SEX = unique(df$SEX),
                              LENGTH_MM = unique(df$LENGTH_MM)))
                        }
          )
  )

## Remove the records in nbs_alk with the missing alk values
for (irow in 1:nrow(x = missing_alk)) 
  nbs_alk <- 
  subset(x = nbs_alk,
         subset = !(YEAR == missing_alk$YEAR[irow] & 
                      SPECIES_CODE == missing_alk$SPECIES_CODE[irow] & 
                      SEX == missing_alk$SEX[irow] & 
                      LENGTH_MM == missing_alk$LENGTH_MM[irow]) ) 

## Fill in missing alk values using the global NBS ALK
nbs_alk_fill <- 
  merge(x = missing_alk,
        y = nbs_alk_global,
        all.x = TRUE,
        by = c("YEAR", "SPECIES_CODE", "SEX", "LENGTH_MM"))

nbs_alk <- rbind(nbs_alk,
                 nbs_alk_fill[, names(nbs_alk)])

## Combine EBS and NBS ALKs
alk <- rbind(ebs_alk, nbs_alk)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Decompose the numerical CPUE (S_ijklm) for a given haul and sex/length 
##   bin across ages. S_ijklm is the numerical CPUE of the ith station in 
##   stratum j for species k, length bin l, and sex m. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ebs_nbs_alk <- merge(x = sizecomp,
                     y = alk,
                     by = c("SURVEY", "YEAR", "SPECIES_CODE", 
                            "SEX", "LENGTH_MM"))
ebs_nbs_alk$AGE_CPUE_NOKM2 <- 
  ebs_nbs_alk$S_ijklm_NOKM2 * ebs_nbs_alk$AGE_FRAC

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Aggregate numerical CPUE across lengths for a given age/sex/haul. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
age_cpue <- rbind(
  ## Aggregate ages younger than the `plus_group`
  stats::aggregate(AGE_CPUE_NOKM2 ~ AGE + HAULJOIN + SPECIES_CODE,
                   data = ebs_nbs_alk,
                   FUN = sum,
                   drop = F,
                   subset = AGE < plus_group),
  ## Aggregate ages at or older than the `plus_group` as one age
  cbind(AGE = plus_group,
        stats::aggregate(AGE_CPUE_NOKM2 ~ HAULJOIN + SPECIES_CODE,
                         data = ebs_nbs_alk,
                         FUN = sum,
                         drop = F,
                         subset = AGE >= plus_group))
)

## Zero-fill CPUEs for missing ages. First we create a grid of all possible
## HAULJOINs and ages. But, there are hauls with positive numerical catches
## but no associated size data. These hauls will be removed.

## Hauls with zero catch
unique_hauls_cpue_zeros <- 
  sort(x = unique(x = ebs_nbs_cpue$HAULJOIN[ebs_nbs_cpue$CPUE_NOKM2 == 0]))
## Hauls with postive count data
unique_hauls_cpue_pos <- 
  sort(x = unique(x = ebs_nbs_cpue$HAULJOIN[ebs_nbs_cpue$CPUE_NOKM2 > 0]))
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

## Append the latitude and longitude information from the `ebs_nbs_cpue` df
## using "HAULJOIN" as a key. 
age_cpue <- merge(x = every_combo_of_ages,
                  y = ebs_nbs_cpue[, c("HAULJOIN", "SURVEY", "YEAR",
                                       "LATITUDE_DD_START",
                                       "LONGITUDE_DD_START")],
                  by = "HAULJOIN")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format VAST Data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
data_geostat_biomass_index <- with(ebs_nbs_cpue,
                                   data.frame(Hauljoin = HAULJOIN,
                                              Region = SURVEY,
                                              Catch_KG = CPUE_KGKM2,
                                              Year = YEAR,
                                              Vessel = "missing", 
                                              AreaSwept_km2 = 1,
                                              Lat = LATITUDE_DD_START,
                                              Lon = LONGITUDE_DD_START,
                                              Pass = 0))

data_geostat_numerical_index <- with(ebs_nbs_cpue,
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
dir_out <- 
  paste0("species_specific_code/BS/", species_name, "/", which_model, "/data/")
if (!dir.exists(paths = dir_out)) dir.create(path = dir_out, recursive = T)
for (ifile in c("data_geostat_biomass_index", 
                "data_geostat_numerical_index",
                "data_geostat_agecomps", 
                "sizecomp", "alk", "ebs_data", "nbs_data")) 
  saveRDS(object = get(x = ifile), 
          file = paste0(dir_out, ifile, ".RDS"))

