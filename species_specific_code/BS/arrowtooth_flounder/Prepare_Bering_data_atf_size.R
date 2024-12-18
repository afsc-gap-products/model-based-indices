##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## EBS/NBS data pull for ModSquad VAST input via gapindex R package
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

species_info <- data.frame(species_name = "arrowtooth_flounder",
                           species_code = 10110,
                           start_year = 1982,
                           current_year = 2024)

## Set constants
start_year <- species_info$start_year
current_year <- species_info$current_year
species_code <- species_info$species_code
species_name <- species_info$species_name

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
                                        channel = channel,
                                        remove_na_strata = TRUE)

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
                                        channel = channel)

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
               channel = channel,
               query = "SELECT CRUISEJOIN, CRUISE
                          FROM RACEBASE.CRUISE 
                          WHERE CRUISE IN (198502, 198808, 199102) 
                          AND REGION = 'BS'"),
             VESSEL_ID = NA,
             VESSEL_NAME = NA,
             DESIGN_YEAR = 2022)

nbs_other$survey <- nbs_other$cruise[, names(x = nbs_standard_data$survey)]

nbs_other$haul <- 
  RODBC::sqlQuery(channel = channel,
                  query = "SELECT *
                             FROM RACEBASE.HAUL 
                             WHERE CRUISE IN (198502, 198808, 199102) 
                             AND HAUL_TYPE = 3 
                             AND PERFORMANCE >= 0 
                             AND ABUNDANCE_HAUL = 'Y'")
nbs_other$catch <- 
  RODBC::sqlQuery(
    channel = channel,
    query = paste("SELECT HAULJOIN, SPECIES_CODE, WEIGHT, NUMBER_FISH
                     FROM RACEBASE.CATCH 
                     WHERE SPECIES_CODE IN", 
                  gapindex::stitch_entries(species_code),
                  "AND HAULJOIN IN",
                  gapindex::stitch_entries(nbs_other$haul$HAULJOIN)))

nbs_other$size <-
  RODBC::sqlQuery(
    channel = channel,
    query = paste("SELECT CRUISEJOIN, HAULJOIN, SPECIES_CODE, LENGTH, 
                     SEX, FREQUENCY 
                     FROM RACEBASE.LENGTH 
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
                                 channel = channel)
nbs18_data$cruise$SURVEY <- nbs18_data$survey$SURVEY <- "NBS"
nbs18_data$cruise$SURVEY_DEFINITION_ID <- 
  nbs18_data$survey$SURVEY_DEFINITION_ID  <- 143

## Combine all the NBS data into one list. 
nbs_data <- 
  list(survey = rbind(nbs_standard_data$survey, 
                      nbs_other$survey, 
                      nbs18_data$survey),
       survey_design = nbs_standard_data$survey_design,
       cruise = rbind(nbs_standard_data$cruise[, 
                                               names(x = nbs18_data$cruise), 
                                               with = F],
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

nbs_data$catch$HAULJOIN <- as.integer(nbs_data$catch$HAULJOIN)
nbs_data$catch$NUMBER_FISH <- as.numeric(nbs_data$catch$NUMBER_FISH)
nbs_data$catch$WEIGHT <- as.numeric(nbs_data$catch$WEIGHT)
nbs_data$size$HAULJOIN <- as.integer(nbs_data$size$HAULJOIN)
nbs_data$size$LENGTH <- as.integer(nbs_data$size$LENGTH)
nbs_data$size$FREQUENCY <- as.integer(nbs_data$size$FREQUENCY)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate CPUE for EBS and NBS, rbind, and reorder columns
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ebs_cpue <- gapindex::calc_cpue(gapdata = ebs_data)
nbs_cpue <- gapindex::calc_cpue(gapdata = nbs_data)

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
ebs_popn_stratum <- gapindex::calc_biomass_stratum(gapdata = ebs_data,
                                                   cpue = ebs_cpue)
nbs_popn_stratum <- gapindex::calc_biomass_stratum(gapdata = nbs_data,
                                                   cpue = nbs_cpue)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate numerical CPUE (S_ijklm) for a given haul and sex/length 
##   bin across ages. S_ijklm is the numerical CPUE of the ith station in 
##   stratum j for species k, length bin l, and sex m. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sizecomp <- rbind(gapindex::calc_sizecomp_stratum(
  gapdata = ebs_data,
  cpue = ebs_cpue, 
  abundance_stratum = ebs_popn_stratum,
  spatial_level = "haul",
  fill_NA_method = "BS"),
  gapindex::calc_sizecomp_stratum(
    gapdata = nbs_data,
    cpue = nbs_cpue,
    abundance_stratum = nbs_popn_stratum,
    spatial_level = "haul",
    fill_NA_method = "BS"))

dat <- dplyr::left_join(sizecomp, 
                        dplyr::select(ebs_nbs_cpue, 
                               HAULJOIN, 
                               LATITUDE_DD_START,
                               LONGITUDE_DD_START,
                               AREA_SWEPT_KM2), 
                        by = "HAULJOIN")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save output
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir_out <- paste0("species_specific_code/BS/", species_name, "/", phase, "/data/")
if (!dir.exists(paths = dir_out)) dir.create(path = dir_out, recursive = T)
for (ifile in c("dat", "sizecomp", "ebs_data", "nbs_data")) 
  saveRDS(object = get(x = ifile), 
          file = paste0(dir_out, ifile, ".RDS"))