##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Prepare model-based product data inputs from the gapindex R package
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Libraries, using gapindex version 3.0.2
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# devtools::install_github(repo = "afsc-gap-products/gapindex@v3.0.2", 
#                          dependencies = TRUE)
library(gapindex)
library(sf)
library(sdmTMB)
library(here)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Connect to Oracle. Make sure you're on the VPN/network. 
##   Set up constants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
channel <- gapindex::get_connected(check_access = F)

## specify whether hindcast or production phase, create data folder
phase <- c("hindcast", "production")[2]
data_dir <- paste0("data/GOA/", phase, "/")
if (!dir.exists(paths = data_dir)) dir.create(path = data_dir, recursive = TRUE)

year_start <- 1990
year_end <- 2025

species_df <- data.frame(SPECIES_CODE = c(310, 10110, 21720, 21740, 30060, 
                                          30420, 30152),
                         START_YEAR = c(rep(year_start, 6), 1996))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull catch and effort data and calculate CPUE (zero-filled)
##   Standard data: ABUNDANCE_HAUL = "Y" hauls
##   
##   Non-standard data: hauls in NMFS areas 519 (Unimak Pass) and 659 (SE
##   Inside) that were previously ABUNDANCE_HAUL = "Y" but are now outside the
##   survey footprint when the survey design was updated in 2025. These hauls
##   are now coded as HAUL_TYPE 24
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_standard_data <- 
  gapindex::get_data(survey_set = "GOA",
                     year_set = year_start:year_end,
                     spp_codes = species_df$SPECIES_CODE,
                     pull_lengths = FALSE,
                     taxonomic_source = "GAP_PRODUCTS.TAXONOMIC_CLASSIFICATION",
                     channel = channel)
goa_standard_cpue <- 
  as.data.frame(x = gapindex::calc_cpue(gapdata = goa_standard_data))

goa_nonstandard_data <-
  gapindex::get_data(survey_set = "GOA",
                     year_set = year_start:year_end,
                     spp_codes = species_df$SPECIES_CODE,
                     abundance_haul = "N", haul_type = 24,
                     pull_lengths = FALSE,
                     taxonomic_source = "GAP_PRODUCTS.TAXONOMIC_CLASSIFICATION",
                     channel = channel)
goa_nonstandard_cpue <-
  as.data.frame(x = gapindex::calc_cpue(gapdata = goa_nonstandard_data))

## Combine standard and non-standard CPUE
goa_all_cpue <- rbind(goa_standard_cpue, goa_nonstandard_cpue)

## Merge the species name from goa_standard_data$species to goa_standard_cpue using
## the SPECIES_CODE field. goa_standard_data$species contains a complex of the 
## Sebastes sp. (dusky/dark) and dusky RF species codes, so the workaround is 
## only subsetting for species-level records in gapindex$species. Note: we
## may need a more elegant solution if we include more complexes... 
goa_all_cpue <- merge(x = goa_all_cpue,
                      y = subset(x = goa_standard_data$species,
                                 subset = ID_RANK == "species",
                                 select = c(SPECIES_CODE, SPECIES_NAME)),
                      by = "SPECIES_CODE")

## For species with later starting years (i.e., dusky rockfish), remove
## records prior to the starting year for those species
for (ispp in 1:nrow(x = species_df)) 
  goa_all_cpue <-  
  subset(x = goa_all_cpue, 
         subset = !(SPECIES_CODE == species_df$SPECIES_CODE[ispp] &
                      YEAR < species_df$START_YEAR[ispp]))

table(goa_all_cpue$YEAR, goa_all_cpue$SPECIES_CODE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create a lat/lon spatial object of the station locations 
##   Transform station spatial object to UTM (zone 5)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cpue_sf <- sf::st_as_sf(x = goa_all_cpue, 
                        coords = c("LONGITUDE_DD_START", "LATITUDE_DD_START"), 
                        crs = "+proj=longlat +datum=WGS84")
goa_all_cpue[, c("E_km_z5", "N_km_z5")] <-
  sf::st_coordinates(sf::st_transform(x = cpue_sf, 
                                      crs = "+proj=utm +zone=5 +units=km")) 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format CPUE data for sdmTMB() and save
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_allspp <- with(goa_all_cpue, 
                   data.frame(hauljoin = HAULJOIN,
                              year = as.integer(YEAR),
                              lon = LONGITUDE_DD_START,
                              lat = LATITUDE_DD_START,
                              X = E_km_z5,
                              Y = N_km_z5,
                              depth_m = DEPTH_M,
                              bottom_temp_c = BOTTOM_TEMPERATURE_C,
                              species_code = SPECIES_CODE,
                              species = SPECIES_NAME,
                              catch_kg = WEIGHT_KG, 
                              catch_n = COUNT,
                              effort_km2 = AREA_SWEPT_KM2,
                              cpue_kg_km2 = CPUE_KGKM2,
                              cpue_n_km2 = CPUE_NOKM2))

## Save catch and effort data
saveRDS(dat_allspp, file = paste0(data_dir, "dat_allspp.RDS"))
