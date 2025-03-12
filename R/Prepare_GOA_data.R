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
phase <- c("hindcast", "production")[1]
data_dir <- paste0("data/GOA/", phase, "/")
if (!dir.exists(paths = data_dir)) dir.create(path = data_dir, recursive = TRUE)

year_start <- 1990
year_end <- 2023

species_df <- data.frame(SPECIES_CODE = c(310, 10110, 21720, 21740, 30060, 
                                          30420, 30152),
                         START_YEAR = c(rep(year_start, 6), 1996))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull catch and effort data and calculate CPUE (zero-filled)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gapindex_data <- 
  gapindex::get_data(survey_set = "GOA",
                     year_set = year_start:year_end,
                     spp_codes = species_df$SPECIES_CODE,
                     pull_lengths = FALSE,
                     taxonomic_source = "GAP_PRODUCTS.TAXONOMIC_CLASSIFICATION",
                     channel = channel)

gapindex_cpue <- as.data.frame(gapindex::calc_cpue(gapdata = gapindex_data))

## Merge the species name from gapindex_data$species to gapindex_cpue using
## the SPECIES_CODE field. gapindex_data$species contains a complex of the 
## Sebastes sp. (dusky/dark) and dusky RF species codes, so the workaround is 
## only subsetting for species-level records in gapindex$species. Note: we
## may need a more elegant solution if we include more complexes... 
gapindex_cpue <- merge(x = gapindex_cpue,
                       y = subset(x = gapindex_data$species,
                                  subset = ID_RANK == "species",
                                  select = c(SPECIES_CODE, SPECIES_NAME)),
                       by = "SPECIES_CODE")

## For species with later starting years (i.e., dusky rockfish), remove
## records prior to the starting year for those species
for (ispp in 1:nrow(x = species_df)) 
  gapindex_cpue <-  
  subset(x = gapindex_cpue, 
         subset = !(SPECIES_CODE == species_df$SPECIES_CODE[ispp] &
                      YEAR < species_df$START_YEAR[ispp]))

table(gapindex_cpue$YEAR, gapindex_cpue$SPECIES_CODE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create a lat/lon spatial object of the station locations 
##   Transform station spatial object to UTM (zone 5)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cpue_sf <- sf::st_as_sf(x = gapindex_cpue, 
                        coords = c("LONGITUDE_DD_START", "LATITUDE_DD_START"), 
                        crs = "+proj=longlat +datum=WGS84")
gapindex_cpue[, c("E_km_z5", "N_km_z5")] <-
  sf::st_coordinates(sf::st_transform(x = cpue_sf, 
                                      crs = "+proj=utm +zone=5 +units=km")) 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format CPUE data for sdmTMB() and save
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_allspp <- with(gapindex_cpue, data.frame(hauljoin = HAULJOIN,
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