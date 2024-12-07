##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare GOA VAST input to output from the gapindex R package
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Libraries, using gapindex version 3.0.2
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
devtools::install_github(repo = "afsc-gap-products/gapindex@v3.0.2", dependencies = TRUE)
library(gapindex)
library(sf)
library(sdmTMB)
library(here)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Connect to Oracle. Make sure you're on the VPN/network. 
##   Set up constants
##   Double-check how dusky rockfish is pulled...
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
channel <- gapindex::get_connected(check_access = F)

species_df <- data.frame(
  SPECIES_CODE = c(310, 10110, 21720, 21740, 30060, 30420, 30150, 30152),
  GROUP_CODE =   c(310, 10110, 21720, 21740, 30060, 30420, 30152, 30152)
)

year_start <- 1990
year_end <- 2023

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull catch and effort data and calculate CPUE (zero-filled)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gapindex_data <- 
  gapindex::get_data(survey_set = "GOA",
                     year_set = year_start:year_end,
                     spp_codes = species_df,
                     pull_lengths = FALSE,
                     channel = channel)

gapindex_cpue <- as.data.frame(gapindex::calc_cpue(gapdata = gapindex_data))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Transform lat/lons to UTM (zone 5)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cpue_sf <- sf::st_as_sf(x = gapindex_cpue, 
                        coords = c("LONGITUDE_DD_START", "LATITUDE_DD_START"), 
                        crs = "+proj=longlat +datum=WGS84")
gapindex_cpue[, c("E_km_z5", "N_km_z5")] <-
  sf::st_coordinates(sf::st_transform(x = cpue_sf, 
                                      crs = "+proj=utm +zone=5 +units=km")) 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format CPUE data for sdmTMB()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_allspp <- with(gapindex_cpue, data.frame(hauljoin = HAULJOIN,
                                             species_code = SPECIES_CODE,
                                             cpue_kg_km2 = CPUE_KGKM2, 
                                             year = as.integer(YEAR),
                                             X = E_km_z5,
                                             Y = N_km_z5))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format mesh
##   I downloaded the mesh from a 2023 GOA PCod production VAST fit, not 
##   sure how we standardize this bit...
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# goa_cod_vast_fit <- readRDS(file = "C:/Users/zack.oyafuso/Downloads/Gadus_macrocephalus_VASTfit.RDS")
# goa_vast_mesh <- goa_cod_vast_fit$spatial_list$MeshList$anisotropic_mesh
# saveRDS(object = goa_vast_mesh, file = "extrapolation_grids/goa_vast_mesh.RDS")
goa_vast_mesh <- readRDS(file = "extrapolation_grids/goa_vast_mesh.RDS")

goa_sdmtmb_mesh <- sdmTMB::make_mesh(data = dat_allspp, 
                                     xy_cols = c("X", "Y"), 
                                     mesh = goa_vast_mesh) #pass same mesh as VAST model

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format prediction grid
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_grid <- read.csv(here("extrapolation_grids", "GOAThorsonGrid_Less700m.csv"))
goa_grid <- data.frame(Lat=goa_grid$Latitude, 
                       Lon=goa_grid$Longitude,
                       Area_km2=goa_grid$Shape_Area/1000000)

sf_grid <- sf::st_as_sf(x = goa_grid,
                        coords = c("Lon", "Lat"),
                        crs = "+proj=longlat +datum=WGS84"
)
sf_grid <- sf::st_transform(sf_grid, crs = "+proj=utm +zone=5 +units=km")
goa_grid[, c("X", "Y")] <- sf::st_coordinates(x = sf_grid)

goa_sdmtmb_grid <- goa_grid[, c("Lon", "Lat", "X", "Y", "Area_km2")]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save data objects somewhere...
##   Things to consider: where do these data objects live? Which, if not all,
##   get pushed to GitHub? Should we save dat_allspp itself or loop through
##   and subset and save individual species cpue records? 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
