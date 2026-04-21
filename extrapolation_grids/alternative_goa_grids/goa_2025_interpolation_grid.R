##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       GOA interpolation grid under 2025 GOA survey footprint
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Import libraries
library(akgfmaps)
library(terra)
library(sf)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create 2-nmi grid within the 2025 GOA footprint
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_footprint_2025 <- 
  ## Import goa strata
  akgfmaps::get_base_layers(
    select.region = "goa",
    design.year = 2025,
    set.crs = "+proj=utm +zone=5 +units=km"
  )$survey.strata |>
  ## Remove strata in the 700 - 1000 m depth zone
  subset(subset = STRATUM < 500) |>
  ## Dissolve boundaries
  terra::vect() |>
  terra::aggregate() 

goa_grid <- 
  ## Create a rectangular grid that overlaps with the 2025 GOA footprint
  sf::st_make_grid(x = goa_footprint_2025, 
                   square = T, 
                   cellsize = 3.704) |> 
  terra::vect() |> 
  ## Intersect with the footprint
  terra::intersect(goa_footprint_2025)

goa_grid_df <- data.frame(terra::centroids(x = goa_grid) |> terra::crds()) 
names(x = goa_grid_df) <- c("X", "Y")

goa_grid_df[, c("lon", "lat")] <- 
  terra::project(x = goa_grid, "+proj=longlat +datum=WGS84") |> 
  terra::centroids() |> 
  terra::crds()

goa_grid_df$area_km2 <- terra::expanse(x = goa_grid) / 1e6

goa_grid_df <- subset(x = goa_grid_df,
                      subset = area_km2 >= 0.001)

## Check to make sure the total area of the interpolation grid matches the 
## actual survey footprint.
sum(goa_grid_df$area_km2) # 301174.5 km2
expanse(goa_footprint_2025) * 1e-6 # 301174.7 km2

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Write to file
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(x = goa_grid_df,
          file = paste0("extrapolation_grids/",
                        "goa_2025_interpolation_grid.csv"),
          row.names = F)
