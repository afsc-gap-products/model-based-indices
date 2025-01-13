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
  terra::vect(x = paste0("extrapolation_grids/alternative_goa_grids/",
                         "goa_strata_2025.gpkg"))|> 
  ## Project to UTM 
  terra::project("+proj=utm +zone=5 +units=km") |>
  ## Remove strata in the 700 - 1000 m depth zone
  subset(DEPTH_MAX_M != 1000, NSE = TRUE) |>
  ## Dissolve boundaries
  terra::aggregate()

goa_grid <- 
  ## Create a rectangular grid that overlaps with the 2025 GOA footprint
  sf::st_make_grid(x = goa_footprint_2025, 
                             square = T, 
                             cellsize = 3.704) |> 
  terra::vect() |> 
  ## Intersect with the footprint
  terra::intersect(goa_footprint_2025)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Compare with prior interpolation grid
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_base <- akgfmaps::get_base_layers(select.region = "goa", 
                                      set.crs = "+proj=utm +zone=5 +units=km")
## Import prior interpolation grid
prior_grid <- 
  read.csv(file = "extrapolation_grids/GOAThorsonGrid_Less700m.csv") |>
  terra::vect(geom = c("Longitude", "Latitude"), 
              crs = "+proj=longlat +datum=WGS84") |>
  terra::project("+proj=utm +zone=5 +units=km")

## Compare total survey areas
sum(terra::expanse(x = goa_grid) / 1e6) ## 300459.4 km2
sum(prior_grid$Shape_Area / 1e6) ## 304953.1 km2
# 1.473571% reduction in total area

## Plot differences
plot(goa_footprint_2025, mar = c(0,0,0,0))
plot(prior_grid, add = TRUE, cex = 0.2, pch = 15, col = "red")
plot(goa_grid, add = TRUE, col = "black", border = TRUE)
plot(st_geometry(goa_base$akland), add = TRUE, col = "tan")
