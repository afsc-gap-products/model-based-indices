# coarsen Bering prediction grid for faster index model computation

library(akgfmaps)

# Load survey area polygon
bs <- akgfmaps::get_base_layers("ebs", set.crs = "EPSG:32602")

ggplot() +
  geom_sf(data = bs$survey.area,
          fill = NA,
          mapping = aes(color = "Survey area"))

# create a coarsened grid using the raster package
resolution_original <- 3704
resolution <- resolution_original * 6 # define what factor to coarsen by

r <- raster::raster(as(bs$survey.area, "Spatial"), resolution = resolution)
rr <- raster::rasterize(as(bs$survey.area, "Spatial"), r, getCover = TRUE)

grid <- as.data.frame(raster::rasterToPoints(rr))
grid$area <- grid$layer * resolution * resolution
grid <- dplyr::filter(grid, area > 0) |>
  dplyr::mutate(X = x / 1000, Y = y / 1000) 
grid$area_km2 <- grid$area/1e6
grid <- dplyr::select(grid, X, Y, area_km2)

ggplot(grid, aes(X, Y, colour = area_km2)) +
  geom_tile(width = 2, height = 2, fill = NA) +
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.5) +
  coord_fixed()

nrow(grid) # check N cells: results in 1539 for 6x resolution

write.csv(grid, file = here::here("extrapolation_grids", "bering_coarse_grid.csv"))

# repeat above for EBS alone ----

ebs <- akgfmaps::get_base_layers("sebs", set.crs = "EPSG:32602")

r_ebs <- raster::raster(as(ebs$survey.area, "Spatial"), resolution = resolution)
rr_ebs <- raster::rasterize(as(ebs$survey.area, "Spatial"), r_ebs, getCover = TRUE)

grid_ebs <- as.data.frame(raster::rasterToPoints(rr_ebs))
grid_ebs$area <- grid_ebs$layer * resolution * resolution
grid_ebs <- dplyr::filter(grid_ebs, area > 0) |>
  dplyr::mutate(X = x / 1000, Y = y / 1000) 
grid_ebs$area_km2 <- grid_ebs$area/1e6
grid_ebs <- dplyr::select(grid_ebs, X, Y, area_km2)

ggplot(grid_ebs, aes(X, Y, colour = area_km2)) +
  geom_tile(width = 2, height = 2, fill = NA) +
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.5) +
  coord_fixed()

write.csv(grid_ebs, file = here::here("extrapolation_grids", "ebs_coarse_grid.csv"))

# repeat above for NBS alone ----

nbs <- akgfmaps::get_base_layers("nbs", set.crs = "EPSG:32602")

r_nbs <- raster::raster(as(nbs$survey.area, "Spatial"), resolution = resolution)
rr_nbs <- raster::rasterize(as(nbs$survey.area, "Spatial"), r_nbs, getCover = TRUE)

grid_nbs <- as.data.frame(raster::rasterToPoints(rr_nbs))
grid_nbs$area <- grid_nbs$layer * resolution * resolution
grid_nbs <- dplyr::filter(grid_nbs, area > 0) |>
  dplyr::mutate(X = x / 1000, Y = y / 1000) 
grid_nbs$area_km2 <- grid_nbs$area/1e6
grid_nbs <- dplyr::select(grid_nbs, X, Y, area_km2)

ggplot(grid_nbs, aes(X, Y, colour = area_km2)) +
  geom_tile(width = 2, height = 2, fill = NA) +
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.5) +
  coord_fixed()

write.csv(grid_nbs, file = here::here("extrapolation_grids", "nbs_coarse_grid.csv"))