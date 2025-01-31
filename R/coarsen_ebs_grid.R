## Title: coarsen Bering prediction grids for index models ----
# date: Jan 31 2025
# author: Lewis Barnett, with code from Sean Rohan

library(akgfmaps)

# Setup grid using dimensions used to create 2022 updated grid
bbox = c(xmin = -2846314.206900,
        ymin = 251495.775400,
        xmax = 2128157.793100,
        ymax = 2859111.775400)

# EBS alone ----
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")

# original resolution is 2x2nm (3704, 3704), coarsened to 12x12nm here
ebs_grid <- akgfmaps:::make_2d_grid(obj = dplyr::filter(ebs_layers$survey.strata, 
                                                        SURVEY_DEFINITION_ID == 98),
                                    resolution = c(3704, 3704)*6, 
                                    bbox = bbox,
                                    output_type = "point",
                                    include_tile_center = TRUE) |>
  sf::st_transform(crs = "EPSG:32602")

ebs_grid[, c('LON_UTM', "LAT_UTM")] <- sf::st_coordinates(ebs_grid)

# change units to km 
ebs_grid <- ebs_grid |>
  dplyr::mutate(X = LON_UTM / 1000,
                Y = LAT_UTM / 1000,
                area_km2 = as.numeric(AREA)/1e6) |> 
  dplyr::select(X, Y, area_km2) |>
  as.data.frame()

ggplot2::ggplot(ebs_grid, aes(X, Y, colour = area_km2)) +
  geom_tile(width = 2, height = 2, fill = NA) +
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.5) +
  coord_fixed()

write.csv(ebs_grid, file = here::here("extrapolation_grids", "ebs_coarse_grid.csv"))

# NBS alone ----
nbs_grid <- akgfmaps:::make_2d_grid(obj = dplyr::filter(ebs_layers$survey.strata, 
                                                        SURVEY_DEFINITION_ID == 143),
                                    resolution = c(3704, 3704)*6, 
                                    bbox = bbox,
                                    output_type = "point",
                                    include_tile_center = TRUE) |>
  sf::st_transform(crs = "EPSG:32602")

nbs_grid[, c('LON_UTM', "LAT_UTM")] <- sf::st_coordinates(nbs_grid)

# change units to km 
nbs_grid <- nbs_grid |>
  dplyr::mutate(X = LON_UTM / 1000,
                Y = LAT_UTM / 1000,
                area_km2 = as.numeric(AREA)/1e6) |> 
  dplyr::select(X, Y, area_km2) |>
  as.data.frame()

ggplot2::ggplot(nbs_grid, aes(X, Y, colour = area_km2)) +
  geom_tile(width = 2, height = 2, fill = NA) +
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.5) +
  coord_fixed()

write.csv(nbs_grid, file = here::here("extrapolation_grids", "nbs_coarse_grid.csv"))

# combine to get full Bering grid ----
grid <- dplyr::bind_rows(nbs_grid, ebs_grid)

ggplot2::ggplot(grid, aes(X, Y, colour = area_km2)) +
  geom_tile(width = 2, height = 2, fill = NA) +
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.5) +
  coord_fixed()

nrow(grid) # check N cells: results in 1908 for 12x12nm resolution

write.csv(grid, file = here::here("extrapolation_grids", "bering_coarse_grid.csv"))