#' Code in development for pollock abundance within 100km of the Pribilof 
#' Islands for research on fur seal foraging success. Code managed and updated 
#' by Sophia Wassermann, based on sdmTMB code adapted from bs_indices.Rmd 
#' (developed by Lewis Barnett)

library(sdmTMB)
library(dplyr)
library(ggplot2)
library(here)
library(devtools)
library(sf)

# install_github("afsc-gap-products/akgfmaps", build_vignettes = TRUE)
library(akgfmaps)

# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
theme_set(theme_sleek())

# Get pollock CPUE data -------------------------------------------------------
phase <- "hindcast" # determines (combined w/ the year) which cycle the data are from

this_year <- as.numeric(format(Sys.Date(), "%Y"))
if(phase == "hindcast") {this_year <- this_year - 1}  

results_wd <- here("species_specific_code", "BS", "pollock", "research", "250kts")

dat <- read.csv(here("species_specific_code", "BS", "pollock", phase,
                     "data", paste0("VAST_ddc_all_", this_year, ".csv")))  
dat <- transmute(dat,
                 cpue = ddc_cpue_kg_ha * 100, # converts cpue from kg/ha to kg/km^2
                 year = as.integer(year),
                 lat = start_latitude,
                 lon = start_longitude)


# detect if any years have occurrences at every haul and fix params as needed
mins <- dat %>% group_by(year) %>% summarize(min = min(cpue))
if(sum(mins$min) == 0){
  control = sdmTMBcontrol()
} else {
  no_zero_yr <- as.integer(mins %>% filter(min > 0) %>% select(year))
  # set up map and fix value of p(occurrence) to slightly less than 1:
  yrs <- sort(unique(factor(dat$year)))
  .map <- seq_along(yrs)
  .map[yrs %in% no_zero_yr] <- NA
  .map <- factor(.map)
  .start <- rep(0, length(yrs))
  .start[yrs %in% no_zero_yr] <- 20
  
  control =  sdmTMBcontrol(
    map = list(b_j = .map),
    start = list(b_j = .start)
  )
}

# Set up cold pool covariate
# devtools::install_github("afsc-gap-products/coldpool")
env_df <- coldpool::cold_pool_index

env <- cbind(env_df, env = scale(coldpool::cold_pool_index$AREA_LTE2_KM2)) %>%
  mutate(year = as.integer(YEAR)) %>%
  select(year, env)

dat <- left_join(dat, env, by = "year") 

# Final data manipulation steps
dat$year_f <- as.factor(dat$year)

dat <- add_utm_columns(dat, ll_names = c("lon", "lat"), utm_crs = 32602, units = "km")

# Fit model (if needed) -------------------------------------------------------
f1 <- here("species_specific_code", "BS", "pollock", phase, "results", "fit_250_knots.RDS")
if (file.exists(f1)) {
  fit <- readRDS(f1)
  } else {
    mesh <-  make_mesh(dat, xy_cols = c("X", "Y"), 
                       mesh = fmesher::fm_as_fm(readRDS(file = here("meshes", "bs_vast_mesh_250_knots.RDS"))))
    fit <- sdmTMB( 
      cpue ~ 0 + year_f,
      spatial_varying = ~ env,
      data = dat, 
      mesh = mesh,
      family = delta_gamma(type = "poisson-link"), 
      time = "year", 
      spatial = "on",
      spatiotemporal = "ar1",
      extra_time = 2020L, 
      silent = FALSE,
      anisotropy = TRUE,
      control = control
    )
  }

# Check fit
sanity(fit)
summary(fit)

# Make predictions and index --------------------------------------------------
load(here("species_specific_code", "BS", "pollock", "research", "Pribilof_polygons.RData"))
# Projected crs = 3338, unprojected crs = 4326
final_combined_hr_polygons_projected_sf$associated_circle_radius_meters[1] <- "Pribilofs"

index_by_area <- function(area) {
  polygon <- final_combined_hr_polygons_projected_sf$geometry[area]
  dir_name <- paste0("radius", final_combined_hr_polygons_projected_sf$associated_circle_radius_meters[area])
  dir.create(here(results_wd, dir_name), recursive = TRUE, showWarnings = FALSE)
  
  grid <- make_2d_grid(obj = polygon,
                       resolution = c(3704, 3704),  # default resolution - 2x2nm
                       output_type = "point",
                       include_tile_center = TRUE) %>%
    st_transform(crs = "EPSG:32602") 
  
  grid[, c('LON_UTM', "LAT_UTM")] <- st_coordinates(grid)
  
  grid <- data.frame(grid) %>%
    select(LON_UTM, LAT_UTM, AREA) %>%
    mutate(X = LON_UTM / 1000,
           Y = LAT_UTM / 1000,
           area_km2 = as.numeric(AREA)/1e6) %>%
    select(X, Y, area_km2)
  grid <- as.data.frame(as.matrix(grid)) # drop attributes
  
  ggplot(grid, aes(X, Y, colour = area_km2)) +
    geom_tile(width = 2, height = 2, fill = NA) +
    scale_colour_viridis_c(direction = -1) +
    geom_point(size = 0.5) +
    coord_fixed()
  ggsave(file = here(results_wd, dir_name, "pred_grid.pdf"),
         height = 4, width = 4.5, units = c("in"))
  
  # replicate prediction grid for each year in data
  pred_grid <- replicate_df(grid, "year_f", unique(dat$year_f))
  pred_grid$year <- as.integer(as.character(factor(pred_grid$year_f)))
  
  # join with environmental covariate (cold pool)
  pred_grid <- left_join(pred_grid, env, by = "year")
  
  # get prediction
  p <- predict(fit, newdata = pred_grid, return_tmb_object = TRUE)
  save(p, file = here(results_wd, dir_name, "pred.Rdata"))
  
  # get index
  ind <- get_index(p, bias_correct = TRUE, area = pred_grid$area_km2)
  ind$stratum <- paste0("radius = ", final_combined_hr_polygons_projected_sf$associated_circle_radius_meters[area])
  write.csv(ind, file = here(results_wd, dir_name, "index.csv"), row.names = FALSE)
  
  return(ind)
}

# Produce index for each of the areas
indices <- data.frame()
for(i in 1:nrow(final_combined_hr_polygons_projected_sf)) {
  rbind.data.frame(indices, index_by_area(i))
}

# Plot index, scaled from kg to Mt
ggplot(indicies, aes(x = year, y = (est / 1e9))) +
  geom_point() +
  geom_line() +
  ylim(0, NA) +
  geom_ribbon(aes(ymin = (lwr / 1e9), ymax = (upr / 1e9)), alpha = 0.4) +
  xlab("") + ylab("Biomass (Mt)") +
  facet_wrap(~stratum)
ggsave(file = here(results_wd, "index.pdf"), 
       height = 3, width = 5, units = "in")

# Plot predicted density maps and fit diagnostics -----------------------------
# q-q plot
pdf(file = here(results_wd, "qq.pdf"),
    width = 5, height = 5)
sims <- simulate(fit, nsim = 500, type = "mle-mvn") 
sims |> dharma_residuals(fit, test_uniformity = FALSE)
dev.off()

#residuals on map plot, by year
resids <- sims |>
  dharma_residuals(fit, test_uniformity = FALSE, return_DHARMa = TRUE)
fit$data$resids <- resids$scaledResiduals

ggplot(subset(fit$data, !is.na(resids) & is.finite(resids)), aes(X, Y, col = resids)) +
  scale_colour_gradient2(name = "residuals", midpoint = 0.5) +
  geom_point(size = 0.7) +
  scale_x_continuous(breaks = c(250, 750)) +
  scale_y_continuous(breaks = c(6000, 6500, 7000)) +
  facet_wrap(~year) +
  coord_fixed() 
ggsave(file = here(results_wd, "residuals_map.pdf"),
       height = 9, width = 6.5, units = "in")

# predictions on map plot, by year
ggplot(p$data, aes(X, Y, fill = exp(est1 + est2))) +
  geom_tile(width = 10, height = 10) +
  scale_fill_viridis_c(trans = "sqrt", name = "") +
  scale_x_continuous(breaks = c(250, 750)) +
  scale_y_continuous(breaks = c(6000, 6500, 7000)) +
  facet_wrap(~year) +
  coord_fixed() +
  ggtitle("Predicted densitites (kg / square km)") 
ggsave(file = here(results_wd, "predictions_map.pdf"),
       height = 9, width = 6.5, units = "in")
