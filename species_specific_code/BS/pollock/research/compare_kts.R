#' Script for comparing indices using different numbers of knots across the 
#' Pribilof polygons. Theoretically, increasing the number of knots should be
#' more important for the smaller areas.

library(here)
library(dplyr)
library(ggplot2)
library(ggsidekick)
theme_set(theme_sleek())

# Set up ----------------------------------------------------------------------
wd <- here("species_specific_code", "BS", "pollock", "research")

# Read in polygons
load(here(wd, "Pribilof_polygons.RData"))
# Clean up labelling for later plotting
final_combined_hr_polygons_projected_sf$associated_circle_radius_meters <- final_combined_hr_polygons_projected_sf$associated_circle_radius_meters / 1000
final_combined_hr_polygons_projected_sf$associated_circle_radius_meters[1] <- "Pribilofs"

# Read in data ----------------------------------------------------------------
read_indices <- function(kts) {
  ind_list <- list()
  for(i in 1:nrow(final_combined_hr_polygons_projected_sf)) {
    dir_name <- paste0("radius", final_combined_hr_polygons_projected_sf$associated_circle_radius_meters[i])
    df <- read.csv(here(wd, paste0(kts, "kts"), dir_name, "index.csv"))
    ind_list[[i]] <- df
  }
  
  ind <- do.call(rbind, ind_list)
  if(!file.exists(here(wd, paste0(kts, "kts"), "index.csv"))) {
    write.csv(ind, here(wd, paste0(kts, "kts"), "index.csv"), row.names = FALSE)
  }
  ind$knots <- kts
  return(ind)
}

indices <- rbind.data.frame(read_indices(250), read_indices(500))

# Final tweaks & plot ---------------------------------------------------------
indices$stratum  <- factor(indices$stratum,
                           levels = c("25", "50", "75", "100", "125", "150", "175", "200", "225", "250", "Pribilofs"),
                           labels = c("25km", "50km", "75km", "100km", "125km", "150km", "175km", "200km", "225km", "250km", "Pribilofs"))
indices$knots <- factor(indices$knots)

# Plot indices (converted to millions of tons)
ggplot(indices, aes(x = year, y = (est / 1e9)), color = knots) +  
  geom_line() +
  ylim(0, NA) +
  geom_ribbon(aes(ymin = (lwr / 1e9), ymax = (upr / 1e9), fill = knots), alpha = 0.4) +
  xlab("") + ylab("Biomass (Mt)") +
  facet_wrap(~stratum, scales = "free")

ggsave(file = here(wd, "index_kt_comare.png"), 
       height = 6, width = 10, units = "in")
