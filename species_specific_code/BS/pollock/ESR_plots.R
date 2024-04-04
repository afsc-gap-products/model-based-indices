# Script for plots of VAST model output for pollock that are used in the Bering 
# Sea Ecosystem Status Reports
# By: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.03.14

library(here)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)
library(sf)
library(rnaturalearth)
library(akgfmaps)
library(cowplot)

# Set ggplot theme
devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
theme_set(theme_sleek())

this_year <- 2023

# Read in VAST results --------------------------------------------------------
workDir <- here("species_specific_code", "BS", "pollock")
results <- readRDS(here(workDir, "results", "VAST Index", "VASTresults.RDS"))  # for COG
full_fit <- readRDS(here(workDir, "results", "VAST Index", "VASTfit_full.RDS"))  # for EAO

# Center of gravity -----------------------------------------------------------
cog <- data.frame(results$Range$COG_Table)
cog$Year <- as.numeric(cog$Year)
cog$COG_hat <- as.numeric(cog$COG_hat)
cog$SE <- as.numeric(cog$SE)
# write.csv(cog, file = here(workDir, "results", "ESR products", "COG.csv"), row.names = FALSE)

cog$m[cog$m == 1] <- "Easting"
cog$m[cog$m == 2] <- "Northing"

# Original plot - UTM
cog_utm <- ggplot(cog %>% filter(Year != 2020), aes(x = Year, y = COG_hat)) +
  geom_line(alpha = 0.4) +
  geom_pointrange(aes(ymin = (COG_hat - SE), ymax = (COG_hat + SE)), size = 0.3) +
  ylab(" ") +
  facet_wrap(~ m, ncol = 1, scales = "free_y")
cog_utm

# Convert to lat/long using akgfmaps package and plot
cog_east <- cog %>% filter(m == "Easting" & Year != 2020)
cog_east <- cog_east[, 2:4]  # get rid of ID column and not sure what to do with SE!
colnames(cog_east)[2] <- "Easting"

cog_north <- cog %>% filter(m == "Northing" & Year != 2020)
cog_north <- cog_north[, 2:4]
colnames(cog_north)[2] <- "Northing"

cog_latlon <- cbind.data.frame(X = cog_east[, 2], Y = cog_north[, 2])
# CRS information for VAST outputs here: 
# https://github.com/James-Thorson-NOAA/FishStatsUtils/blob/main/R/project_coordinates.R
cog_latlon <- transform_data_frame_crs(cog_latlon, 
                                       coords = c("X", "Y"), 
                                       in.crs = "+proj=utm +datum=WGS84 +units=km +zone=2",
                                       out.crs = "+proj=longlat +datum=WGS84")
cog_latlon$Year <- cog_east$Year

# Include error in COG estimate before transformation to get min & max values
cog_min <- cbind.data.frame(X = cog_east[, 2] - cog_east[, 3], 
                            Y = cog_north[, 2] - cog_north[, 3])
cog_min <- transform_data_frame_crs(cog_min, 
                                    coords = c("X", "Y"), 
                                    in.crs = "+proj=utm +datum=WGS84 +units=km +zone=2",
                                    out.crs = "+proj=longlat +datum=WGS84")

cog_max <- cbind.data.frame(X = cog_east[, 2] + cog_east[, 3], 
                            Y = cog_north[, 2] + cog_north[, 3])
cog_max <- transform_data_frame_crs(cog_max, 
                                    coords = c("X", "Y"), 
                                    in.crs = "+proj=utm +datum=WGS84 +units=km +zone=2",
                                    out.crs = "+proj=longlat +datum=WGS84")
cog_error <- cbind.data.frame(cog_latlon, 
                              xmin = cog_min$X, xmax = cog_max$X,
                              ymin = cog_min$Y, ymax = cog_max$Y)

# Plot on a map
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)  # turn off spherical geometry
cog_map <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = cog_latlon,
             aes(x = X, y = Y, color = Year), size = 1) +
  # geom_errorbar(data = cog_error,
  #               aes(x = X, y = Y, ymin = ymin,ymax = ymax, color = Year), alpha = 0.8) +
  # geom_errorbarh(data = cog_error,
  #                aes(x = X, y = Y, xmin = xmin,xmax = xmax, color = Year), alpha = 0.8) +
  coord_sf(xlim = c(-179, -157), ylim = c(54, 65), expand = FALSE) +
  scale_color_viridis(option = "plasma", discrete = FALSE, end = 0.9) +
  xlab(" ") + ylab(" ") 
cog_map

# Plot as scatter (sparkleplot)
cog_scatter <- ggplot(cog_error, aes(x = X, y = Y, color = Year)) +
  geom_point() +
  geom_errorbar(aes(ymin = ymin,ymax = ymax, color = Year), alpha = 0.7) +
  geom_errorbarh(aes(xmin = xmin,xmax = xmax, color = Year), alpha = 0.7) +
  scale_color_viridis(option = "plasma", discrete = FALSE, end = 0.9) +
  xlab("Longitude (°W)") + ylab("Latitude (°N)")
cog_scatter

# Add map as insert to sparkleplot
cog_inset <- ggdraw() +
  draw_plot(cog_scatter) +
  draw_plot(cog_map + 
              theme(legend.position = "none") +
              guides(x = "none", y = "none") +
              theme(plot.background = element_rect(fill = "transparent")), 
            x = 0.77, y = 0.7, width = 0.25, height = 0.25)
cog_inset

# Combine sparkleplot and time-series plot
cog_all <- plot_grid(cog_inset, cog_utm)
cog_all

# Effective area occupied -----------------------------------------------------
options(scipen = 999)

# Get EAO estimate from VAST fit object
report <- TMB::summary.sdreport(full_fit$parameter_estimates$SD)
area <- report[which(rownames(report) == "log_effective_area_ctl"), c('Estimate', 'Std. Error')]
Year <- sort(unique(full_fit$year_labels))
area <- as.data.frame(cbind(area, Year))
area <- area[which(area$Year %in% unique(full_fit$data_frame$t_i)), ]
# write.csv(area, file = here(workDir, "results", "ESR products", "ln_effective_area.csv"),
#           row.names = FALSE)

area$Region <- c(rep("Both", length(1982:this_year) - 1), 
                 rep("EBS", length(1982:this_year) - 1),
                 rep("NBS", length(1982:this_year) - 1))
colnames(area)[2] <- "error"
area$Estimate <- as.numeric(area$Estimate)
area$error <- as.numeric(area$error)
area$Year <- as.numeric(area$Year)

area$Estimate <- exp(area$Estimate)
area$error <- area$error

area_plot <- ggplot(area %>% filter(Year != 2020), 
                    aes(x = Year, y = Estimate, 
                        color = Region, shape = Region)) +
  geom_line(alpha = 0.4) +
  geom_pointrange(aes(ymin = (Estimate - (Estimate * error)), 
                      ymax = (Estimate + (Estimate * error))),
                  alpha = 0.8) +
  scale_y_continuous(labels = scales::comma, limits = c(0, NA)) +
  scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.7) +
  ylab(expression(paste("Effective Area Occupied (", km^2, ")")))
area_plot

# Save plots ------------------------------------------------------------------
ggsave(cog_utm, filename = here(workDir, "results", "ESR products", "2023_pollock_COG_utm.png"),
       width = 150, height = 180, unit = "mm", dpi = 300)
ggsave(cog_map, filename = here(workDir, "results", "ESR products", "2023_pollock_COG_map.png"),
       width = 110, height = 90, unit = "mm", dpi = 300)
ggsave(cog_scatter, filename = here(workDir, "results", "ESR products", "2023_pollock_COG_scatter.png"),
       width = 130, height = 100, unit = "mm", dpi = 300)
ggsave(cog_inset, filename = here(workDir, "results", "ESR products", "2023_pollock_COG_inset.png"),
       width = 130, height = 100, unit = "mm", dpi = 300)
ggsave(cog_all, filename = here(workDir, "results", "ESR products", "2023_pollock_COG_all.png"),
       width = 240, height = 100, unit = "mm", dpi = 300)
ggsave(area_plot, filename = here(workDir, "results", "ESR products", "2023_pollock_area.png"),
       width = 150, height = 100, unit = "mm", dpi = 300)
