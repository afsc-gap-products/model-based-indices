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

# Set ggplot theme
devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
theme_set(theme_sleek())

this_year <- 2023

# Center of gravity -----------------------------------------------------------
# Original plot - UTM
cog <- read.csv(here("VAST_results", "COG.csv")) 
cog$m[cog$m == 1] <- "Eastings"
cog$m[cog$m == 2] <- "Northings"

cog_plot <- ggplot(cog %>% filter(Year != 2020), aes(x = Year, y = COG_hat)) +
  geom_line(alpha = 0.4) +
  geom_pointrange(aes(ymin = (COG_hat - SE), ymax = (COG_hat + SE))) +
  ylab("Center of Gravity") +
  facet_wrap(~ m, ncol = 1, scales = "free_y")
cog_plot

# Convert to lat/long and plot
cog_east <- cog %>% filter(m == "Eastings")
cog_east <- cog_east[, 2:3]  # get rid of ID column and not sure what to do with SE!
colnames(cog_east)[2] <- "Eastings"

cog_north <- cog %>% filter(m == "Northings")
cog_north <- cog_north[, 2:3]
colnames(cog_north)[2] <- "Northings"

# Convert to lat & lon with sf
# TODO: figure out correct CRS
cog_latlon <- cbind.data.frame(cog_east, Northings = cog_north[, 2]) %>%
  st_as_sf(coords = c("Eastings", "Northings"), crs = 4326) %>%  # convert to an sf object
  st_transform(4326) %>%  # transform or convert coordinates of sf
  st_coordinates() %>%  # retrieve coordinates in matrix form
  as_tibble() 
cog_latlon$Year <- cog_east$Year

world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)  # turn off spherical geometry

ggplot(data = world) +
  geom_sf() +
  geom_point(data = cog_latlon,
             aes(x = X, y = Y, color = Year)) +
  coord_sf(xlim = c(-179, -157), ylim = c(54, 66), expand = FALSE) +
  scale_color_viridis(option = "plasma", discrete = FALSE, end = 0.9) +
  xlab(" ") + ylab(" ")

# Area occupied ---------------------------------------------------------------
options(scipen = 999)
area <- read.csv(here("VAST_results", "ln_effective_area.csv")) 
area$Region <- c(rep("Both", length(1982:this_year) - 1), 
                 rep("EBS", length(1982:this_year) - 1),
                 rep("NBS", length(1982:this_year) - 1))
colnames(area)[2] <- "error"
area$Estimate <- exp(area$Estimate)
area$error <- area$error

area_plot <- ggplot(area %>% filter(Year != 2020), aes(x = Year, y = Estimate, color = Region, shape = Region)) +
  geom_line(alpha = 0.4) +
  geom_pointrange(aes(ymin = (Estimate - (Estimate * error)), 
                      ymax = (Estimate + (Estimate * error))),
                  alpha = 0.8) +
  scale_y_continuous(labels = scales::comma, limits = c(0, NA)) +
  scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
  ylab("Effective area occupied (km^2)")
area_plot

# Save plots ------------------------------------------------------------------
# ggsave(cog_plot, filename = here("VAST_results", "2023_pollock_COG.png"),
#        width = 150, height = 180, unit = "mm", dpi = 300)
# ggsave(area_plot, filename = here("VAST_results", "2023_pollock_area.png"),
#        width = 150, height = 100, unit = "mm", dpi = 300)
