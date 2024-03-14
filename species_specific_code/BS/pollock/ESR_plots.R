# Script for plots of VAST model output for pollock that are used in the Bering 
# Sea Ecosystem Status Reports
# By: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.03.14

library(here)
library(ggplot2)
library(dplyr)

# Set ggplot theme
devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
theme_set(theme_sleek())

this_year <- 2023

# Center of gravity -----------------------------------------------------------
cog <- read.csv(here("VAST_results", "COG.csv")) 
cog$m[cog$m == 1] <- "Eastings"
cog$m[cog$m == 2] <- "Northings"

cog_plot <- ggplot(cog, aes(x = Year, y = COG_hat)) +
  geom_line(alpha = 0.4) +
  geom_pointrange(aes(ymin = (COG_hat - SE), ymax = (COG_hat + SE))) +
  ylab("Center of Gravity") +
  facet_wrap(~ m, ncol = 1, scales = "free_y")
cog_plot

# Area occupied ---------------------------------------------------------------
options(scipen = 999)
area <- read.csv(here("VAST_results", "ln_effective_area.csv")) 
area$Region <- c(rep("Both", length(1982:this_year) - 1), 
                 rep("EBS", length(1982:this_year) - 1),
                 rep("NBS", length(1982:this_year) - 1))
colnames(area)[2] <- "error"
area$Estimate <- exp(area$Estimate)
area$error <- area$error

area_plot <- ggplot(area, aes(x = Year, y = Estimate, color = Region, shape = Region)) +
  geom_line(alpha = 0.4) +
  geom_pointrange(aes(ymin = (Estimate - (Estimate * error)), 
                      ymax = (Estimate + (Estimate * error))),
                  alpha = 0.8) +
  scale_y_continuous(labels = scales::comma, limits = c(0, NA)) +
  scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
  ylab("Effective area occupied (km^2)")
area_plot

# Save plots ------------------------------------------------------------------
ggsave(cog_plot, filename = here("VAST_results", "2023_pollock_COG.png"),
       width = 150, height = 180, unit = "mm", dpi = 300)
ggsave(area_plot, filename = here("VAST_results", "2023_pollock_area.png"),
       width = 150, height = 100, unit = "mm", dpi = 300)
