# Standardized plots of center of gravity and effective area occupied produced
# from the VAST models of abundance and age composition.
# By: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.04.04

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

# SET-UP ----------------------------------------------------------------------
this_year <- 2024

# Read in VAST results - update for each species 
workDir <- here("species_specific_code", "BS", "pollock", "results") 
# TODO: structure all species folders to have explicit hindcast and production
# folders each containing a "results" folder for the index model, e.g.,
# workDir <- here("species_specific_code", "BS", "pacific_cod", "production", "results") 
# workDir <- here("VAST_results", "BS", "nrs")

# Reading these in will take a moment! ...note that these subfolders should be
# standardized to results (for the index), and results_age (for comps), and 
# we want the former
VAST_results <- readRDS(here(workDir, "VASTresults.RDS"))  # for COG 
VAST_fit <- readRDS(here(workDir, "VASTfit.RDS"))  # for EAO

# Make a results object 
# TODO: test this a couple more times. Doesn't seem to work consistently
# VAST_results <- plot_results(VAST_fit, 
#                              zrange = c(-3,3),
#                              n_cells = 600, 
#                              strata_names = c("Both", "EBS", "NBS"), 
#                              check_residuals=TRUE,
#                              n_samples=0)

# Create directory for saving output (should we change to underscores and no caps,
# for the folder maybe? Maybe just change it to "esp")
saveDir <- here(workDir, paste0(this_year, " ESR products"))
dir.create(saveDir, showWarnings = FALSE)


# Center of gravity -----------------------------------------------------------
cog <- function(results = VAST_results, dir = saveDir, save_data = TRUE, save_plots = TRUE) {
  cog <- data.frame(results$Range$COG_Table)
  cog$Year <- as.numeric(cog$Year)
  cog$COG_hat <- as.numeric(cog$COG_hat)
  cog$SE <- as.numeric(cog$SE)
  
  if(grepl("nrs", dir) == TRUE) {
    # Remove years before 1996 if stock is Northern rock sole (due to change in IDs)
    cog <- cog %>% filter(Year >= 1996)
  }
  
  cog$m[cog$m == 1] <- "Easting (km)"
  cog$m[cog$m == 2] <- "Northing (km)"
  
  # Original plot - UTM
  ts <- ggplot(cog %>% filter(Year != 2020), aes(x = Year, y = COG_hat)) +
    geom_line() +
    geom_ribbon(aes(ymin = (COG_hat - SE), ymax = (COG_hat + SE)), alpha = 0.2) +
    xlab("") + ylab("") +
    facet_wrap(~ m, ncol = 1, scales = "free_y") +
    theme(plot.background = element_rect(fill = "transparent"))
  
  # Convert to lat/long using akgfmaps package and plot
  cog_east <- cog %>% filter(m == "Easting (km)" & Year != 2020)
  cog_east <- cog_east[, 2:4]  # get rid of ID column and not sure what to do with SE!
  colnames(cog_east)[2] <- "Easting (km)"
  
  cog_north <- cog %>% filter(m == "Northing (km)" & Year != 2020)
  cog_north <- cog_north[, 2:4]
  colnames(cog_north)[2] <- "Northing (km)"
  
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
  map <- ggplot(data = world) +
    geom_sf() +
    geom_point(data = cog_latlon,
               aes(x = X, y = Y, color = Year), size = 1) +
    # geom_errorbar(data = cog_error,
    #               aes(x = X, y = Y, ymin = ymin,ymax = ymax, color = Year), alpha = 0.8) +
    # geom_errorbarh(data = cog_error,
    #                aes(x = X, y = Y, xmin = xmin,xmax = xmax, color = Year), alpha = 0.8) +
    coord_sf(xlim = c(-179, -157), ylim = c(54, 65), expand = FALSE) +
    scale_color_viridis(option = "plasma", discrete = FALSE, end = 0.9) +
    scale_x_continuous(breaks = c(-178, -158)) +
    scale_y_continuous(breaks = c(55, 64)) +
    labs(x = NULL, y = NULL) +
    theme(plot.background = element_rect(fill = "transparent"))
  
  cog_error2 <- cog_error[-1, c(1:2)]
  cog_error2[nrow(cog_error2) + 1, ] <- NA
  colnames(cog_error2) <- c("X2", "Y2")
  cog_error2 <- cbind.data.frame(cog_error, cog_error2) 
  
  # Plot as scatter (sparkleplot)
  sparkle <- ggplot(data = cog_error, aes(x = X, y = Y, color = Year)) +
    geom_point() +
    # With arrow
    # geom_segment(data = cog_error2 %>% filter(Year >= this_year - 10), 
    #              aes(x = X, y = Y, xend = X2, yend = Y2), 
    #              alpha = 0.8, arrow = arrow(length = unit(0.03, "npc"))) +
    # Without arrow
    geom_segment(data = cog_error2 %>% filter(Year >= this_year - 10), 
                 aes(x = X, y = Y, xend = X2, yend = Y2), 
                 alpha = 0.8) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax, color = Year), alpha = 0.4) +
    geom_errorbarh(aes(xmin = xmin, xmax = xmax, color = Year), alpha = 0.4) +
    scale_color_viridis(option = "plasma", discrete = FALSE, end = 0.9) +
    xlab("Longitude (°W)") + ylab("Latitude (°N)") +
    theme(plot.background = element_rect(fill = "transparent"))
  
  # Inset map into sparkleplot
  inset <- ggdraw() +
    draw_plot(plot = sparkle) +
    draw_plot(plot = map +
                theme(legend.position = "none") +
                guides(x = "none", y = "none"), 
              x = 0.81, y = 0.75, width = 0.21, height = 0.21) 
  
  # Combine sparkleplot and time-series plot
  all <- plot_grid(inset, ts)
  all
  
  # Another inset for using independently
  inset2 <- ggdraw() +
    draw_plot(plot = sparkle) +
    draw_plot(plot = map +
                theme(legend.position = "none") +
                guides(x = "none", y = "none"), 
              x = 0.60, y = 0.76, width = 0.21, height = 0.21) 
  
  if(save_data == TRUE) {
    # Save COG as UTM (easting/northing)
    write.csv(cog, file = here(dir, "COG_utm.csv"), row.names = FALSE)
    
    # Save COG as lat/long (without error)
    cog_latlon <- cog_latlon[, c(3, 2, 1)]
    colnames(cog_latlon) <- c("Year", "Latitude", "Longitude")
    write.csv(cog_latlon, file = here(dir, "COG_latlong.csv"), row.names = FALSE)
  }
  
  if(save_plots == TRUE) {
    ggsave(ts, filename = here(dir, "COG_utm.png"),
           width = 150, height = 180, unit = "mm", dpi = 300, bg = "white")
    ggsave(map, filename = here(dir, "COG_map.png"),
           width = 110, height = 90, unit = "mm", dpi = 300,  bg = "white")
    ggsave(sparkle, filename = here(dir, "COG_scatter.png"),
           width = 130, height = 100, unit = "mm", dpi = 300,  bg = "white")
    ggsave(all, filename = here(dir, "COG_all.png"),
           width = 250, height = 100, unit = "mm", dpi = 300,  bg = "white")
    ggsave(inset2, filename = here(dir, "COG_inset.png"),
           width = 130, height = 100, unit = "mm", dpi = 300,  bg = "white")
  }
  
  return(list(table = cog_latlon, table_error = cog_error, 
              ts = ts, map = map, sparkle = sparkle, all = all, inset = inset2))
}

cog_plots <- cog()
# Map insert may look funny here because of the dimensions of the Rstudio plotting window. Check saved plot!
cog_plots$all  
cog_plots$map

# Effective area occupied -----------------------------------------------------
options(scipen = 999)

eao <- function(fit = VAST_fit, dir = saveDir, save_data = TRUE, save_plot = TRUE) {
  # Get EAO estimate from VAST fit object
  report <- TMB::summary.sdreport(fit$parameter_estimates$SD)
  area <- report[which(rownames(report) == "log_effective_area_ctl"), c('Estimate', 'Std. Error')]
  Year <- sort(unique(fit$year_labels))
  area <- as.data.frame(cbind(area, Year))
  area <- area[which(area$Year %in% unique(fit$data_frame$t_i)), ]
  
  # Reshape and plot EAO
  area$Region <- c(rep("Both", length(1982:this_year) - 1), 
                   rep("EBS", length(1982:this_year) - 1),
                   rep("NBS", length(1982:this_year) - 1))
  colnames(area)[2] <- "error"
  area$Estimate <- as.numeric(area$Estimate)
  area$error <- as.numeric(area$error)
  area$Year <- as.numeric(area$Year)
  
  area$Estimate <- exp(area$Estimate)
  area$error <- area$error
  
  area <- area %>%
    mutate(ymin = Estimate - (Estimate * error),
           ymax = Estimate + (Estimate * error))
  
  if(grepl("nrs", dir) == TRUE) {
    # Remove years before 1996 if stock is Northern rock sole (due to change in IDs)
    area <- area %>% filter(Year >= 1996)
  }
  
  plot <- ggplot(area %>% filter(Year != 2020), aes(x = Year, y = Estimate)) +
    geom_line(alpha = 0.8, aes(color = Region, linetype = Region)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = Region), alpha = 0.3) +
    scale_y_continuous(labels = scales::comma, limits = c(0, NA)) +
    scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.7) +
    scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.7) +
    ylab(expression(paste("Effective Area Occupied (", km^2, ")")))
  
  if(save_data == TRUE) {
    area <- area[, c(3, 4, 1, 2)]
    write.csv(area, file = here(dir, "Ln_effective_area.csv"), row.names = FALSE)
  }
  
  if(save_plot == TRUE) {
    ggsave(plot, filename = here(dir, "effective_area.png"),
           width = 150, height = 100, unit = "mm", dpi = 300)
  }
  
  return(plot)
}

eao_plot <- eao()
eao_plot

############ Stop here for single-species plotting! ##############

# Combine regional COG into one map -------------------------------------------
# TODO: run after standardizing folder structures as noted above
# Read in each species in the region and combine into a list
bs_pol_res <- readRDS(here("species_specific_code", "BS", "pollock", "results", "VAST Index", "VASTresults.RDS"))
bs_cod_res <- readRDS(here("VAST_results", "BS", "pcod", "VASTresults.RDS"))
bs_yel_res <- readRDS(here("VAST_results", "BS", "yfs", "VASTresults.RDS"))
bs_nrs_res <- readRDS(here("VAST_results", "BS", "nrs", "VASTresults.RDS"))

bs_res <- list(bs_pol_res, bs_cod_res, bs_yel_res, bs_nrs_res) 
species_bs <- c("Walleye pollock", "Pacific cod", "Yellowfin sole", "Northern rock sole")
sp2 <- c("Walleye pollock", "Pacific cod", "Yellowfin sole", "nrs")  # hack to get NRS filtering to work

# Run the COG function for each species and combine together
cog_bs <- data.frame()
cog_bs_error <- data.frame()
for(i in 1:length(bs_res)) {
  out <- cog(results = bs_res[[i]], save_data = FALSE, save_plots = FALSE)
  df <- out$table
  df$Species <- species_bs[i]
  cog_bs <- rbind.data.frame(cog_bs, df)
  
  df_error <- out$table_error
  df_error$Species <- species_bs[i]
  cog_bs_error <- rbind.data.frame(cog_bs_error, df_error)
}

# Plot all species on one plot
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)  # turn off spherical geometry
cog_bs_map <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = cog_bs %>% filter(Year == this_year),
             aes(x = X, y = Y, color = Species, shape = Species)) +
  geom_errorbar(data = cog_bs_error %>% filter(Year == this_year),
                aes(x = X, y = Y, ymin = ymin,ymax = ymax, color = Species), 
                alpha = 0.8, width = 0) +
  geom_errorbarh(data = cog_bs_error %>% filter(Year == this_year),
                 aes(x = X, y = Y, xmin = xmin,xmax = xmax, color = Species), 
                 alpha = 0.8, height = 0) +
  coord_sf(xlim = c(-179, -157), ylim = c(54, 65), expand = FALSE) +
  scale_color_viridis(option = "plasma", discrete = TRUE, end = 0.9) +
  xlab(" ") + ylab(" ") + 
  ggtitle("Center of Gravity for 2023") +
  theme(legend.title=element_blank())
cog_bs_map

ggsave(cog_bs_map, filename = here("VAST_results", "BS", "COG_bs_map.png"), 
       width = 110, height = 70, unit = "mm", dpi = 300)

# Combine all COG plots together ----------------------------------------------
cog_bs_plots <- list()
for(i in 1:length(bs_res)) {
  out <- cog(results = bs_res[[i]], dir = sp2[i], save_data = FALSE, save_plots = FALSE)
  cog_bs_plots[[i]] <- out$all
}

cog_combined <- plot_grid(plotlist = cog_bs_plots, labels = species_bs,
                          ncol = 1, 
                          label_size = 10, label_x = 0.5, hjust = 0.5, label_fontface = "plain")
cog_combined

ggsave(cog_combined, filename = here("VAST_results", "BS", "COG_bs_all.png"), 
       width = 240, height = 300, unit = "mm", dpi = 300)

# Combine all EAO plots together ----------------------------------------------
# Read in each species in the region and combine into a list
bs_pol_fit <- readRDS(here("species_specific_code", "BS", "pollock", "results", "VAST Index", "VASTfit.RDS"))
bs_cod_fit <- readRDS(here("VAST_results", "BS", "pcod", "VASTfit.RDS"))
bs_yel_fit <- readRDS(here("VAST_results", "BS", "yfs", "VASTfit.RDS"))
bs_nrs_fit <- readRDS(here("VAST_results", "BS", "nrs", "VASTfit.RDS"))

bs_fit <- list(bs_pol_fit, bs_cod_fit, bs_yel_fit, bs_nrs_fit) 

eao_bs_plots <- list()
for(i in 1:length(bs_fit)) {
  out <- eao(fit = bs_fit[[i]], dir = sp2[i], save_data = FALSE, save_plot = FALSE)
  plot <- out + ggtitle(species_bs[i])
  eao_bs_plots[[i]] <- plot
}

eao_combined <- plot_grid(plotlist = eao_bs_plots)
eao_combined

ggsave(eao_combined, filename = here("VAST_results", "BS", "EAO_bs.png"), 
       width = 220, height = 130, unit = "mm", dpi = 300)
