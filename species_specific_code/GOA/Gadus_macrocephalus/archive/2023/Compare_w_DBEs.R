##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Create output and diagnostic plots from VAST fits
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(RODBC)
library(getPass)
# library(terra)
# library(RColorBrewer)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import VAST fit
##   Import location of grid cells
##   MBE is the model-based index file
##   DBE is the design-based query from Oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit <- readRDS(paste0("species_specific_code/GOA/Gadus_macrocephalus/",
#                       "hindcast/results/Gadus_macrocephalus_VASTfit.RDS"))
# 
# load(paste0("species_specific_code/GOA/Gadus_macrocephalus/hindcast/",
#             "results/Kmeans_extrapolation-2000.RData"))
# grid_locs <- as.data.frame(Kmeans$centers)
# xrange <- range(grid_locs$E_km)
# xrange_diff <- diff(x = range(grid_locs$E_km))
# yrange <- range(grid_locs$N_km)
# yrange_diff <- diff(x = range(grid_locs$N_km))

MBE <- read.csv(file = paste0("species_specific_code/GOA/Gadus_macrocephalus/",
                              "results/Index.csv"))
MBE <- subset(x = MBE, Estimate > 0)
names(MBE) <- c("Category", "Time", "Stratum", "Units",
                "Estimate", "SE", "CV")

source("R/get_connected.R")

DBE <- RODBC::sqlQuery(
  channel = channel, 
  query = paste("SELECT YEAR, SPECIES_CODE, TOTAL_BIOMASS, BIOMASS_VAR",
                "FROM GOA.BIOMASS_TOTAL WHERE",
                "SPECIES_CODE = 21720 AND",
                "YEAR >= 1990"))
DBE <- DBE[order(DBE$YEAR), ]
DBE$BIOMASS_SD <- sqrt(DBE$BIOMASS_VAR)
DBE$BIOMASS_CV <- DBE$BIOMASS_SD / DBE$TOTAL_BIOMASS

# Using gapindex instead of GOA schema:
dat <- gapindex::get_data(year_set = unique(DBE$YEAR),
                          survey_set = "GOA",
                          spp_codes = 21720,
                          abundance_haul = "Y",
                          pull_lengths = TRUE)
cpue <- gapindex::calc_cpue(racebase_tables = dat)
biomass_stratum <- gapindex::calc_biomass_stratum(racebase_tables = dat,cpue = cpue)
biomass_subareas <- gapindex::calc_biomass_subarea(racebase_tables = dat,
                                                   biomass_strata = biomass_stratum)
DBE <- subset(biomass_subareas, AREA_ID==99903)
DBE <- DBE[,c('YEAR', 'SPECIES_CODE', 'BIOMASS_MT', 'BIOMASS_VAR')]
names(DBE)[names(DBE) == "BIOMASS_MT"] <- "TOTAL_BIOMASS"
DBE <- DBE[order(DBE$YEAR),]
DBE$BIOMASS_SD <- sqrt(DBE$BIOMASS_VAR)
DBE$BIOMASS_CV <- DBE$BIOMASS_SD / DBE$TOTAL_BIOMASS


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Compare MBEs vs DBE of Index and CV
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(2, 1), mar = c(2, 5, 1, 1), oma = c(2, 0,0,0))

plot(TOTAL_BIOMASS/1000 ~ YEAR, data = DBE, type = "p", pch = 16,
     ylim = c(0, 1000), las = 1, 
     xlab = "Year", ylab = "Biomass (000s mt)")
lines(TOTAL_BIOMASS/1000 ~ YEAR, data = DBE)
segments(x0 = DBE$YEAR,
         x1 = DBE$YEAR,
         y0 = DBE$TOTAL_BIOMASS/1000 - DBE$BIOMASS_SD/1000,
         y1 = DBE$TOTAL_BIOMASS/1000 + DBE$BIOMASS_SD/1000)

points(Estimate/1e6 ~ Time,
       data = MBE,
       pch = 16,
       col = "red")
lines(Estimate/1e6 ~ Time,
      data = MBE,
      pch = 16,
      col = "red")
segments(x0 = MBE$Time,
         x1 = MBE$Time,
         y0 = MBE$Estimate/1e6 - MBE$SE/1e6,
         y1 = MBE$Estimate/1e6 + MBE$SE/1e6,
         col = "red")

legend("topleft", 
       legend = c("DBE (+/- 1 SD)", "MBE (+/- 1 SD)"), 
       pch = 16, lty = 1, col = c("black", "red"),
       bty = "n")

plot(BIOMASS_CV ~ YEAR, data = DBE, type = "p", pch = 16,
     ylim = c(0, 0.33), las = 1, 
     xlab = "Year", ylab = "Coefficient of Variation (CV)")
lines(BIOMASS_CV ~ YEAR, data = DBE)

points(CV ~ Time, data = MBE, pch = 16, col = "red")
lines(CV ~ Time, data = MBE, col = "red")
legend("topleft", 
       legend = c("DBE CV", "MBE CV"), 
       pch = 16, lty = 1, col = c("black", "red"),
       bty = "n")
mtext(side = 1, text = "Year", outer = TRUE, line = 0.5)

# ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##   Average Spatial Effects 
# ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ######################################
# ## Spatial Effects
# ######################################
# colors = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral"))
# 
# for (omegatype in 1:2){ #Two types for the 0/1 and pos components
#   
#   ## Extract spatial component
#   scaled_var = list(fit$Report$Omega1_gc[, 1],
#                     fit$Report$Omega2_gc[, 1])[[omegatype]]
#   
#   ## Plot spatial effect
#   goa <- terra::vect(x = cbind(grid_locs, 
#                                Omega1 = fit$Report$Omega1_gc[, 1],
#                                Omega2 = fit$Report$Omega2_gc[, 1]), 
#                      geom =  c("E_km", "N_km"), 
#                      crs = "+proj=utm +zone=5N" )
#   
#   goa_ras <- terra::rast(x = goa, 
#                          resolution = 20)
#   goa_ras <- terra::rasterize(x = goa, 
#                               y = goa_ras, 
#                               field = c("Omega1", "Omega2")[omegatype])
#   
#   zlim_ <- max(abs(scaled_var))
#   plot(goa_ras, 
#        col = colors, 
#        asp = 1,
#        legend = F, 
#        las = 1)
#   
#   # Legend
#   plotrix::color.legend(
#     xl = xrange[1] + xrange_diff * 0.4,
#     xr = xrange[1] + xrange_diff * 0.9,
#     yb = yrange[1] + yrange_diff * 0.05,
#     yt = yrange[1] + yrange_diff * 0.15,
#     legend = pretty(((-ceiling(zlim_)):(ceiling(zlim_))), n = 3),
#     rect.col = colorRampPalette(colors)(1000) ,
#     gradient = "x",
#     align = "rb",
#     cex = 0.75)
# }
# 
# colors = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral"))
# years_idx <- sort(unique(fit$data_frame$t_i))
# 
# par(mfrow = c(5, 3), mar = c(0, 0, 0, 0))
# 
# for (eps_type in 1:2){ #Two types for the 0/1 and pos components
#   
#   ## Extract spatial component
#   scaled_var = list(fit$Report$Epsilon1_gct[, 1, paste0(years_idx)],
#                     fit$Report$Epsilon2_gct[, 1, paste0(years_idx)])[[eps_type]]
#   zlim_ <- max(abs(scaled_var))
#   
#   ## Plot spatial effect
#   goa <- terra::vect(x = cbind(grid_locs, 
#                                as.data.frame(scaled_var)), 
#                      geom =  c("E_km", "N_km"), 
#                      crs = "+proj=utm +zone=5N" )
#   
#   for (iyear in paste0(years_idx)) {
#     goa_ras <- terra::rast(x = goa, 
#                            resolution = 20)
#     goa_ras <- terra::rasterize(x = goa, 
#                                 y = goa_ras, 
#                                 field = iyear)
#     
#     image(goa_ras, 
#           col = colors, 
#           asp = 1,
#           las = 1, axes = F)
#     
#     text(x = xrange[1] + xrange_diff * 0.7,
#          y = yrange[1] + yrange_diff * 0.6, 
#          labels = iyear)
#   }
#   # Legend
#   plotrix::color.legend(
#     xl = xrange[1] + xrange_diff * 0.45,
#     xr = xrange[1] + xrange_diff * 0.95,
#     yb = yrange[1] + yrange_diff * 0.15,
#     yt = yrange[1] + yrange_diff * 0.25,
#     legend = pretty(((-ceiling(zlim_)):(ceiling(zlim_))), n = 3),
#     rect.col = colorRampPalette(colors)(1000) ,
#     gradient = "x",
#     align = "rb",
#     cex = 0.75)
# }
# 
# par(mfrow = c(5, 3), mar = c(0, 0, 0, 0))
# ## Density Plots
# vals  = fit$Report$D_gct[, 1, paste0(years_idx)]
# val_cuts = quantile(x = vals,
#                     probs = seq(0, 1, length = 11))
# val_cuts_legend <- round(x = val_cuts, digits = 0)
# 
# ## Plot spatial effect
# goa <- terra::vect(x = cbind(grid_locs, 
#                              as.data.frame(scaled_var)), 
#                    geom =  c("E_km", "N_km"), 
#                    crs = "+proj=utm +zone=5N" )
# 
# for (iyear in paste0(years_idx)) {
#   goa_ras <- terra::rast(x = goa, 
#                          resolution = 20)
#   goa_ras <- terra::rasterize(x = goa, 
#                               y = goa_ras, 
#                               field = iyear)
#   
#   colors = c(brewer.pal(n = 9, name = "Oranges"), "black")
#   
#   image(goa_ras, 
#         col = colors, 
#         asp = 1,
#         las = 1, axes = F)
#   
#   text(x = xrange[1] + xrange_diff * 0.7,
#        y = yrange[1] + yrange_diff * 0.6, 
#        labels = iyear)
# }
# # Legend
# legend(x = xrange[1] + xrange_diff * 0.5,
#        y = yrange[1] + yrange_diff * 0.5,
#        fill = colors,
#        bty = "n",
#        ncol = 3,
#        cex = 0.7,
#        legend = c(paste0("< ", val_cuts_legend[2]),
#                   paste0(val_cuts_legend[2:(length(val_cuts_legend)-1)], "-",
#                          val_cuts_legend[3:length(val_cuts_legend)])) )
