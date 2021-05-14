###############################################################################
## Project:       VAST Model for Center of Gravity (COG) and Effective Area
##                for Gulf of Alaska walleye pollock and Pacific cod. 
##
## Author:        Zack Oyafuso, zack.oyafuso@noaa.gov with code adapted from
##                Cecilia O'Leary, cecilia.oleary@noaa.gov
##                Lewis Barnett, lewis.barnett@noaa.gov
##                James Thorson, james.thorson@noaa.gov
##
## Description:   Code to produce model-based indices using VAST software 
##                implemented in R, with adjusted COG output specific to ESP 
##                needs (data west of 140W, area occupied and COG with NE 
##                rotation, 1984-2019), originally created Oct 2020
##
## Pkgs Versions: 
                  rm(list = ls())
                  r_version <- "R version 4.0.2 (2020-06-22)"
                  vast_version <- "3.6.1"
                  fishstatsutils_version <- "2.8.0"
                  vast_cpp_version <- "VAST_v12_0_0"
###############################################################################
###############################################################################


##################################################
####  Set directories
##################################################
wd <- "G:/Oyafuso/ModSquad_2021_Oyafuso/"

##################################################
####  Load VAST package
##################################################
library(VAST)

##################################################
####  Make sure versions are correct
##################################################
ifelse(test = sessionInfo()$R.version$version.string == r_version &
         packageVersion("VAST") == vast_version &
         packageVersion("FishStatsUtils") == fishstatsutils_version &
         FishStatsUtils::get_latest_version() == vast_cpp_version,
       yes = "versions are good to go", 
       no = "check your versions")

##################################################
####  Load the data for VAST
##################################################
data_geostat <- read.csv(file = paste0(wd, 
                                       "data/AK_BTS_GOA_PCod_pollock_140.csv"))
spp_names <- read.csv(file = paste0(wd, 
                                    "data/species.csv") )

##################################################
####  Filter data to species of interest
##################################################
which_spp <- c("Pacific cod", "walleye pollock")[1]
spp_code <- spp_names$SPECIES_CODE[spp_names$COMMON_NAME == which_spp]
data_geostat <- subset(x = data_geostat, 
                       subset = SPECIES_CODE == spp_code)

##################################################
####  Save data_geostat
##################################################
saveRDS(object = data_geostat,
        file = paste0(wd, 
                      "data/data_geostat_W140_", 
                      gsub(x = which_spp,
                           pattern = " ", 
                           replacement = "_"), 
                      ".RDS"))

##################################################
#### Load extrapolation grid and create extrapolation grid in the VAST format
#### Extrapolation grid area is in m^2 and is converted to km^2
#### Subset input grid to those west of 140 longitude
##################################################
goa_grid <- read.csv(paste0(wd, "data/GOAThorsonGrid_Less700m.csv") )
input_grid <- cbind(Lat = goa_grid$Lat,
                    Lon = goa_grid$Lon,
                    Area_km2 = goa_grid$Shape_Area / 1000000)  
input_grid <- input_grid[input_grid[, "Lon"] <= -140, ]

##################################################
####  Customize strata limits to limit extrapolation grid to west of 140
##################################################
strata_limits <- data.frame(STRATA = "west_of_140W",
                            west_border = -Inf,
                            east_border = -140)

##################################################
####  Make VAST settings
##################################################
settings <- FishStatsUtils::make_settings(Version = vast_cpp_version,
                                          n_x = 750, 
                                          Region = "User", 
                                          purpose = "index2", 
                                          ObsModel = c(2, 1), 
                                          strata.limits = strata_limits)

##################################################
#### Build object but don't run
##################################################
fit <- fit_model(settings = settings, 
                 working_dir = paste0(wd, "output/", which_spp, "/"),
                 Lat_i = data_geostat[, "LATITUDE"], 
                 Lon_i = data_geostat[, "LONGITUDE"], 
                 t_i = data_geostat[, "YEAR"], 
                 b_i = data_geostat[, "BIOMASS"], 
                 a_i = data_geostat[, "EFFORT"], 
                 input_grid = input_grid, 
                 knot_method = "grid",
                 optimize_args = list("lower" = -Inf, "upper" = Inf),
                 run_model = FALSE)

##################################################
####  Rotation is counter-clockwise from east, where you specify the amount 
####  of rotation: 90 degrees or pi/4 radians; or as requested (Martin Dorn
####  requested for GOA pollock: 40 degrees or 2*pi/9 radians). I used -pi/8
####  radians and centered the coordinates as a test for now.
####  see https://en.wikipedia.org/wiki/Rotation_matrix or 
####  https://mathworld.wolfram.com/RotationMatrix.html for quick math refresh
####  
####  Plot the original and rotated coordinates
##################################################
rotation_radians <-  -pi / 8
Z_gm <- fit$data_list$Z_gm

Zprime_gm <- cbind(
  Z_gm, 
  x_rotated_km = scale(cos(rotation_radians) * Z_gm[, "E_km"] -
                         sin(rotation_radians) * Z_gm[, "N_km"], 
                       scale = FALSE),
  y_rotated_km =  scale(sin(rotation_radians) * Z_gm[, "E_km"] +
                          cos(rotation_radians) * Z_gm[, "N_km"], 
                        scale = FALSE))

par(mfrow = c(3, 1))
plot(x = Zprime_gm[, 1], 
     y = Zprime_gm[, 2], 
     main = "original Z_gm", 
     asp = 1)
plot(x = Zprime_gm[, 3], 
     y = Zprime_gm[, 4], 
     main = "rotated Z_gm", 
     asp = 1)

##################################################
####  Re-build and run
##################################################
fit = fit_model( settings = settings, 
                 working_dir = paste0(wd, "output/", which_spp, "/"),
                 Lat_i = data_geostat[, "LATITUDE"], 
                 Lon_i = data_geostat[, "LONGITUDE"], 
                 t_i = data_geostat[, "YEAR"], 
                 b_i = data_geostat[, "BIOMASS"], 
                 a_i = data_geostat[, "EFFORT"], 
                 input_grid = input_grid, 
                 knot_method = "grid", 
                 optimize_args = list("lower" = -Inf, "upper" = Inf),
                 Z_gm = Zprime_gm )

##################################################
####  Save and plot Center Of Gravity 
##################################################
results <- plot( fit, 
                 working_dir = paste0(wd, "output/", which_spp, "/"),
                 check_residuals = FALSE)

write.csv( results$Range$COG_Table, 
           file = paste0(wd, "output/", which_spp, "/COG.csv"), 
           row.names = FALSE )

par(mfrow = c(2, 1), mar = c(3, 3, 1, 1))
plot(Zprime_gm[, 1:2], asp = 1)
with(results$Range,
     points(COG_Table[COG_Table[, "m"] == 1, "COG_hat"],
            COG_Table[COG_Table[, "m"] == 2, "COG_hat"],
            pch = 16, col = "red"))

plot(Zprime_gm[, 3:4], asp = 1)
with(results$Range,
     points(COG_Table[COG_Table[, "m"] == 3, "COG_hat"],
            COG_Table[COG_Table[, "m"] == 4, "COG_hat"],
            pch = 16, col = "red"))

##################################################
#### Save effective log-area occupied 
##################################################
report <- TMB::summary.sdreport(fit$parameter_estimates$SD)
ln_km2 <- report[which(rownames(report) == "log_effective_area_ctl"),
                 c("Estimate", "Std. Error")]
year <- sort(unique(fit$year_labels))
ln_km2 <- as.data.frame(cbind(ln_km2, year))
ln_km2 <- ln_km2[which(ln_km2$year %in% unique(fit$data_frame$t_i)), ]

write.csv(x = ln_km2, 
          file = paste0(wd, "output/", which_spp, "/ln_effective_area.csv"), 
          row.names = FALSE )

##################################################
#### Save the VAST model fit object 
##################################################
save(fit, file = paste0(wd, "output/", which_spp, "/VASTfit.RData"))
