# Run all Bering Sea age compositions 
# Author: Sophia Wassermann
# Date: 26-3-2025

library(tinyVAST)
library(fmesher)
library(sf)
library(here)
library(dplyr)


# Set up ----------------------------------------------------------------------
phase <- c("hindcast", "production")[1] # specify analysis phase

sp <- 2 # specify species from species vector
species <- c("yellowfin_sole", "pollock", "pacific_cod")[sp]

# Set year
this_year <- as.numeric(format(Sys.Date(), "%Y"))
if(phase == "hindcast") {this_year <- this_year - 1}  

# Set working directory specific to species & phase
workDir <- here("species_specific_code", "BS", species, phase)

# Read in and format data -----------------------------------------------------
if(species == "pollock"){
  dat <- read.csv(here(workDir, "data", paste0("VAST_ddc_alk_", this_year, ".csv")))  
  dat <- transmute(dat,
                   cpue = CPUE_num * 100, # converts cpue from kg/ha to kg/km^2
                   year = as.integer(Year),
                   lat = Lat,
                   lon = Lon,
                   age = Age
  )
}

if(species %in% c("yellowfin_sole", "pacific_cod")){
  dat <- readRDS(here(workDir, "data", "data_geostat_agecomps.RDS"))  
  dat <- rename(dat, cpue = cpue_n_km2)
}

dat <- dat[!is.na(dat$cpue),]

# Add Year-Age interaction
ages <- unique(dat$age)
dat$age_f <- factor(paste0("age_", dat$age))
dat$year_age <- interaction(dat$year, dat$age_f)

# Project data to UTM
dat <- st_as_sf(dat,
                 coords = c('lon','lat'),
                 crs = st_crs(4326))
dat <- st_transform(dat,
                     crs = st_crs("+proj=utm +zone=2 +units=km"))
# Add UTM coordinates as columns X & Y
dat <- cbind(st_drop_geometry(dat), st_coordinates(dat))

# For yellowfin sole and Pacific cod, drop data for levels with all zeros 
if(species %in% c("yellowfin_sole", "pacific_cod")) {
  N_z <- tapply(dat$cpue, 
                INDEX = dat$year_age, 
                FUN = \(x) sum(x > 0))
  year_age_to_drop <- names(which(N_z == 0))
  dat <- subset(dat, !(year_age %in% year_age_to_drop))
  dat <- droplevels(dat)
  
  # Check results
  tapply(dat$cpue,
         INDEX = list(dat$age, dat$year), 
         FUN = \(x) sum(x > 0))
}

# Inputs to tinyVAST ----------------------------------------------------------
# Constant AR1 spatio-temporal term across ages & different variances for each age
dsem <- "\n  "
for(i in min(ages):max(ages)) {
  dsem <- paste0(dsem,  "age_", i, " -> age_", i, ", 1, lag1\n")
}

# # TinyVAST mesh
# mesh <- fm_mesh_2d(loc = Data[,c("X","Y")],
#                    cutoff = 50)
# 
# Format mesh for tinyVAST (using a function from sdmTMB; same method as sdmTMB index bridging)
old_mesh <- sdmTMB::make_mesh(dat, 
                              xy_cols = c("X", "Y"), 
                              mesh = readRDS(file = here("meshes/bs_vast_mesh_50_knots.RDS")),
                              fmesher_func = fm_mesh_2d()) # unsure if this is needed or what it is doing?

control <- tinyVASTcontrol(getsd = FALSE,
                           profile = c("alpha_j"),
                           trace = 0)

family <- setNames(lapply(ages, function(x) delta_gamma(type = "poisson-link")), 
                   paste0("age_", ages))

# Fit model -------------------------------------------------------------------
fit <- tinyVAST(
  data = dat,
  formula = cpue ~ 0 + year_age,  
  sem = "",
  dsem = dsem,
  family = family,
  delta_options = list(delta_formula = ~ 0 + factor(year_age)),  # 2nd linear predictor
  space_column = c("X", "Y"), 
  variable_column = "age_f",
  time_column = "year",
  distribution_column = "age_f",
  spatial_graph = old_mesh,
  control = control
)

# Save fit object
saveRDS(fit, here(workDir, "results_age", paste0("tinyVAST_fit.RDS")))

# Model diagnostics -----------------------------------------------------------
# working off of this vignette: https://vast-lib.github.io/tinyVAST/articles/mgcv.html
sim <- replicate(n = 100, expr = fit$obj$simulate()$y_i)

res <- DHARMa::createDHARMa(simulatedResponse = sim, 
                            observedResponse = dat$cpue, 
                            fittedPredictedResponse = fitted(fit))

res_data <- residuals(res)

# q-q plot
library(ggplot2)
library(qqplotr)

qqplot <- ggplot(data.frame(resid = res_data), aes(sample = resid)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  theme_bw()

ggsave(qqplot, filename = here(workDir, "results_age", "qq.png"),
       width=130, height=130, units="mm", dpi=300)
