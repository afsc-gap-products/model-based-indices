#' Test of an age composition model using sdmTMB, following a vignette for 
#' multispecies models: https://github.com/pbs-assess/sdmTMB/blob/main/vignettes/web_only/multispecies.Rmd
#' and a population composition example: https://github.com/pbs-assess/sdmTMB/blob/comps/scratch/composition-example.R
# By: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.12.09

library(here)
library(sdmTMB)
library(dplyr)
library(sp)

# Set species -----------------------------------------------------------------
species <- 21740
this_year <- lubridate::year(Sys.Date())
# this_year <- 2022  # different year for debugging
Species <- "pollock"
speciesName <- paste0("Walleye_Pollock_age_", as.character(this_year), "_EBS-NBS")
workDir <- here::here("species_specific_code", "BS", Species)
Data <- read.csv(here("species_specific_code", "BS", Species, "data", 
                      paste0("VAST_ddc_alk_", this_year, ".csv")))

# Set up sdmTMB model ---------------------------------------------------------
# Re-format data for sdmTMB
dat <-  transmute(Data,
                  cpue = CPUE_num, # converts cpue from kg/ha to kg/km^2
                  year = as.integer(Year),
                  age = as.integer(Age),
                  vessel = "missing",
                  effort = 1, # area swept is 1 when using CPUE instead of observed weight
                  Y = Lat,
                  X = Lon,
                  pass = 0) %>% 
  as.data.frame() %>%  # ensure not a tibble
  filter(year != 2020)  # drop dummy 2020 data

dat$year_f <- as.factor(dat$year)

# Transform lat & lon to UTM
coordinates(dat) <- ~ X + Y
proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")
dat <- as.data.frame(spTransform(dat, CRS("+proj=utm +zone=2")))
# scale to km so values don't get too large
dat$X <- dat$coords.x1 / 1000
dat$Y <- dat$coords.x2 / 1000

# Read in VAST fit object to use mesh
VASTfit <- readRDS(here("VAST_results", "BS", "pollock", "VASTfit_age.RDS"))
mesh <- make_mesh(dat, xy_cols = c("X", "Y"), mesh = VASTfit$spatial_list$MeshList$anisotropic_mesh)

# Fit sdmTMB model ------------------------------------------------------------
fit_sdmTMB <- sdmTMB( 
  cpue ~ year_f * age,
  spatial_varying = ~ 0 + factor(age),
  data = dat, 
  mesh = mesh,
  family = delta_gamma(type = "poisson-link"), 
  time = "year", 
  spatial = "off",
  spatiotemporal = "iid",
  extra_time = 2020L, #  omit if dummy 2020 included in data
  silent = FALSE,
  anisotropy = TRUE,
  do_fit = TRUE
)
fit_sdmTMB
saveRDS(fit_sdmTMB, file = here("species_specific_code", "BS", Species, "results", "fit_sdmTMB_age.RDS"))
