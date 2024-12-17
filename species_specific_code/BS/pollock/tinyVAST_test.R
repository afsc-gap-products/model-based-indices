#' Test of the age composition expansion for BS pollock in tinyVAST. Original
#' vignette by Jim Thorson here: 
#' https://vast-lib.github.io/tinyVAST/articles/web_only/age_composition_expansion.html
#' Updates by Sophia Wassermann

# Install tinyVAST
# library(devtools)
# install_github("vast-lib/tinyVAST", dependencies = TRUE)

library(tinyVAST)
library(fmesher)
library(sf)
library(here)

#' To start, we load sampling data that has undergone first-stage expansion. 
#' This arises when each primary sampling unit includes secondary subsampling 
#' of ages, and the subsampled proportion-at-age in each primary unit has been 
#' expanded to the total abundance in that primary sample:
#' ----------------------------------------------------------------------------
species <- 21740
this_year <- lubridate::year(Sys.Date())
# this_year <- 2022  # different year for debugging
Species <- "pollock"
speciesName <- paste0("Walleye_Pollock_age_", as.character(this_year), "_EBS-NBS")
workDir <- here::here("species_specific_code", "BS", Species)
Data <- read.csv(here("species_specific_code", "BS", Species, "data", 
                      paste0("VAST_ddc_alk_", this_year, ".csv")))

# Add Year-Age interaction
Data$Age <- factor(paste0("Age_", Data$Age))
Data$Year_Age <- interaction(Data$Year, Data$Age)

# Project data to UTM
Data <- st_as_sf(Data,
                 coords = c('Lon','Lat'),
                 crs = st_crs(4326))
Data <- st_transform(Data,
                     crs = st_crs("+proj=utm +zone=2 +units=km"))
# Add UTM coordinates as columns X & Y
Data <- cbind(st_drop_geometry(Data), st_coordinates(Data))

#' Next, we construct the various inputs to tinyVAST
#' ----------------------------------------------------------------------------
# adds different variances for each age
sem <- ""

# Constant AR1 spatio-temporal term across ages
# and adds different variances for each age
dsem <- "
  Age_1 -> Age_1, 1, lag1
  Age_2 -> Age_2, 1, lag1
  Age_3 -> Age_3, 1, lag1
  Age_4 -> Age_4, 1, lag1
  Age_5 -> Age_5, 1, lag1
  Age_6 -> Age_6, 1, lag1
  Age_7 -> Age_7, 1, lag1
  Age_8 -> Age_8, 1, lag1
  Age_9 -> Age_9, 1, lag1
  Age_10 -> Age_10, 1, lag1
  Age_11 -> Age_11, 1, lag1
  Age_12 -> Age_12, 1, lag1
  Age_13 -> Age_13, 1, lag1
  Age_14 -> Age_14, 1, lag1
  Age_15 -> Age_15, 1, lag1
"

# TinyVAST mesh
mesh <- fm_mesh_2d(loc = Data[,c("X","Y")],
                   cutoff = 50)

control <- tinyVASTcontrol(getsd = FALSE,
                           profile = c("alpha_j"),  
                           trace = 0)

# SNW: use mesh from the VAST model for bridging ------------------------------
# Read in VAST model & get mesh
VASTfit_age <- readRDS(here("VAST_results", "BS", "pollock", "VASTfit_age.RDS"))
mesh_vast <- VASTfit_age$spatial_list$MeshList$anisotropic_mesh

# Format mesh for tinyVAST (using a function from sdmTMB; same method as sdmTMB index bridging)
old_mesh <- sdmTMB::make_mesh(dat, xy_cols = c("X", "Y"), 
                              mesh = VASTfit_age$spatial_list$MeshList$anisotropic_mesh,
                              fmesher_func = fm_mesh_2d())

#' We could run the model with a log-linked Tweedie distribution and a single 
#' linear predictor:
#' ----------------------------------------------------------------------------
family <- list(
  Age_1 = tweedie(),
  Age_2 = tweedie(),
  Age_3 = tweedie(),
  Age_4 = tweedie(),
  Age_5 = tweedie(), 
  Age_6 = tweedie(),
  Age_7 = tweedie(),
  Age_8 = tweedie(),
  Age_9 = tweedie(),
  Age_10 = tweedie(),
  Age_11 = tweedie(),
  Age_12 = tweedie(),
  Age_13 = tweedie(),
  Age_14 = tweedie(),
  Age_15 = tweedie()     
)

# Data$Year = factor(Data$Year)
start.time <- Sys.time() 
myfit <- tinyVAST(
  data = Data,
  formula = CPUE_num ~ 0 + Year_Age,  
  sem = sem,
  dsem = dsem,
  family = family,
  space_column = c("X", "Y"), 
  variable_column = "Age",
  time_column = "Year",
  distribution_column = "Age",
  spatial_graph = old_mesh,
  control = control
)
stop.time <- Sys.time()

# Save fit object
saveRDS(myfit, here(workDir, "results", "tinyVAST_fit_mesh.RDS"))

#' ----------------------------------------------------------------------------
#' SNW: this step takes longer than the model fitting and is memory-intensive. 
#' It's a good idea to re-start R before this point. You'll need to re-run lines
#' 10-40, too.

# Load back fit object if you restarted R
myfit <- readRDS(here(workDir, "results", "tinyVAST_fit_mesh.RDS"))
# Get shapefile for survey extent
data(bering_sea)

# Make extrapolation grid based on shapefile
bering_sea <- st_transform(bering_sea, 
                           st_crs("+proj=utm +zone=2 +units=km"))
grid <- st_make_grid(bering_sea, n=c(50,50))
grid <- st_intersection(grid, bering_sea)
grid <- st_make_valid(grid)
loc_gz <- st_coordinates(st_centroid(grid))

# Get area for extrapolation grid
library(units)
areas <- set_units(st_area(grid), "hectares") #  / 100^2 # Hectares

# Get abundance
N_jz <- expand.grid(Age = myfit$internal$variables, Year = sort(unique(Data$Year)))
N_jz <- cbind(N_jz, "Biomass" = NA, "SE" = NA)
for(j in seq_len(nrow(N_jz))){
  if(N_jz[j, 'Age'] == 1){
    message("Integrating ", N_jz[j,'Year'], " ", N_jz[j,'Age'], ": ", Sys.time())
  }
  if(is.na(N_jz[j,'Biomass'])){
    newdata = data.frame(loc_gz, Year = N_jz[j, 'Year'], Age = N_jz[j, 'Age'])
    newdata$Year_Age = paste(newdata$Year, newdata$Age, sep=".")
    # Area-expansion
    index1 = integrate_output(myfit,
                              area = areas,
                              newdata = newdata,
                              apply.epsilon = TRUE,
                              bias.correct = TRUE,
                              intern = TRUE )
    N_jz[j, 'Biomass'] = index1[3] / 1e9
  }
}
N_ct <- array(N_jz$Biomass, dim=c(length(myfit$internal$variables),length(unique(Data$Year))),
              dimnames=list(myfit$internal$variables,sort(unique(Data$Year))))
N_ct <- N_ct / outer(rep(1, nrow(N_ct)), colSums(N_ct))

#' Finally, we can compare these estimates with those from package VAST. 
#' Estimates differ somewhat because VAST used a delta-gamma distribution with 
#' spatio-temporal variation in two linear predictors, and also used a different mesh.
#' ----------------------------------------------------------------------------
# Load VAST results for same data
data(bering_sea_pollock_vast)
myvast <- bering_sea_pollock_vast
rownames(myvast) <- 1:15

# Reformat tinyVAST output with same dimnames
mytiny <- N_ct
rownames(mytiny) <- 1:15

longvast <- cbind(expand.grid(dimnames(myvast)), 
                  p = as.numeric(myvast),
                  method = "VAST")
longtiny <- cbind(expand.grid(dimnames(mytiny)), 
                  p = as.numeric(mytiny), 
                  method = "tinyVAST")
long <- rbind(longvast, longtiny)

library(ggplot2)
ggplot(data = long, aes(x = Var2, y = p, col = method)) +
  facet_grid(rows = vars(Var1), scales = "free") +
  geom_point() +
  scale_y_log10()

# SNW: Reformat and save a version for comparison plots -----------------------
tiny_out <- tibble::rownames_to_column(data.frame(t(mytiny)), "VALUE")
tiny_out <- tiny_out[, c(2:16, 1)]
colnames(tiny_out) <- c("age_1", "age_2", "age_3", "age_4", "age_5", "age_6",
                        "age_7", "age_8", "age_9", "age_10", "age_11", "age_12",
                        "age_13", "age_14", "age_15", "Year")
tiny_out$Region <- "EBS"
write.csv(tiny_out, here(workDir, "results", "tinyVAST_props.csv"), row.names = FALSE)
