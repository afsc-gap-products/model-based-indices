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

# Load data -------------------------------------------------------------------
species <- 21740
# this_year <- lubridate::year(Sys.Date())
this_year <- 2024  # different year for debugging
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


# Inputs to tinyVAST ----------------------------------------------------------
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

# # TinyVAST mesh
# mesh <- fm_mesh_2d(loc = Data[,c("X","Y")],
#                    cutoff = 50)
# 
control <- tinyVASTcontrol(getsd = FALSE,
                           profile = c("alpha_j"),
                           trace = 0)

# Use mesh from the VAST model for bridging -----------------------------------
# Format mesh for tinyVAST (using a function from sdmTMB; same method as sdmTMB index bridging)
old_mesh <- sdmTMB::make_mesh(Data, 
                              xy_cols = c("X", "Y"), 
                              mesh = readRDS(file = here("meshes/bs_vast_mesh_50_knots.RDS")),
                              fmesher_func = fm_mesh_2d()) # unsure if this is needed or what it is doing?

#' Run the model with a Poisson-linked delta gamma distribution with time- and 
#' age- varying intercepts in the second linear predictor
#' ----------------------------------------------------------------------------
family <- list(
  Age_1 = delta_gamma(type = "poisson-link"),
  Age_2 = delta_gamma(type = "poisson-link"),
  Age_3 = delta_gamma(type = "poisson-link"),
  Age_4 = delta_gamma(type = "poisson-link"),
  Age_5 = delta_gamma(type = "poisson-link"), 
  Age_6 = delta_gamma(type = "poisson-link"),
  Age_7 = delta_gamma(type = "poisson-link"),
  Age_8 = delta_gamma(type = "poisson-link"),
  Age_9 = delta_gamma(type = "poisson-link"),
  Age_10 = delta_gamma(type = "poisson-link"),
  Age_11 = delta_gamma(type = "poisson-link"),
  Age_12 = delta_gamma(type = "poisson-link"),
  Age_13 = delta_gamma(type = "poisson-link"),
  Age_14 = delta_gamma(type = "poisson-link"),
  Age_15 = delta_gamma(type = "poisson-link")     
)

start.time <- Sys.time() 
myfit <- tinyVAST(
  data = Data,
  formula = CPUE_num ~ 0 + Year_Age,  
  sem = sem,
  dsem = dsem,
  family = family,
  delta_options = list(delta_formula = ~ 0 + factor(Year_Age)),  # 2nd linear predictor
  space_column = c("X", "Y"), 
  variable_column = "Age",
  time_column = "Year",
  distribution_column = "Age",
  spatial_graph = old_mesh,
  control = control
)
stop.time <- Sys.time()

# Save fit object
saveRDS(myfit, here(workDir, "results", paste0("tinyVAST_fit_", as.character(this_year), ".RDS")))


# Index expansion -------------------------------------------------------------
# Load fit object if needed
myfit <- readRDS(here(workDir, "results", paste0("tinyVAST_fit_", as.character(this_year), ".RDS")))


# Read in coarsened extrapolation grid for desired region 
# coarse_grid <- read.csv(here("extrapolation_grids", "ebs_coarse_grid.csv"))  # EBS
# coarse_grid <- read.csv(here("extrapolation_grids", "nbs_coarse_grid.csv"))  # NBS
coarse_grid <- read.csv(here("extrapolation_grids", "bering_coarse_grid.csv"))  # both regions

# Get abundance
N_jz <- expand.grid(Age = myfit$internal$variables, Year = sort(unique(Data$Year)))
N_jz <- cbind(N_jz, "Biomass" = NA, "SE" = NA)
for(j in seq_len(nrow(N_jz))){
  if(N_jz[j, 'Age'] == 1){
    message("Integrating ", N_jz[j,'Year'], " ", N_jz[j,'Age'], ": ", Sys.time())
  }
  if(is.na(N_jz[j,'Biomass'])){
    newdata = data.frame(coarse_grid, Year = N_jz[j, 'Year'], Age = N_jz[j, 'Age'])
    newdata$Year_Age = paste(newdata$Year, newdata$Age, sep=".")
    # Area-expansion
    index1 = integrate_output(myfit,
                              area = coarse_grid$area_km2,
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

# Save abundance estimate
write.csv(N_ct, here(workDir, "results", paste0("tinyVAST_natage_", as.character(this_year), ".csv")), row.names = FALSE)


# Reformat and calculate proportions ------------------------------------------
rownames(N_ct) <- 1:15
tiny_out <- tibble::rownames_to_column(data.frame(t(N_ct)), "VALUE")
tiny_out <- tiny_out[, c(2:16, 1)]
colnames(tiny_out) <- c("age_1", "age_2", "age_3", "age_4", "age_5", "age_6",
                        "age_7", "age_8", "age_9", "age_10", "age_11", "age_12",
                        "age_13", "age_14", "age_15", "Year")

# Save proportions
write.csv(tiny_out, here(workDir, "results", paste0("tinyVAST_props_", as.character(this_year), ".csv")), row.names = FALSE)


# Diagnostic plots ------------------------------------------------------------
# working off of this vignette: https://vast-lib.github.io/tinyVAST/articles/mgcv.html
sim <- replicate(n = 500, expr = myfit$obj$simulate()$y_i)

res <- DHARMa::createDHARMa(simulatedResponse = sim, 
                            observedResponse = Data$CPUE_num, 
                            fittedPredictedResponse = fitted(myfit))

res_data <- residuals(res)

# q-q plot
library(ggplot2)
library(qqplotr)

ggplot(data.frame(resid = res_data), aes(sample = resid)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  theme_bw()
