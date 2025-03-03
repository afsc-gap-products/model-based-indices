#' Test of the age composition expansion for BS cod in tinyVAST. Original
#' vignette for pollock by Jim Thorson here: 
#' https://vast-lib.github.io/tinyVAST/articles/web_only/age_composition_expansion.html
#' Updates by Sophia Wassermann and Lewis Barnett

# Install tinyVAST
# library(devtools)
# install_github("vast-lib/tinyVAST", dependencies = TRUE)

library(tinyVAST)
library(fmesher)
library(sf)
library(here)

# Load data -------------------------------------------------------------------
species <- 21720
# this_year <- lubridate::year(Sys.Date())
this_year <- 2024  # different year for debugging
Species <- "pacific_cod"
speciesName <- paste0("pacific_cod_age_", as.character(this_year), "_EBS-NBS")
workDir <- here::here("species_specific_code", "BS", Species)
Data <- readRDS(here::here("species_specific_code", "BS", Species, "hindcast", "data", 
                   "data_geostat_agecomps.RDS"))

# Add year-age interaction
Data$age <- factor(paste0("age_", Data$age))
Data$year_age <- interaction(Data$year, Data$age)

# Project data to UTM
Data <- st_as_sf(Data,
                 coords = c('lon','lat'),
                 crs = st_crs(4326))
Data <- st_transform(Data,
                     crs = st_crs("+proj=utm +zone=2 +units=km"))
# Add UTM coordinates as columns X & Y
Data <- cbind(st_drop_geometry(Data), st_coordinates(Data))

#  DROP DATA FOR LEVELS WITH ALL ZEROS
N_z = tapply( Data$cpue_n_km, INDEX = Data$year_age, FUN = \(x) sum(x>0) )
year_age_to_drop = names(which( N_z==0 ))
Data = subset( Data, !(year_age %in% year_age_to_drop) )
Data = droplevels(Data)

# Check results
tapply( Data$cpue_n_km,
        INDEX = list( Data$age,
                      Data$year),
        FUN = \(x) sum(x>0) )

# Npool = 100 in VAST
N_a = tapply( Data$cpue_n_km, INDEX = Data$age, FUN=\(x) sum(x>0) )
which(N_a < 100)

# Inputs to tinyVAST ----------------------------------------------------------
sem <- ""

# Constant AR1 spatio-temporal term across ages
# and adds different variances for each age
dsem <- "
  age_0 -> age_0, 1, lag1
  age_1 -> age_1, 1, lag1
  age_2 -> age_2, 1, lag1
  age_3 -> age_3, 1, lag1
  age_4 -> age_4, 1, lag1
  age_5 -> age_5, 1, lag1
  age_6 -> age_6, 1, lag1
  age_7 -> age_7, 1, lag1
  age_8 -> age_8, 1, lag1
  age_9 -> age_9, 1, lag1
  age_10 -> age_10, 1, lag1
  age_11 -> age_11, 1, lag1
  age_12 -> age_12, 1, lag1
"

# # TinyVAST mesh
# mesh <- fm_mesh_2d(loc = Data[,c("X","Y")],
#                    cutoff = 50)
# 
# Use mesh from the VAST model for bridging -----------------------------------
# Format mesh for tinyVAST (using a function from sdmTMB; same method as sdmTMB index bridging)
old_mesh <- sdmTMB::make_mesh(Data, 
                              xy_cols = c("X", "Y"), 
                              mesh = readRDS(file = "meshes/bs_vast_mesh_50_knots.RDS"))

#' Run the model with a Poisson-linked delta gamma distribution with time- and 
#' age- varying intercepts in the second linear predictor
#' ----------------------------------------------------------------------------
family <- list(
  age_0 = delta_gamma(type = "poisson-link"),
  age_1 = delta_gamma(type = "poisson-link"),
  age_2 = delta_gamma(type = "poisson-link"),
  age_3 = delta_gamma(type = "poisson-link"),
  age_4 = delta_gamma(type = "poisson-link"),
  age_5 = delta_gamma(type = "poisson-link"), 
  age_6 = delta_gamma(type = "poisson-link"),
  age_7 = delta_gamma(type = "poisson-link"),
  age_8 = delta_gamma(type = "poisson-link"),
  age_9 = delta_gamma(type = "poisson-link"),
  age_10 = delta_gamma(type = "poisson-link"),
  age_11 = delta_gamma(type = "poisson-link"),
  age_12 = delta_gamma(type = "poisson-link")  
)

control <- tinyVASTcontrol(profile = c("alpha_j","alpha2_j"),
                           trace = 1,
                           #nlminb_loops = 0,
                           #newton_loops = 0,
                           getsd = FALSE,
                           calculate_deviance_explained = FALSE)

start.time <- Sys.time() 
myfit <- tinyVAST(
  data = Data,
  formula = cpue_n_km2 ~ 0 + year_age,  
  sem = sem,
  dsem = dsem,
  family = family,
  delta_options = list(delta_formula = ~ 0 + factor(year_age)),  # 2nd linear predictor
  space_column = c("X", "Y"), 
  variable_column = "age",
  time_column = "year",
  distribution_column = "age",
  spatial_graph = old_mesh,
  control = control
)
stop.time <- Sys.time()

if( FALSE ){
  myfit$obj$fn( myfit$obj$par )
  Gr = myfit$obj$gr( myfit$obj$par )
  
  #
  P_ct = tapply( Data$cpue_n_km,
                 INDEX = list( Data$age,
                               Data$year),
                 FUN = \(x) mean(x>0) )
  P_z = tapply( Data$cpue_n_km, INDEX = Data$year_age, FUN = \(x) mean(x>0) )
  
  #
  
}

# Save fit object
saveRDS(myfit, here(workDir, "hindcast", "results_age", paste0("tinyVAST_fit_", as.character(this_year), ".RDS")))


# Index expansion -------------------------------------------------------------
# Load fit object if needed
#myfit <- readRDS(here(workDir, "hindcast", "results_age",  paste0("tinyVAST_fit_", as.character(this_year), ".RDS")))


# Read in coarsened extrapolation grid for desired region 
# coarse_grid <- read.csv(here("extrapolation_grids", "ebs_coarse_grid.csv"))  # EBS
# coarse_grid <- read.csv(here("extrapolation_grids", "nbs_coarse_grid.csv"))  # NBS
coarse_grid <- read.csv(here("extrapolation_grids", "bering_coarse_grid.csv"))  # both regions

# Get abundance after dropping age:year combinations with all 0s (MAY NEED TO IMPUTE 0s for these levels after?)
N_jz <- expand.grid(age = myfit$internal$variables, year = sort(unique(Data$year)))
#  DROP LEVELS WITH ALL ZEROS
N_jz$age_ <- factor(paste0("age_", N_jz$age))
N_jz$year_age <- interaction(N_jz$year, N_jz$age)
N_jz <- subset( N_jz, !(year_age %in% year_age_to_drop) )
N_jz <- dplyr::select(N_jz, age, year)
N_jz <- cbind(N_jz, "Biomass" = NA, "SE" = NA)

for(j in seq_len(nrow(N_jz))){
  if(N_jz[j, 'age'] == 1){
    message("Integrating ", N_jz[j,'year'], " ", N_jz[j,'age'], ": ", Sys.time())
  }
  if(is.na(N_jz[j,'Biomass'])){
    newdata = data.frame(coarse_grid, year = N_jz[j, 'year'], age = N_jz[j, 'age'])
    newdata$year_age = paste(newdata$year, newdata$age, sep=".")
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
N_ct <- array(N_jz$Biomass, dim=c(length(myfit$internal$variables),length(unique(Data$year))),
              dimnames=list(myfit$internal$variables,sort(unique(Data$year))))
N_ct <- N_ct / outer(rep(1, nrow(N_ct)), colSums(N_ct))

# Save abundance estimate
write.csv(N_ct, here(workDir, "hindcast", "results_age", paste0("tinyVAST_natage_", as.character(this_year), ".csv")), row.names = FALSE)


# Reformat and calculate proportions ------------------------------------------
rownames(N_ct) <- 1:13
tiny_out <- tibble::rownames_to_column(data.frame(t(N_ct)), "VALUE")
tiny_out <- tiny_out[, c(2:14, 1)]
colnames(tiny_out) <- c("age_0", "age_1", "age_2", "age_3", "age_4", "age_5", "age_6",
                        "age_7", "age_8", "age_9", "age_10", "age_11", "age_12", 
                        "year")

# Save proportions
write.csv(tiny_out, here(workDir, "hindcast", "results_age", paste0("tinyVAST_props_", as.character(this_year), ".csv")), row.names = FALSE)
