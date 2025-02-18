# Generate VAST-type meshes with lower knots for experimental runs

library(VAST)
library(sp)
library(sdmTMB)
library(dplyr)
library(ggplot2)
library(here)

phase <- c("hindcast", "production")[1] # specify analysis phase
species <- "kamchatka_flounder"

dat <- readRDS(here("species_specific_code", "BS", species, phase, 
                      "data", "data_geostat_index.RDS"))
dat <- rename(dat, cpue = cpue_n_km2)
dat <- dat[!is.na(dat$cpue),]
dat$year_f <- as.factor(dat$year)
dat <- add_utm_columns(dat, ll_names = c("lon", "lat"), utm_crs = 32602, units = "km")

FieldConfig <- c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig <- c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 4, "Epsilon2" = 4)
OverdispersionConfig <- c("Eta1"=0, "Eta2"=0)

for(knots in seq(250, 500, 50)){
  settings <- FishStatsUtils::make_settings(
    n_x = knots, # number of vertices in the SPDE mesh
    Region = "Eastern_Bering_Sea",
    purpose = "index2", # index of abundance with Gamma for positive catches
    fine_scale = TRUE, # use bilinear interpolation from the INLA 'A' matrix
    #zone = 2,
    FieldConfig = FieldConfig,
    RhoConfig = RhoConfig,
    OverdispersionConfig = c("Eta1" = 0, "Eta2" = 0),
    ObsModel = c(10, 2), #delta-Gamma; (2,4) if there are years with 100% encounter rate; (10, 2) for Tweedie
    use_anisotropy = TRUE,
    max_cells = Inf, # use all grid cells from the extrapolation grid, production model used 2000
    knot_method = "grid"
  )
  
  fit <- fit_model(
    settings = settings,
    Lat_i = dat[, "lat"],
    Lon_i = dat[, "lon"],
    t_i = dat[, "year"],
    b_i = dat[, "cpue"],
    a_i = rep(1, nrow(dat)),
    create_strata_per_region = TRUE,
    run_model = FALSE,
    working_dir = paste0(here("species_specific_code", "BS", species, "index_comparison"), "/")
  )
  
  mesh <-  sdmTMB::make_mesh(dat, xy_cols = c("X", "Y"), mesh = fit$spatial_list$MeshList$anisotropic_mesh) 
  plot(mesh)
  
  saveRDS(fit$spatial_list$MeshList$anisotropic_mesh, file = here("meshes", paste0("ebs_vast_mesh_", knots, "_knots.RDS")))
}