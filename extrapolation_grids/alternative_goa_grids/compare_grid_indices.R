## Scratch/draft script to run all GOA sdmTMB indices ----

#remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE)
library(sdmTMB)
library(dplyr)
library(ggplot2)
library(here)
library(gapindex)

## List of species to loop over
species_list <- 
  c("Gadus_macrocephalus", "Sebastes_alutus", "Sebastes_polyspinis", 
    "Sebastes_variabilis", "Squalus_suckleyi", "Atheresthes_stomias")

for (species in species_list) { ## Loop over species -- start
  
  dir <- here("species_specific_code", "GOA", species, "/grid_comparison/")
  if (!dir.exists(paths = dir)) dir.create(path = dir, recursive = TRUE)
  
  # load and process data ----
  sp <- gsub(pattern = "_", replacement = " ", x = species)
  dat <- readRDS(file = here(paste0("data/GOA/hindcast/dat_allspp.RDS"))) %>%
    filter(species == sp)
  dat$year_f <- as.factor(x = dat$year)
  
  # fit model ----
  f1 <- here("species_specific_code", "GOA", species, 
             "grid_comparison", "fit_sdmTMB.RDS")
  if (!file.exists(f1)) {
    #pass same mesh as prior model
    mesh <- sdmTMB::make_mesh(
      data = dat, 
      xy_cols = c("X", "Y"), 
      mesh = readRDS(file = here("meshes/goa_vast_mesh.RDS"))
    )
    
    fit_sdmTMB <- sdmTMB::sdmTMB( 
      cpue_kg_km2 ~ 0 + year_f,
      data = dat, 
      mesh = mesh,
      family = delta_gamma(type = "poisson-link"), 
      time = "year", 
      spatial = "on",
      spatiotemporal = "iid",
      silent = FALSE,
      anisotropy = TRUE,
      do_fit = TRUE
    )
    fit_sdmTMB
    saveRDS(object = fit_sdmTMB, file = f1)
  } else {
    fit_sdmTMB <- readRDS(file = f1)
  }
  
  # load grid and process prediction grid for all years desired
  grid <- readRDS(file = "extrapolation_grids/goa_sdmtmb_grid.RDS")
  pred_grid <- sdmTMB::replicate_df(dat = grid, 
                                    time_name = "year_f", 
                                    time_values = unique(x = dat$year_f))
  pred_grid$year <- 
    as.integer(x = as.character(x = factor(x = pred_grid$year_f)))
  
  # load new grid and process prediction grid for all years desired
  grid_new <- read.csv(file = "extrapolation_grids/alternative_goa_grids/goa_2025_interpolation_grid.csv")
  pred_grid_new <- sdmTMB::replicate_df(dat = as.data.frame(grid_new), 
                                        time_name = "year_f",
                                        time_values = unique(x = dat$year_f))
  pred_grid_new$year <- 
    as.integer(x = as.character(x = factor(x = pred_grid_new$year_f)))
  
  # make predictions, calculate index, save
  p <- predict(fit_sdmTMB, 
               newdata = pred_grid, 
               return_tmb_object = TRUE)
  ind <- sdm get_index(p, 
                       bias_correct = TRUE, 
                       area = p$data$area_km2)
  
  saveRDS(p, file = here("species_specific_code", "GOA", species, 
                         "grid_comparison", "predictions.RDS"))
  saveRDS(ind, file = here("species_specific_code", "GOA", species, 
                           "grid_comparison", "index.RDS"))
  
  # Make predictions on new interpolation grid, calculate index, save
  p_new <- predict(fit_sdmTMB, 
                   newdata = pred_grid_new, 
                   return_tmb_object = TRUE)
  ind_new <- get_index(p_new, 
                       bias_correct = TRUE, 
                       area = p_new$data$area_km2)
  
  saveRDS(object = p_new, file = here("species_specific_code", "GOA", species, 
                                      "grid_comparison", "new_predictions.RDS"))
  saveRDS(ind_new, file = here("species_specific_code", "GOA", species, 
                               "grid_comparison", "index_new.RDS"))
  
  # rbind indices
  both_i <- bind_rows(ind_new %>% mutate(index = "new") %>% 
                        select(index, year, est, lwr, upr), 
                      ind %>% mutate(index = "old") %>% 
                        select(index, year, est, lwr, upr)) %>% 
    filter(est > 0)
  both_i[both_i < 0] <- 0
  
  ## Save Plot
  ggplot(both_i, aes(x = year, y = est, ymin = lwr, ymax = upr, 
                     colour = index)) + 
    geom_ribbon(alpha = 0.1) +
    geom_line(alpha = 0.8) + 
    ylim(0, max(both_i$upr)) +
    ggtitle(species) +
    coord_cartesian(expand = FALSE) + 
    ylab("Biomass (kg)") +
    theme_bw()
  
  ggsave(file = here("species_specific_code", "GOA", species, 
                     "grid_comparison", "grid_comparison.pdf"), 
         height = 4, width = 6, units = c("in"))
  
} ## Loop over species -- end
