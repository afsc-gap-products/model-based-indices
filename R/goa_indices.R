## Scratch/draft script to run all GOA sdmTMB indices ----

#remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE)
library(sdmTMB)
library(dplyr)
library(ggplot2)
library(here)

phase <- c("hindcast", "production")[1] # specify analysis phase

channel <- gapindex::get_connected(check_access = FALSE) # enter credentials in pop-out window

species_list <- c("Gadus_macrocephalus", 
                  "Sebastes_alutus", "Sebastes_polyspinis", 
                  "Sebastes_variabilis", "Squalus_suckleyi",
                  "Atheresthes_stomias")

for (i in species_list){
  species <- i
  
  dir <- here("species_specific_code", "GOA", species, phase, "/")
  if (!dir.exists(paths = dir)) dir.create(path = dir, recursive = TRUE)
  
  # load and process data ----
  sp <- gsub("_", " ", species)
  dat <- readRDS(here(paste0("data/GOA/", phase, "/dat_allspp.RDS"))) %>%
    filter(species == sp)
  dat$year_f <- as.factor(dat$year)
  
  # fit model ----
  f1 <- here("species_specific_code", "GOA", species, 
             "index_comparison", "fit_sdmTMB.RDS")
  if (!file.exists(f1)) {
    #pass same mesh as prior model
    mesh <-  make_mesh(dat, xy_cols = c("X", "Y"), 
                       mesh = readRDS(file = here("meshes/goa_vast_mesh.RDS")))
    #or use coarser mesh for experimentation
    #mesh <-  make_mesh(dat, xy_cols = c("X", "Y"), n_knots = 50) 
    
    fit_sdmTMB <- sdmTMB( 
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
      #, do_index = TRUE (to compute index at same time, requires passing args)
    )
    fit_sdmTMB
    saveRDS(fit_sdmTMB, file = f1)
  } else {
    fit_sdmTMB <- readRDS(f1)
  }
  
  # load grid and process prediction grid for all years desired
  grid <- readRDS(file = "extrapolation_grids/goa_sdmtmb_grid.RDS")
  pred_grid <- replicate_df(grid, "year_f", unique(dat$year_f))
  pred_grid$year <- as.integer(as.character(factor(pred_grid$year_f)))
  
  # make predictions ----
  f2 <- here("species_specific_code", "GOA", species, 
             "index_comparison", "predictions.RDS")
  if (!file.exists(f2)) {
    p <- predict(fit_sdmTMB, newdata = pred_grid, return_tmb_object = TRUE)
    saveRDS(p, file = f2)
  } else {
    p <- readRDS(f2)
  }
  
  # Plot predicted density maps and fit diagnostics ----
  # q-q plot
  pdf(file = here("species_specific_code", "GOA", species, phase, "qq.pdf"), 
      width = 5, height = 5)
    resids <- residuals(fit_sdmTMB, type ="mle-mvn") 
    qqnorm(resids);abline(0, 1)
  dev.off()
  
  # residuals on map plot, by year
  dat$resids <- resids
  ggplot(dat, aes(X, Y, col = resids)) + 
    scale_colour_gradient2() +
    geom_point() + 
    facet_wrap(~year) + 
    coord_fixed()
  ggsave(file = here("species_specific_code", "GOA", species, phase, 
                     "residuals_map.pdf"), 
         height = 9, width = 6, units = c("in"))
  
  # predictions on map plot, by year
  ggplot(p$data, aes(X, Y, fill = exp(est))) +
    geom_raster() +
    scale_fill_viridis_c(trans = "sqrt") +
    facet_wrap(~year) +
    coord_fixed() +
    ggtitle("Predicted densities")
  ggsave(file = here("species_specific_code", "GOA", species, phase, 
                     "predictions_map.pdf"), 
         height = 9, width = 6, units = c("in"))
  
  # compute index ----
  f3 <- here("species_specific_code", "GOA", species, 
             "index_comparison", "index.RDS")
  if (!file.exists(f3)) {
    ind <- get_index(p, bias_correct = TRUE, area = p$data$area_km2)
    saveRDS(ind, file = f3)
  } else {
    ind <- readRDS(f3)
  }
  
  # compare indices plot (change loading of old index after hindcast) ----
  # query db index
  query <- paste0('
  SELECT 
  \'db\' AS "index",
  YEAR AS "year",
  BIOMASS_MT * 1000 AS "est",
  (BIOMASS_MT - SQRT(BIOMASS_VAR) * (1.959964)) * 1000 as "lwr", -- qnorm(0.025) in R
  (BIOMASS_MT + SQRT(BIOMASS_VAR) * (1.959964)) * 1000 as "upr" -- qnorm(0.975) in R
  
  -- Identify what tables to pull data from
  FROM GAP_PRODUCTS.BIOMASS
  WHERE AREA_ID = 99903 -- GOA REGION
  AND YEAR < 2025 -- REMOVE THIS LINE AFTER 2025 GOA MOCK DATA HAVE BEEN REMOVED
  AND SPECIES_CODE = ', unique(dat$species_code)
  )
  db_i <- RODBC::sqlQuery(channel = channel, query = query)
  
  new_i <- ind %>% mutate(index = "new") %>% select(index, year, est, lwr, upr)
  old_i <- read.csv(here("species_specific_code", "GOA", species, 
                         "index_comparison", "Index.csv")) %>%
    mutate(index = "old", year = as.numeric(Time), est = Estimate, 
           se = Std..Error.for.ln.Estimate.) %>% 
    filter(year %in% unique(new_i$year)) %>%
    mutate(lwr = exp(log(est) + qnorm(0.025) * se),
           upr = exp(log(est) + qnorm(0.975) * se)) %>%
    select(index, year, est, lwr, upr)
  both_i <- bind_rows(new_i, old_i, db_i) %>% 
    filter(est > 0)
  both_i[both_i<0] <- 0
    
  
  ggplot(both_i, aes(x = year, y = est, ymin = lwr, ymax = upr, 
                     colour = index)) + 
    geom_ribbon(alpha = 0.1) +
    geom_line(alpha = 0.8) + 
    ylim(0, max(both_i$upr)) +
    ggtitle(species) +
    coord_cartesian(expand = FALSE) + 
    ylab("Biomass (kg)") +
    theme_bw()
  ggsave(file = here("species_specific_code", "GOA", species, phase, 
                     "index_comparison.pdf"), 
         height = 4, width = 6, units = c("in"))
  
  # ESP products ----
  # TODO: check if there is different trend without bias correction
  cog <- get_cog(p, bias_correct = FALSE, 
                 area = p$data$area_km2, format = "wide")
  saveRDS(cog, file = here("species_specific_code", "GOA", 
                           species, phase, "cog.RDS"))
  
  eao <- get_eao(p, bias_correct = FALSE, area = p$data$area_km2)
  saveRDS(eao, file = here("species_specific_code", "GOA", 
                           species, phase, "area_occupied.RDS"))
}