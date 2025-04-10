# Plot all GOA sdmTMB indices compared to design-based indices ----
# Author: Lewis Barnett
# Date: 12-3-2025

library(dplyr)
library(ggplot2)
library(here)

phase <- c("hindcast", "production")[1] # specify analysis phase

channel <- gapindex::get_connected(check_access = FALSE) # enter credentials 

species_df <- data.frame(species_name = c("Gadus_macrocephalus", 
                                          "Sebastes_alutus", 
                                          "Sebastes_polyspinis", 
                                          "Squalus_suckleyi", 
                                          "Atheresthes_stomias", 
                                          "Sebastes_variabilis"), 
                         species_code = c("21720",
                                          "30060",
                                          "30420",
                                          "310",
                                          "10110",
                                          "30152")
                         )

for (i in 1:nrow(species_df)){
  species <- species_df$species_name[i]
  
  dir <- here("species_specific_code", "GOA", species, phase, "/")
  if (!dir.exists(paths = dir)) dir.create(path = dir, recursive = TRUE)

  # query db index
  if(species == "Gadus_macrocephalus"){
    query <- paste0('
          SELECT 
          \'db\' AS "index",
          YEAR AS "year",
          POPULATION_COUNT AS "est",
          (POPULATION_COUNT - SQRT(POPULATION_VAR) * (1.959964)) as "lwr", -- qnorm(0.025) in R
          (POPULATION_COUNT + SQRT(POPULATION_VAR) * (1.959964)) as "upr" -- qnorm(0.975) in R
          
          -- Identify what tables to pull data from
          FROM GAP_PRODUCTS.BIOMASS
          WHERE AREA_ID = 99903 -- GOA REGION
          AND SPECIES_CODE = ', species_df$species_code[i]
    )
  } else {
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
          AND SPECIES_CODE = ', species_df$species_code[i]
    )
  }
  db_i <- RODBC::sqlQuery(channel = channel, query = query)

  # load latest model-based index estimates  
  new_i <- readRDS(here("species_specific_code", "GOA", 
                         species, phase, "index.RDS")) %>% 
    mutate(index = "mb") %>% 
    select(index, year, est, lwr, upr)
  both_i <- bind_rows(new_i, db_i) %>% 
    filter(est > 0)
  both_i[both_i < 0] <- 0
  
  if(species == "Gadus_macrocephalus"){
    ylab <- "Abundance (n)"
  } else {
    ylab <- "Biomass (kg)"
  }
  ggplot(both_i, aes(x = year, y = est, ymin = lwr, ymax = upr, 
                     colour = index)) + 
    geom_ribbon(alpha = 0.1) +
    geom_line(alpha = 0.8) + 
    ylim(0, max(both_i$upr)) +
    ggtitle(species) +
    coord_cartesian(expand = FALSE) + 
    ylab(ylab) +
    theme_bw()
  ggsave(file = here("species_specific_code", "GOA", species, phase, 
                     "index_comparison_db.pdf"), 
         height = 4, width = 6, units = c("in"))
}