library(dplyr)
library(ggplot2)
library(RODBC)

source("R/get_connected.R")

db <- RODBC::sqlQuery(
  channel = channel, 
  query = paste("SELECT YEAR, SPECIES_CODE, TOTAL_BIOMASS, BIOMASS_VAR",
                "FROM GOA.BIOMASS_TOTAL WHERE",
                "SPECIES_CODE = 30060 AND",
                "YEAR >= 1990"))
names(db) <- tolower(names(db))
db$sd_mt <- sqrt(db$biomass_var)
db$method <- "design-based"
db <- db %>% rename(estimate_metric_tons = total_biomass) %>% 
  select(year, estimate_metric_tons, sd_mt, method)

mb <- read.csv(file = paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/results/Table_for_SS3.csv"))
names(mb) <- tolower(names(mb))
mb$method <- "model-based"
mb <- select(mb, year, estimate_metric_tons, sd_mt, method)
mb <- mb[complete.cases(mb),]

dat <- bind_rows(mb, db)

ggplot(dat, aes(x = year, y = estimate_metric_tons, group = method, color = method)) +
  geom_pointrange(aes(ymin=estimate_metric_tons-sd_mt, ymax=estimate_metric_tons+sd_mt),
                  position = position_dodge(1.1))

ggsave(paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/results/GOA_POP_compare_DBE.png"), 
       width = 6,
       height = 3,
       units = "in")