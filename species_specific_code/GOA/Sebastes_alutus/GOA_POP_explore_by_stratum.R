library(dplyr)
library(rgdal)
library(rgeos)
library(sf)
library(ggplot2)

dat <- readRDS(paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/data/Data_Geostat_Sebastes_alutus_stratum.rds"))

##################################################
####  Import Data
####  AK_land: alaska land polygons
####  goa_strata: Current GoA strata polygons
##################################################
AK_land <- rgdal::readOGR("C:/Users/lewis.barnett/Work/AFSC/Rprojects/Optimal_Allocation_GoA/data/shapefiles/AKland.shp")
# AK_land <- sp::spTransform(x = AK_land,
#                            CRSobj = "+proj=utm +units=km +zone=5")

goa_strata <- rgdal::readOGR("C:/Users/lewis.barnett/Work/AFSC/Rprojects/Optimal_Allocation_GoA/data/shapefiles/goa_strata.shp")
names(goa_strata@data) <- tolower(names(goa_strata@data))
# goa_strata <- sp::spTransform(x = goa_strata,
#                               CRSobj = "+proj=utm +units=km +zone=5")

plot(goa_strata)
plot(AK_land, add = TRUE)

## Plot strata, colored by catch rate
dat_strat <- dat %>% group_by(stratum) %>% summarise(mean_cpue = mean(Catch_KG)) # or specific to years diverging?

goa_strata@data <- left_join(goa_strata@data, dat_strat)

# harder to control coloring correctly in base R
#plot(goa_strata, col = goa_strata$mean_cpue, border = NA) 

goa_sf <- st_as_sf(goa_strata)

ggplot() + 
  geom_sf(data = goa_sf, aes(fill = mean_cpue), color = NA) + 
  scale_fill_viridis()
