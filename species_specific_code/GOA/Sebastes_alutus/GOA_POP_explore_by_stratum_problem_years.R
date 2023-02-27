library(dplyr)
library(rgdal)
library(rgeos)
library(sp)
library(sf)
library(terra)
library(ggplot2)
library(viridis)

dat <- readRDS(paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/data/Data_Geostat_Sebastes_alutus_stratum.rds"))

# design-based ----
#AK_land <- rgdal::readOGR("C:/Users/lewis.barnett/Work/AFSC/Rprojects/Optimal_Allocation_GoA/data/shapefiles/AKland.shp")
# AK_land <- sp::spTransform(x = AK_land,
#                            CRSobj = "+proj=utm +units=km +zone=5")

goa_strata <- rgdal::readOGR("C:/Users/lewis.barnett/Work/AFSC/Rprojects/Optimal_Allocation_GoA/data/shapefiles/goa_strata.shp")
names(goa_strata@data) <- tolower(names(goa_strata@data))
# goa_strata <- sp::spTransform(x = goa_strata,
#                               CRSobj = "+proj=utm +units=km +zone=5")

#plot(goa_strata)
#plot(AK_land, add = TRUE)

## Plot strata, colored by catch rate
dat_strat <- dat %>% filter(Year %in% c(2013,2015,2017,2019,2021)) %>% group_by(stratum) %>% summarise(mean_cpue = mean(Catch_KG)) # or specific to years diverging?

goa_strata@data <- left_join(goa_strata@data, dat_strat)

# harder to control coloring correctly in base R
#plot(goa_strata, col = goa_strata$mean_cpue, border = NA) 

goa_sf <- st_as_sf(goa_strata)

ggplot() + 
  geom_sf(data = goa_sf, aes(fill = mean_cpue), color = NA) + 
  scale_fill_viridis() +
  ggtitle("design-based: problem years")
ggsave(paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/results/design_based_mean_density_problem_years.png"),
       width = 6,
       height = 2,
       units = "in")

# model-based ----
m <- readRDS(file = paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/results/Sebastes_alutusVASTfit.RDS"))

# extract predictions and aggregate within polygons of goa_strata
d <- rowMeans(m$Report$D_gct[,1,match(c(2013,2015,2017,2019,2021), min(dat$Year):max(dat$Year))]) # or specific to years diverging?
d_df <- cbind(as.data.frame(m$extrapolation_list$Data_Extrap)[,1:2], d)

dgeo <- SpatialPoints(d_df[c("Lon", "Lat")], CRS("+proj=longlat +ellps=WGS84 +datum=NAD83 +no_defs")) 
dgeo <- spTransform(dgeo, CRS("+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))

#goageo <- spTransform(goa_strata, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

p <- vect(dgeo)
v <- vect(goa_strata)

values(p) <- data.frame(val=d_df[,3]) #attach predicted values to points

r <- relate(v, p, "intersects")
a <- apply(r, 1, function(i) mean(p$val[i]))
v$mean_cpue <- a

#plot(v, 'pmean')

v_sf <- st_as_sf(v)

ggplot() + 
  geom_sf(data = v_sf, aes(fill = mean_cpue), color = NA) + 
  scale_fill_viridis() +
  ggtitle("model-based: problem years")
ggsave(paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/results/model_based_mean_density_problem_years.png"),
       width = 6,
       height = 2,
       units = "in")
