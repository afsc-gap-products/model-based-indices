library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(viridis)

# load data passed to fit_model, with stratum field from RACEBASE
dat <- readRDS(paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/data/Data_Geostat_Sebastes_alutus_stratum.rds"))

# specify years of interest and remaining years to compare with
focal_years <- c(2013,2015,2017,2019,2021)
other_years <- unique(dat$Year)[!unique(dat$Year) %in% focal_years]
years <- other_years

# design-based ----

# shapefile source on google drive: https://drive.google.com/drive/folders/1pxq5Z_TDHMiZ_0WU0Wx6o3lJGJmuqe_L?usp=share_link
goa_strata <- st_read("C:/Users/lewis.barnett/Work/AFSC/Rprojects/Optimal_Allocation_GoA/data/shapefiles/goa_strata.shp")
names(goa_strata) <- tolower(names(goa_strata))

## Plot strata, colored by catch rate
dat_strat <- dat %>% filter(!Year %in% years) %>% group_by(stratum) %>% summarise(mean_cpue = mean(Catch_KG)) 
goa_strata <- left_join(goa_strata, dat_strat, by = "stratum")

ggplot() + 
  geom_sf(data = goa_strata, aes(fill = mean_cpue), color = NA) + 
  scale_fill_viridis() +
  ggtitle("design-based")
ggsave(paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/results/design_based_mean_density_",paste0(min(years),"_",max(years)),".png"),
       width = 6,
       height = 2,
       units = "in")

# model-based ----

# load model fit
m <- readRDS(file = paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/results/Sebastes_alutusVASTfit.RDS"))

# extract predictions and aggregate within strata polygons
d <- rowMeans(m$Report$D_gct[,1,match(years, min(dat$Year):max(dat$Year))]) 
d_df <- cbind(as.data.frame(m$extrapolation_list$Data_Extrap)[,1:2], d)

dgeo <- st_as_sf(d_df, coords = c("Lon", "Lat"), crs = "+proj=longlat +ellps=WGS84 +datum=NAD83 +no_defs")
dgeo <- st_transform(dgeo, st_crs(goa_strata))

p <- vect(dgeo)
v <- vect(goa_strata)

values(p) <- data.frame(val=d_df[,3]) #attach predicted values to points

r <- relate(v, p, "intersects")
a <- apply(r, 1, function(i) mean(p$val[i]))
v$mean_cpue <- a

v_sf <- st_as_sf(v)

ggplot() + 
  geom_sf(data = v_sf, aes(fill = mean_cpue), color = NA) + 
  scale_fill_viridis() +
  ggtitle("model-based")
ggsave(paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/results/model_based_mean_density_",paste0(min(years),"_",max(years)),".png"),
       width = 6,
       height = 2,
       units = "in")

# compute and plot contrast between db and mb (db-mb)
joined_sf <- st_join(goa_strata, v_sf) %>% 
  mutate(diff_cpue = mean_cpue.x - mean_cpue.y, lr_cpue = log(mean_cpue.x / mean_cpue.y))
  
ggplot() + geom_sf(data = joined_sf, aes(fill = diff_cpue), color = NA) + 
  scale_fill_gradient2() +
  ggtitle("difference in predictions (db-mb)") +
  labs(fill = "cpue difference")
ggsave(paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/results/difference_mean_density_",paste0(min(years),"_",max(years)),".png"),
       width = 6,
       height = 2,
       units = "in")

ggplot() + geom_sf(data = joined_sf, aes(fill = lr_cpue), color = NA) + 
  scale_fill_gradient2() +
  ggtitle("log-ratio of predictions log(db/mb)") +
  labs(fill = "log-ratio cpue")
ggsave(paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/results/logratio_mean_density_",paste0(min(years),"_",max(years)),".png"),
       width = 6,
       height = 2,
       units = "in")