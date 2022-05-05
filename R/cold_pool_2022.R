# Method of prepping cold pool covariate for GAP index VAST runs
# author: Caitlin Allen Akselrud
# date: 2020.05.05
# contact: caitlin.allen_akselrud@noaa.gov

# Cold pool covariate -----------------------------------------------------

# devtools::install_github("afsc-gap-products/coldpool")#, lib = .libPaths()[2])
coldpool:::cold_pool_index

cold_pool <- coldpool:::cold_pool_index %>% 
  as_tibble() %>% 
  clean_names #%>%

# cold_pool <- read_csv(here("output", "cold_pool.csv"))

cp <- cold_pool %>% 
  dplyr::select(-year, -last_update) %>% 
  scale(center = TRUE, scale = TRUE) %>% 
  as_tibble() 
missing_yr <- bind_cols(year = 2020, 
                        area_lte2_km2 = 0,
                        area_lte1_km2 = 0,
                        area_lte0_km2 = 0,
                        area_lteminus1_km2 = 0,
                        mean_gear_temperature = 0,
                        mean_bt_lt100m = 0,
                        mean_surface_temperature = 0) 

cp_data <- cold_pool %>% 
  dplyr::select(year) %>%
  bind_cols(cp) %>% 
  dplyr::filter(year %in% unique(Data_Geostat[,'Year'])) %>% 
  full_join(missing_yr) %>% 
  arrange(year)

covariate_data <- data.frame( "Year"=cp_data[,"year"],
                              "Lat"=mean(Data_Geostat[,'Lat']),
                              "Lon"=mean(Data_Geostat[,'Lon']),
                              "area_lte2_km2" = cp_data[,"area_lte2_km2"]) %>%
  rename(Year = year) %>% 
  data.frame()