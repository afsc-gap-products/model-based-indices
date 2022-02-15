CP <- read.csv(here::here("data","cpa_out_ste_simplified.csv"))
covariate_data <- CP[ which(as.integer(CP[,'YEAR']) %in% unique(Data_Geostat[,'Year'])), ]
covariate_data <- data.frame( "Year"=covariate_data[,"YEAR"], 
                              "Lat"=mean(Data_Geostat[,'Lat']), 
                              "Lon"=mean(Data_Geostat[,'Lon']), 
                              apply(covariate_data[c(-1,-3)],
                                    MARGIN=2,
                                    FUN=function(vec){(vec-mean(vec))/sd(vec)}) ) %>% 
  bind_rows(c(Year=2020L, 
              "Lat"=mean(Data_Geostat[,'Lat']), 
              "Lon"=mean(Data_Geostat[,'Lon']),
              AREA_LTE2_KM2=0)) %>%
  data.frame()

saveRDS(covariate_data, file = "Data_ColdPool.RDS")
