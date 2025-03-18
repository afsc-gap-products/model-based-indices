###############################################################################
## Project:         ESP input data formatting 
##
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
##
## Description:     Import COG and log-Effective Area csv output files from 
##                  VAST and reformat to fit the format of the ESP submission
##                  tool
##
## Notes:           ESP contact: Kalei Shotwell (kalei.shotwell@noaa.gov)
###############################################################################

#rm(list = ls())

## Constants
species <- c("Pacific cod", "walleye pollock")
years_to_include <- c(seq(from = 1984, to = 1999, by = 3),
                      seq(from = 2001, to = 2021, by = 2))


## Loop over species, modify output names, save csv

for (ispecies in species) {
  #Create result directory if not created yet
  if(!dir.exists(paste0("ESP_output/", ispecies, "/"))) 
    dir.create(paste0("ESP_output/", ispecies, "/"), recursive = TRUE)
  
  #Import COG csv
  temp_cog <- read.csv(paste0(getwd(), "/output/", 
                              ispecies, " rotated axis/", "COG.csv"))
  
  #Reformat temp_cog to fit ESP format
  temp_cog_esp <- rbind(
    data.frame(YEAR = temp_cog$Year,
               INDICATOR_NAME = paste0("Summer_", 
                                       c("Pacific cod" = "Pacific Cod",
                                         "walleye pollock" = "Pollock")[ispecies],
                                       "_Center_Gravity_Axis", 
                                       temp_cog$m, "_WCGOA_Model"),
               DATA_VALUE = temp_cog$COG_hat
    ),
    data.frame(YEAR = temp_cog$Year,
               INDICATOR_NAME = paste0("Summer_", 
                                       c("Pacific cod" = "Pacific Cod",
                                         "walleye pollock" = "Pollock")[ispecies],
                                       "_Center_Gravity_Axis", 
                                       temp_cog$m, "_STD_ERROR_WCGOA_Model"),
               DATA_VALUE = temp_cog$SE
    )
  )
  
  #Sort by year and write csv
  temp_cog_esp <- subset(x = temp_cog_esp[order(temp_cog_esp$YEAR), ],
                         subset = YEAR %in% years_to_include)
  
  write.csv(x = temp_cog_esp, 
            file = paste0("ESP_output/", ispecies, "/Summer_", 
                          c("Pacific cod" = "Pacific Cod",
                            "walleye pollock" = "Pollock")[ispecies],
                          "_Center_Gravity_WCGOA_Model_2021.csv"),
            row.names = FALSE)
  
  #Import log-effective area occupied csv
  temp_area_occ <- read.csv(paste0(getwd(), "/output/", ispecies, 
                                   " rotated axis/", "ln_effective_area.csv"))
  
  #Reformat temp_area_occ to fit ESP format
  temp_area_occ_esp <- rbind(
    data.frame(YEAR = temp_area_occ$year,
               INDICATOR_NAME = paste0("Summer_", 
                                       c("Pacific cod" = "Pacific Cod",
                                         "walleye pollock" = "Pollock")[ispecies],
                                       "_Area_Occupied_WCGOA_Model"),
               DATA_VALUE = temp_area_occ[, "Estimate"]
    ),
    
    data.frame(YEAR = temp_area_occ$year,
               INDICATOR_NAME = paste0("Summer_", 
                                       c("Pacific cod" = "Pacific Cod",
                                         "walleye pollock" = "Pollock")[ispecies],
                                       "_Area_Occupied_Std_Error_WCGOA_Model"),
               DATA_VALUE = temp_area_occ[, "Std..Error"]
    )
  )
  
  #Sort by year and write csv
  temp_area_occ_esp <- 
    subset(x = temp_area_occ_esp[order(temp_area_occ_esp$YEAR), ],
           subset = YEAR %in% years_to_include)
  
  write.csv(x = temp_area_occ_esp, 
            file = paste0("ESP_output/", ispecies, "/Summer_", 
                          c("Pacific cod" = "Pacific Cod",
                            "walleye pollock" = "Pollock")[ispecies],
                          "_Area_Occupied_WCGOA_Model_2021.csv"),
            row.names = FALSE)
  
}
