# Process GOA groundfish bottom trawl survey data for 2021 ModSquad analyses
# Result is cleaned cpue (kg/km^2) by haul, with zeros included
# Requires VPN connection
# Author: Lewis Barnett
# Modified by Madison Hall 05/2021

library(dplyr)

setwd("//nmfs.local/akc-public/Dropbox/mod_squad")

# import flat file of CPUE (kg / km^2) exported from AFSC database
dat <- read.csv("data_raw/cpue_GOA_selected_spp.csv", stringsAsFactors = FALSE) 

# join haul data to get coordinates, and any haul data you may want: depth, bottom and surface temperature
haul <- read.csv("data_raw/haul.csv", stringsAsFactors = FALSE)
haul$DATE <- as.Date(haul$START_TIME, format = "%d-%b-%y")
haul$MONTH <- lubridate::month(haul$DATE)
haul$DAY <- lubridate::day(haul$DATE)
haul$YEAR <- lubridate::year(haul$DATE)

# filter to satisfactory performing standard survey hauls
haul <- haul %>% 
  filter(REGION == "GOA", PERFORMANCE >= 0 & HAUL_TYPE == 3 & ABUNDANCE_HAUL == 'Y') %>% 
  select(HAULJOIN, GEAR_DEPTH, SURFACE_TEMPERATURE, GEAR_TEMPERATURE, START_LATITUDE, START_LONGITUDE, DATE, DAY, MONTH)

dat <- left_join(dat, haul)

# join species names
species_codes <-  read.csv("data_raw/species.csv", stringsAsFactors = FALSE)
species_codes <- select(species_codes, -YEAR_ADDED)
dat <- left_join(dat, species_codes)

# select and rename columns, dropping rows with missing coordinates
dat <- dat %>% select(YEAR, SURVEY, STRATUM, BOTTOM_DEPTH = GEAR_DEPTH, SURFACE_TEMPERATURE, GEAR_TEMPERATURE,
                        CPUE = WGTCPUE, START_LATITUDE, START_LONGITUDE, DATE, DAY, MONTH, SPECIES_NAME, COMMON_NAME, SPECIES_CODE) %>%
                 filter(!is.na(START_LATITUDE), !is.na(START_LONGITUDE))

# write.csv(dat, paste0("data/AK_BTS_data_GOA_modsquad_", Sys.Date(), ".csv"))
# saveRDS(dat, paste0("data/AK_BTS_data_GOA_modsquad_", Sys.Date(), ".rds"))

pop_data <- dat[which(dat$SPECIES_CODE == 30060),]
northern_data<-dat[which(dat$SPECIES_CODE == 30420,)]

saveRDS(pop_data, paste0(getwd(),"/species_specific_code/GOA/Sebastes_alutus/data/Data_Geostat_Sebastes_alutus.rds"))
saveRDS(northern_data, paste0(getwd(),"/species_specific_code/GOA/Sebastes_polyspinis/data/Data_Geostat_Sebastes_polyspinis.rds"))

# Example of filtering to your species of interest (you can use any of species name, common name, or species code)
# readRDS("//nmfs.local/akc-public/Dropbox/mod_squad/data/AK_BTS_data_GOA_modsquad_2021.rds") # or use most recent dated version
# dat_pcod <- filter(dat, SPECIES_NAME == "Gadus macrocephalus")
# after filtering, check that you have 11605 rows for all years through 2019

## CAVEATS
## 1) species included in the GOA CPUE table are those requested for assessments and may not be the most data-rich
## 2) some species have been lumped/split over the years, so make sure you discuss this with ModSquad and assessment author 
## (e.g., rock soles, dusky/dark rockfishes, rougheye and blackspotted rockfishes)
## 3) you may need to conduct further filtering to the years requested by assessment authors (e.g., remove 1980s)