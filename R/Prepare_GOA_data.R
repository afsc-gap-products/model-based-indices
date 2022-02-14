library('TMB')
library("devtools")
library('VAST')
library('FishData')
library(dplyr)

#species_code <- c(21720,30152,30150)[1]#[2:3] #c(P.cod, dusky rockfish, dusky and dark rockfishes unid.)
species_code <- c(30420,30060,21740,10110,10261,10262,10130)[2] #c(northern rockfish, POP, pollock, arrowtooth, northern rock sole, southern rock sole, flathead sole)

#species_code <- 10260 ## rock sole unidentified.
PKG <- c("RODBC")
for (p in PKG) {
  if(!require(p,character.only = TRUE)) {  
    install.packages(p)
    require(p,character.only = TRUE)}
}

# channel<-odbcConnect(dsn = "AFSC",
#                      uid = " ", # change
#                      pwd = " ", #change
#                      believeNRows = FALSE)

# odbcGetInfo(channel)
source("R/get_connected.R")

##need to connect to VPN for this to work
cruise <- RODBC::sqlQuery(channel, "SELECT * FROM RACEBASE.CRUISE")
cruise$YEAR = as.numeric(substring(cruise$START_DATE, 8,9))
cruise$YEAR <- ifelse(cruise$YEAR < 21, cruise$YEAR + 2000, cruise$YEAR + 1900)
cruise = cruise[which(cruise$YEAR >= 1984),]
#cruise = cruise[which(cruise$YEAR >= 1996),]
#cruise = cruise[which(cruise$YEAR >= 1990),]

haul <- RODBC::sqlQuery(channel, "SELECT * FROM RACEBASE.HAUL")
haul = haul[which(haul$REGION =='GOA'),]
haul = haul[which(haul$ABUNDANCE_HAUL=='Y'),]
haul = haul[which(haul$HAUL_TYPE == 3),]
haul = haul[which(haul$PERFORMANCE >= 0),]
haul$YEAR = as.numeric(substring(haul$START_TIME, 1,4))
haul = haul[which(haul$YEAR >= 1984),]
#haul = haul[which(haul$YEAR >= 1996),]
#haul = haul[which(haul$YEAR >= 1990),]

catch <- RODBC::sqlQuery(channel, "SELECT * FROM RACEBASE.CATCH")
catch = catch[which(catch$REGION == "GOA"),]

#dplyr::filter(COMMON_NAME %in% c("dusky and dark rockfishes unid.", "dusky rockfish"))
if(species_code == 21720){
  species_name <- "Gadus_macrocephalus"  #Pcod
}
if(species_code == 30152){
  species_name <- "Sebastes_variabilis" #dusky (and dark_ rockfish)
}
if(species_code == 30150){
  species_name <- "Sebastes_variabilis" #dusky and dark rockfish
}
if(species_code == 30420){
  species_name <- "Sebastes_polyspinis" #northern rockfish
}
if(species_code == 30060){
  species_name <- "Sebastes_alutus" #POP
}
if(species_code == 21740){
  species_name <- "Gadus_chalcogrammus" #pollock
}
if(species_code == 10110){
  species_name <- "Atheresthes_stomias" #arrowtooth
}
if(species_code == 10261){
  species_name <- "Lepidopsetta_polyxystra" #northern rock sole
}
if(species_code == 10262){
  species_name <- "Lepidopsetta_bilineata" #southern rock sole
}
if(species_code == 10130){
  species_name <- "Hippoglossoides_elassodon" #flathead sole
}
if(species_code == 10260){
  species_name <- "Lepidopsetta_sp." #rock sole unid.
}

# Set up folder to store species specific results
folder <- paste0(getwd(),"/species_specific_code/GOA/",species_name)
dir.create(folder)
folder <- paste0(getwd(),"/species_specific_code/GOA/",species_name,"/data")
dir.create(folder)
folder <- paste0(getwd(),"/species_specific_code/GOA/",species_name,"/results")
dir.create(folder)

##for rock sole unidentified (fill in zeros only for years that there is positive catch for this specific species)
if(species_code == 10260){
  catch_by_species <- subset(catch, catch$SPECIES_CODE==species_code)
  catch_by_species <- data.frame(HAULJOIN = catch_by_species$HAULJOIN,
                                 CRUISEJOIN = catch_by_species$CRUISEJOIN,
                                 WEIGHT = catch_by_species$WEIGHT,
                                 NUMBER_FISH = catch_by_species$NUMBER_FISH,
                                 SPECIES_CODE = catch_by_species$SPECIES_CODE,
                                 CATCHJOIN = catch_by_species$CATCHJOIN,
                                 SUBSAMPLE_CODE = catch_by_species$SUBSAMPLE_CODE,
                                 VOUCHER = catch_by_species$VOUCHER)
  data_sub <- left_join(haul, catch_by_species)
  positive_catch_years <- unique(data_sub$YEAR[which(data_sub$WEIGHT != 0)]) ##IDs years with positive catch (for design-based unidenf. rock sole )
  data_sub <- data_sub[which(data_sub$YEAR %in% positive_catch_years),]##subset to only include years with positive catch
}


##for everything that isn't duskies/dark rockfish
if(species_code == 21720 | species_code == 30420 | species_code == 30060 | species_code == 21740 | species_code == 10110| species_code == 10261 | species_code == 10262| species_code == 10130 )
{
  catch_by_species <- subset(catch, catch$SPECIES_CODE==species_code)
  catch_by_species <- data.frame(HAULJOIN = catch_by_species$HAULJOIN,
                               CRUISEJOIN = catch_by_species$CRUISEJOIN,
                               WEIGHT = catch_by_species$WEIGHT, 
                               NUMBER_FISH = catch_by_species$NUMBER_FISH,
                               SPECIES_CODE = catch_by_species$SPECIES_CODE,
                               CATCHJOIN = catch_by_species$CATCHJOIN,
                               SUBSAMPLE_CODE = catch_by_species$SUBSAMPLE_CODE,
                               VOUCHER = catch_by_species$VOUCHER)
  data_sub <- merge(haul, catch_by_species, by="HAULJOIN", all.x=TRUE)
}

if(species_code == 30152 | species_code == 30150){  #this bit deals with duskies & dark
  catch_by_species1 <- subset(catch, catch$SPECIES_CODE==species_code[1])
  catch_by_species2 <- subset(catch, catch$SPECIES_CODE==species_code[2])
  
  catch_by_species1 <- data.frame(HAULJOIN = catch_by_species1$HAULJOIN,
                                  CRUISEJOIN = catch_by_species1$CRUISEJOIN,
                                  WEIGHT = catch_by_species1$WEIGHT, 
                                  NUMBER_FISH = catch_by_species1$NUMBER_FISH,
                                  SPECIES_CODE = catch_by_species1$SPECIES_CODE,
                                  CATCHJOIN = catch_by_species1$CATCHJOIN,
                                  SUBSAMPLE_CODE = catch_by_species1$SUBSAMPLE_CODE,
                                  VOUCHER = catch_by_species1$VOUCHER)
  
  
  catch_by_species2 <- data.frame(HAULJOIN = catch_by_species2$HAULJOIN,
                                  CRUISEJOIN = catch_by_species2$CRUISEJOIN,
                                  WEIGHT = catch_by_species2$WEIGHT, 
                                  NUMBER_FISH = catch_by_species2$NUMBER_FISH,
                                  SPECIES_CODE = catch_by_species2$SPECIES_CODE,
                                  CATCHJOIN = catch_by_species2$CATCHJOIN,
                                  SUBSAMPLE_CODE = catch_by_species2$SUBSAMPLE_CODE,
                                  VOUCHER = catch_by_species2$VOUCHER)
  # ##USING FishStatsUtils
  data_sub1 <- merge(haul[which(haul$YEAR %in% c(1990, 1996:2019)),], catch_by_species1, by="HAULJOIN", all.x=TRUE)
  data_sub1$SPECIES_CODE[which(is.na(data_sub1$SPECIES_CODE))] <- 30152
  data_sub1$WEIGHT[which(is.na(data_sub1$WEIGHT))] <- 0
  data_sub2 <- merge(haul[which(haul$YEAR %in% c(1984, 1987, 1990, 1993)),], catch_by_species2, by="HAULJOIN", all.x=TRUE)
  data_sub2$SPECIES_CODE[which(is.na(data_sub2$SPECIES_CODE))] <- 30150
  data_sub2$WEIGHT[which(is.na(data_sub2$WEIGHT))] <- 0
  
  
  data_sub <-  rbind(data_sub1, data_sub2)
}

################################################################################################################
# ##USING FishStatsUtils
data_sub$WEIGHT[which(is.na(data_sub$WEIGHT))] <- 0
data_sub$NUMBER_FISH[which(is.na(data_sub$NUMBER_FISH))] <- 0
nrow(data_sub[which(data_sub$WEIGHT == 0),])
nrow(data_sub[which(data_sub$NUMBER_FISH == 0),])

data_sub2 <- data_sub %>% 
dplyr::mutate(EFFORT = DISTANCE_FISHED * (NET_WIDTH * 0.001) * 100) %>% 
dplyr::mutate(wCPUE = WEIGHT/EFFORT,
              nCPUE = NUMBER_FISH/EFFORT)

################################################################################################################
GOA_DF <- data_sub2

Data_Geostat <-  transmute(GOA_DF,
                           Catch_KG = wCPUE*100, # sumfish calculates CPUE in kg/ha this converts it to kg/km^2
                           Year = YEAR,
                           Vessel = "missing",
                           AreaSwept_km2 = 1, # area swept is 1 when using CPUE instead of observed weight
                           Lat = START_LATITUDE,
                           Lon = START_LONGITUDE,
                           Pass = 0
)

Data_Geostat$Catch_KG[which(is.na(Data_Geostat$Catch_KG))] <- 0

saveRDS(Data_Geostat,paste0(getwd(),"/species_specific_code/GOA/",species_name,"/data/Data_Geostat_",species_name,".rds"))
