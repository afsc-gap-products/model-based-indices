library(sumfish)
# library(here)
# library(FishData)

species_code <- c(21740, 21720, 10210)[2]

# Get data from RACEBASE

# If using sumfish, set user name and password for Oracle
getSQL()

# Query database for survey data, alternatively, use getSQL() to directly query RACEBASE
EBS <- sumfish::getRacebase(1982:2021, 'EBS_SHELF', speciesCode = species_code) 
EBS_nonstandard_hauls <- EBS$haul_other %>%
  filter(CRUISE %in% c(200101, 199401, 200501, 200601) 
         & PERFORMANCE >= 0
         & HAUL_TYPE == 3
         & !is.na(STRATUM)
  )
EBS_nonstandard_catch <- EBS$catch_other %>%
  filter(HAULJOIN %in% EBS_nonstandard_hauls$HAULJOIN
         )
EBS$catch <- bind_rows(EBS$catch, EBS_nonstandard_catch)
EBS$haul <- bind_rows(EBS$haul, EBS_nonstandard_hauls)

ggplot(EBS_nonstandard_hauls, aes(x=START_LONGITUDE, y=START_LATITUDE)) + 
  geom_point(aes(color = CRUISE))    

NBS <- sumfish::getRacebase(1982:2021, 'NBS_SHELF', speciesCode = species_code) 

# Build NBS 18
NBS18.cruise <- filter(EBS$cruise, CRUISE==201801)
NBS18.haul <- getSQL("select * from racebase.haul where region = 'BS' and cruise = 201801 and haul_type = 13 and performance>=0")
NBS18.haul$STRATUM <- '99'
NBS18.catch <- getSQL(
  paste0("select * from racebase.catch where hauljoin in (select hauljoin from racebase.haul where region = 'BS' and cruise = 201801 and haul_type = 13 and performance>=0) and species_code = ",
         species_code)
         )
NBS18.stratum <- data.frame(SURVEY = 'NBS_18', STRATUM = '99', STRATUM_AREA = 158286)
NBS18.length <- getSQL(
  paste0("select * from racebase.length where hauljoin in (select hauljoin from racebase.haul where region = 'BS' and cruise = 201801 and haul_type = 13 and performance>=0) and species_code = ",
         species_code)
)
NBS18.specimen <- getSQL(
  paste0("select * from racebase.specimen where hauljoin in (select hauljoin from racebase.haul where region = 'BS' and cruise = 201801 and haul_type = 13 and performance>=0) and species_code = ",
         species_code)
  )
NBS18 <- list(cruise = NBS18.cruise,
              haul = NBS18.haul,
              catch = NBS18.catch,
              species = NBS$species,
              stratum = NBS18.stratum,
              length = NBS18.length,
              specimen = NBS18.specimen
) 

strata <- bind_rows(EBS$stratum, NBS$stratum, NBS18$stratum) 

sum18 <- sumHaul(NBS18) %>%
  mutate(STRATUM = as.character(STRATUM),
         REGION = 'NBS')
sumNBS <- sumHaul(NBS) %>%
  mutate(REGION = 'NBS')
sumEBS <- sumHaul(EBS) %>%
  mutate(REGION = 'EBS')


weightAll <- sumAll <- bind_rows(sumEBS, sumNBS,sum18) %>%
  dplyr::filter(!is.na(EFFORT))

# Test for missing values
check_hauls <- c(EBS$haul$HAULJOIN, NBS$haul$HAULJOIN, NBS18$haul$HAULJOIN)

spp_sum <- data.frame()
for (s in species_code) {
  test_spp <- filter(sumAll, SPECIES_CODE == s)
  test_spp <- data.frame(species = s, missing_hauls = length(test_spp[!test_spp$HAULJOIN %in% check_hauls,"HAULJOIN"]) )
  spp_sum <- bind_rows(spp_sum, test_spp)
}



size18 <- sumSize(NBS18) %>%
  mutate(STRATUM = as.character(STRATUM),
         REGION = 'NBS')
sizeNBS <- sumSize(NBS) %>%
  mutate(REGION = 'NBS')
sizeEBS <- sumSize(EBS) %>%
  mutate(REGION = 'EBS')

sizeAll <- sizeComp <- bind_rows(sizeEBS, sizeNBS,size18) %>%
  dplyr::filter(!is.na(EFFORT))

# Test for missing values - hauls with no zeros for species
spp_sum2 <- data.frame()
for (s in species_code) {
  test_spp <- filter(sizeComp, SPECIES_CODE == s)
  test_spp <- data.frame(species = s, missing_hauls = length(test_spp[!test_spp$HAULJOIN %in% check_hauls,"HAULJOIN"]) )
  spp_sum <- bind_rows(spp_sum2, test_spp)
}

# Age-length keys
# For cod, use NBS, then fill by EBS in same year, then NBS all years
# NBS18 has no otoliths, use EBS
# Generate ALKs from sumfish
EBS_4alk <- subsetRacebase(EBS, cruise = 199400:202299)
NBS_4alk <- subsetRacebase(NBS, cruise = 199400:202299)
alk_ebs <- sumALK(EBS_4alk)

alk_nbs_prime <- sumALK(NBS_4alk, global_fill = FALSE)    
alk_nbs_fill <- sumALK(NBS_4alk, global_fill = TRUE)

# Find missing observations of age-at-length for the NBS
missing_test <- alk_nbs_prime %>%
  group_by(SPECIES_CODE, SEX, YEAR, LENGTH) %>%
  summarize(total_ratio = sum(probability) ) %>%
  filter(total_ratio == 0)

# Filter ALK to only lengths that have some probability for age
alk_nbs_prime2 <- anti_join(alk_nbs_prime, missing_test,
                            by = c("YEAR","SPECIES_CODE","SEX","LENGTH")
)

# Append EBS age-at-length for missing values
alk_nbs_yearFill <- alk_ebs %>%
  inner_join(missing_test,
             by = c("YEAR","SPECIES_CODE","SEX","LENGTH")
  ) %>%
  select(-total_ratio) %>%
  bind_rows(alk_nbs_prime2)

# Find missing observations of age-at-length for the filled NBS key   
missing_test2 <- alk_nbs_yearFill %>%
  group_by(SPECIES_CODE, YEAR, SEX, LENGTH) %>%
  summarize(total_ratio = sum(probability) ) %>%
  filter(total_ratio == 0)

# Again, filter ALK to only lengths that have some probability for age 
alk_nbs_yearFill2 <- anti_join(alk_nbs_yearFill, missing_test2,
                               by = c("SPECIES_CODE","SEX","LENGTH")
)

# Append NBS global (pooled years) age-at-length for missing values
### NOTE: we should look at filling EBS ALK with age-at-length from NBS where missing
alk_nbs_globalFill <- alk_nbs_fill %>%
  inner_join(missing_test2,
             by = c("YEAR","SPECIES_CODE","SEX","LENGTH")
  ) %>%
  select(-total_ratio) %>%
  bind_rows(alk_nbs_yearFill)

# Test that all probabilities add to 1
probabilities_test <- alk_nbs_globalFill %>%
  group_by(SPECIES_CODE, SEX, YEAR, LENGTH) %>%
  summarize(total_ratio = sum(probability) )

unique(probabilities_test$total_ratio)

alk <- list(EBS=alk_ebs, NBS=alk_nbs_globalFill)



write_rds(weightAll, "F:/R/VAST2021/data/EBS_NBS_Index.RDS")
write_rds(sizeAll, "F:/R/VAST2021/data/EBS_NBS_SizeComp.RDS")
write_rds(alk, here::here("species_specific_code", "BS", "Pacific_cod_age","data","unstratified_alk_2022_oldcode.RDS"))
write_rds(strata, "F:/R/VAST2021/data/EBS_NBS_strata.RDS")
write_rds(list(EBS=EBS,NBS=NBS,NBS18=NBS18), "F:/R/VAST2021/data/raw_data.RDS")