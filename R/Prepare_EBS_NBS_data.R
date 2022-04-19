##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Prepare EBS/NBS weight and size data
## Author:        Jason Conner (jason.conner@noaa.gov)
## Description:   Template for preparing haul-level CPUE and age composition 
##                CPUE for the species of interest
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import sumfish ----
##   Installation instructions: https://github.com/afsc-gap-products/sumfish
##   Setup Oracle Account
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
library(sumfish)
sumfish::getSQL()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Select species and years ----
##   start_yer and current_year are used when querying RACEBASE via sumfish
##
##   For some species, we want to use all the years from 
##   start_year:current_year when calculating the ALKs but only use 
##   data from after a particular year when fitting VAST (e.g., PCod). Thus, 
##   we declare another variable called min_year to filter years at and after 
##   min_year for that purpose. By default, we set min_year <- start_year.
##
##   plus_group is used for the age composition calculations, where ages at or
##   older than the plus group are grouped as one group. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_name <- "yellowfin_sole"
species_code <- 10210
start_year <- 1982
current_year <- 2021
plus_group <- 20
min_year <- start_year

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create directory to store data products
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
res_dir <- paste0("species_specific_code/BS/", species_name, "/hindcast/data/")
if(!dir.exists(res_dir)) dir.create(res_dir)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull EBS database ----
##   Query database for the standard EBS survey haul and catch data.
##   Append nonstandard well-performing EBS hauls that occurred in 1994, 2001, 
##   2005, and 2006. These tows are not used in the design-based estimates but
##   becasue they follow standard procedures, can be included in a model-
##   based estiamte.
##
##   Notes: haul_type == 3 = standard bottom sample (preprogrammed station)
##          PERFORMANCE 0 is a good performance code, > 0 are satisfactory and
##          < 0 are bad performance codes for a particular haul.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
EBS <- sumfish::getRacebase(year = c(start_year, current_year), 
                            survey = 'EBS_SHELF', 
                            speciesCode = species_code) 

EBS_nonstandard_hauls <- EBS$haul_other %>%
    filter(CRUISE %in% c(200101, 199401, 200501, 200601) 
           & PERFORMANCE >= 0
           & HAUL_TYPE == 3
           & !is.na(STRATUM))

EBS_nonstandard_catch <- EBS$catch_other %>%
    filter(HAULJOIN %in% EBS_nonstandard_hauls$HAULJOIN)

EBS$catch <- bind_rows(EBS$catch, EBS_nonstandard_catch)
EBS$haul <- bind_rows(EBS$haul, EBS_nonstandard_hauls)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull NBS database ----
##   Query database for the standard NBS survey haul and catch data.
##   Append nonstandard well-performing NBS hauls that occurred in 2018 that
##   are in the EBS dataset. The way that this is queried is by selecting
##   records from the Bering Sea (BS), CRUISE is 201801 (Year 2018),
##   and haul_type == 13, meaning it is an index sample tow.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

## Standard NBS Hauls
NBS <- sumfish::getRacebase(year = c(start_year, current_year), 
                            survey = 'NBS_SHELF', 
                            speciesCode = species_code) 

## 2018 NBS cruise
NBS18.cruise <- dplyr::filter(EBS$cruise, CRUISE == 201801)
NBS18.haul <- sumfish::getSQL("select * from racebase.haul where region = 'BS' and cruise = 201801 and haul_type = 13 and performance>=0")
NBS18.haul$STRATUM <- '99'
NBS18.catch <- sumfish::getSQL("select * from racebase.catch where hauljoin in (select hauljoin from racebase.haul where region = 'BS' and cruise = 201801 and haul_type = 13 and performance>=0)")
NBS18.stratum <- data.frame(SURVEY = 'NBS_18', 
                            STRATUM = '99', 
                            STRATUM_AREA = 158286)
NBS18.length <- sumfish::getSQL("select * from racebase.length where hauljoin in (select hauljoin from racebase.haul where region = 'BS' and cruise = 201801 and haul_type = 13 and performance>=0)")
NBS18.specimen <- sumfish::getSQL("select * from racebase.specimen where hauljoin in (select hauljoin from racebase.haul where region = 'BS' and cruise = 201801 and haul_type = 13 and performance>=0)")

## Combine all the NBS18 data types into a list
NBS18 <- list(cruise = NBS18.cruise,
              haul = NBS18.haul,
              catch = NBS18.catch,
              species = NBS$species,
              stratum = NBS18.stratum,
              length = NBS18.length,
              specimen = NBS18.specimen) 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate haul-level CPUEs ----
##   Fill in zeros for hauls that did not encounter the species, 
##   Bind EBS, NBS and NBS18 haul data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
sumEBS <- sumfish::sumHaul(EBS) %>%
    dplyr::mutate(REGION = "EBS")
sumNBS <- sumfish::sumHaul(NBS) %>%
    dplyr::mutate(REGION = "NBS")
sum18 <- sumfish::sumHaul(NBS18) %>%
    dplyr::mutate(STRATUM = as.character(STRATUM),
                  REGION = "NBS")

weightAll <- sumAll <- dplyr::bind_rows(sumEBS, sumNBS,sum18) %>%
    dplyr::filter(SPECIES_CODE %in% species_code, 
                  YEAR >= min_year,
                  !is.na(EFFORT))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Test for missing values
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
check_hauls <- c(EBS$haul$HAULJOIN, NBS$haul$HAULJOIN, NBS18$haul$HAULJOIN)
ifelse(test = length(sumAll$HAULJOIN[!sumAll$HAULJOIN %in% check_hauls]) == 0,
       yes = "No missing hauls",
       no = paste0(length(sumAll$HAULJOIN[!sumAll$HAULJOIN %in% check_hauls]),
                   " missing hauls. Check code."))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Expand Length Freqs to Size Comps ----
##   Filter to choose years >= min_year
##   Bind EBS, NBS and NBS18 size data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
size18 <- sumfish::sumSize(NBS18) %>%
    dplyr::mutate(STRATUM = as.character(STRATUM),
                  REGION = 'NBS')
sizeNBS <- sumfish::sumSize(NBS) %>%
    dplyr::mutate(REGION = 'NBS')
sizeEBS <- sumfish::sumSize(EBS) %>%
    dplyr::mutate(REGION = 'EBS')

sizeAll <- sizeComp <- bind_rows(sizeEBS, sizeNBS,size18) %>%
    dplyr::filter(YEAR >= min_year, 
                  !is.na(EFFORT))

## Test for missing values - hauls with no zeros for species
ifelse(
    test = length(sizeComp$HAULJOIN[!sizeComp$HAULJOIN %in% check_hauls]) == 0,
    yes = "No missing hauls",
    no = paste0(length(sizeComp$HAULJOIN[!sizeComp$HAULJOIN %in% check_hauls]),
                " missing hauls. Check code.")
)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Age-length keys (ALKs) ----
##   These are some notes for PCod, hopefully they extend to other species:
##   
##   An age-length key (in this respect) is the probability of a fish being
##      some integer age given some length, e.g., for a fish 20 cm, there's 
##      a 60% probability it is age 1, 30% probability it is age 2, and a 10%
##      probability it is age 3. sumfish::sumALK calculates this ALK for every
##      combination of year, sex, and length bin category using the lengths 
##      and ages collected on the survey.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
##   EBS ALK ----
##   In the EBS, global_fill == TRUE is the default argument used. This 
##        fills a missing year/sex/length ALK by pooling the data for 
##        that sex/length combination across years and then filling in the ALK. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
alk_ebs <- sumfish::sumALK(EBS)
alk_ebs$REGION <- "EBS"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
##   NBS ALK ----
##   Because the NBS data is more spotty than the EBS, we first specify a 
##      "prime NBS ALK" (alk_nbs_prime) where global_fill == FALSE, meaning 
##      the year-pooling does not occur for missing combos of sex/length ALKs.
##
##   Missing ALKs are first queried from alk_nbs_prime, then are filled in
##      from the EBS (global_fill == TRUE) ALK
##
##   We then specify an NBS ALK where global_fill == TRUE (alk_nbs_fill), 
##   and we use this ALK to fill in any remaining missing NBS ALKs 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
alk_nbs_prime <- sumfish::sumALK(NBS, global_fill = FALSE)    

## Find missing observations of age-at-length for the NBS
missing_test <- alk_nbs_prime %>%
    ## For each combination of sex, year, and length bin, append a 
    ## column called total_ratio that sums probabilities (should either be
    ##  zero for missing or one for existing)
    dplyr::group_by(SEX, YEAR, LENGTH) %>%
    dplyr::summarize(total_ratio = sum(probability) ) %>%
    ## then filter the dataframe to only those records where total_ratio == 0. 
    ## These are the ALKs that we need to fill in 
    dplyr::filter(total_ratio == 0)

## Filter ALK to only lengths that have some probability for age by using the
## dplyr::anti_join() to return all rows from alk_nbs_prime without a match in
## missing_test.
alk_nbs_prime2 <- dplyr::anti_join(x = alk_nbs_prime, 
                                   y = missing_test,
                                   by = c("YEAR", "SEX", "LENGTH"))

## For those missing NBS ALKs, use the ALKs from the EBS, then append to the 
## available NBS ALKs
alk_nbs_yearFill <- alk_ebs %>%
    dplyr::inner_join(missing_test,
                      by = c("YEAR", "SEX", "LENGTH")) %>%
    dplyr::select(-total_ratio) %>%
    dplyr::bind_rows(alk_nbs_prime2)

## For the remaining missing ALKs, we apply the NBS ALK with global_fill == T
## But first, find missing observations of age-at-length for the filled NBS key   
missing_test2 <- alk_nbs_yearFill %>%
    group_by(YEAR, SEX, LENGTH) %>%
    summarize(total_ratio = sum(probability) ) %>%
    filter(total_ratio == 0)

## Again, filter ALK to only lengths that have some probability for age 
alk_nbs_yearFill2 <- anti_join(x = alk_nbs_yearFill, 
                               y = missing_test2,
                               by = c("SEX", "LENGTH"))

## Append NBS global (pooled years) age-at-length for missing values
alk_nbs_fill <- sumfish::sumALK(NBS, global_fill = TRUE)
alk_nbs_globalFill <- alk_nbs_fill %>%
    inner_join(missing_test2,
               by = c("YEAR", "SEX", "LENGTH")) %>%
    select(-total_ratio) %>%
    bind_rows(alk_nbs_yearFill) %>%
    ## We use EBS ALK for 2018 because the 2018 NBS specimen data are not filtered 
    ## out of the ALK generation for the EBS.
    bind_rows(filter(alk_ebs, YEAR == 2018)) %>% 
    mutate(REGION = "NBS")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Test that all probabilities add to 1
##   There are still some sex/length combos that were not observed in either the 
##   EBS or NBS, these are usually the first and/or last length bins. and a lot
##   of the unsexed/length combos
##
##   NOTE from Jason: This is likely an artifact that has been present 
##   throughout the time series, and it will be addressed when we research and
##   adopt a more robust methodology.
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
probabilities_test <- alk_nbs_globalFill %>%
    group_by(SEX, YEAR, LENGTH) %>%
    summarize(total_ratio = sum(probability) )

unique(probabilities_test$total_ratio)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge EBS and NBS ALKs ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
alk <- rbind(alk_ebs, alk_nbs_globalFill)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate Age Comps ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Create a template of hauls and ages to make sure down the line, each 
## unique haul has the full age composition
allCats <- expand.grid(HAULJOIN = unique(weightAll$HAULJOIN), 
                       AGE = unique(alk$AGE[alk$AGE<=plus_group]), 
                       noAge = 0) %>%
    inner_join(weightAll, by = c("HAULJOIN")) 

## Aggregate by Age key
Data <- sizeComp %>%
    ## append the alks to the sizeComp, adding columns `AGE`` and `probabilties``
    dplyr::left_join(alk, by = c("YEAR", "REGION", "LENGTH", 
                                 "SEX", "SPECIES_CODE")) %>%
    ## Calculate age CPUE and truncate ages > plus group to the max age
    dplyr::mutate(ageCPUE = nSizeCPUE * probability,
                  AGE = ifelse(test = AGE > plus_group, 
                               yes = plus_group, no = AGE)) %>% 
    ## sum cpues of length bins for a given ages for each unique haul
    dplyr::group_by(YEAR, REGION, HAULJOIN, STRATUM, 
                    START_LONGITUDE, START_LATITUDE, nCPUE, AGE) %>%
    dplyr::summarize(ageCPUE = sum(ageCPUE),
                     count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(HAULJOIN, AGE, ageCPUE, count) %>%
    ## make sure that all age categories are available for each unique haul and 
    ## if there is no cpue for a given age, fill with zero. 
    dplyr::right_join(allCats, by= c("HAULJOIN", "AGE")) %>%
    dplyr::mutate(ageCPUE = ifelse(test = is.na(ageCPUE), 
                                   yes = noAge, 
                                   no = ageCPUE)) 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format VAST Data ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
data_geostat_index <- 
    with(sumAll,
         data.frame(Region = REGION,
                    Catch_KG = wCPUE, #cpue units: kg per hectare
                    Year = YEAR,
                    Vessel = "missing", 
                    AreaSwept_km2 = 0.01, # converts cpue units to: kg per km^2
                    Lat = START_LATITUDE,
                    Lon = START_LONGITUDE,
                    Pass = 0 ))

data_geostat_agecomps <-  dplyr::transmute(
    Data,
    Catch_KG = ageCPUE,
    Year = YEAR,
    Vessel = "missing",
    Age = AGE,
    AreaSwept_km2 = .01, # Converts CPUE to km^2
    Lat = START_LATITUDE,
    Lon = START_LONGITUDE,
    Pass = 0) %>%
    data.frame()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save output ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

## Index data
write_rds(x = weightAll, 
          file = paste0(res_dir, "EBS_NBS_Index.RDS"))
write_rds(x = data_geostat_index, 
          file = paste0(res_dir, "data_geostat_index.RDS"))

## Age comp data
write_rds(x = sizeAll, 
          file = paste0(res_dir, "EBS_NBS_SizeComp.RDS"))
write_rds(x = data_geostat_agecomps, 
          file = paste0(res_dir, "data_geostat_agecomps.RDS"))

## Strata data
strata <- dplyr::bind_rows(EBS$stratum, NBS$stratum, NBS18$stratum) 
write_rds(x = strata, 
          file = paste0(res_dir, "EBS_NBS_strata.RDS"))

## ALK
write_rds(x = alk, 
          file = paste0(res_dir, "unstratified_alk_2021.RDS"))

## Raw data
write_rds(x = list(EBS = EBS, NBS = NBS, NBS18 = NBS18), 
          file = paste0(res_dir, "raw_data.RDS"))
