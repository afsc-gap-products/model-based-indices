##### ORGANIZE AND FORMAT OBSERVATION DATA
# This script organizes and formats the observation data for the spatiotemporal
# weight-at-age model.
# By: Sophia N. Wassermann; modified from Indivero et al. (2023) 
#     https://github.com/jindivero/bs-pollock-weight
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.10.17

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

### CPUE ----------------------------------------------------------------------
## Load data files
wd <- here("species_specific_code", "BS", "pollock", "weight-at-age")
cruise <- read.csv(here(wd, "data", "cruise_data.csv"))

# These files need to be updated every year
alk <- read.csv(here(wd, "data", "age_length_key_full_densdep_corrected_2025.csv"))
specimen <- read.csv(here(wd, "data", "raw_data_pollock_specimen_2025-09-30.csv"))
haul <- read.csv(here(wd, "data", "raw_data_hauls_survey_2025-09-30.csv"))
catch <- read.csv(here(wd, "data", "raw_data_pollock_catch_2025-09-30.csv"))

names(alk) <- tolower(names(alk))
names(specimen) <- tolower(names(specimen))
names(haul) <- tolower(names(haul))
names(catch) <- tolower(names(catch))
names(cruise) <- tolower(names(cruise))

## Merging CPUE with *all* cruise data, retaining zero catches

# Streamline haul data data_sub to just relevant columns
hauls <- merge(catch, haul, by = c("hauljoin", "cruisejoin")) %>%
  select(hauljoin, cruisejoin, cruise.x, haul.x, start_latitude, start_longitude,
         year, weight, number_fish)
colnames(hauls)[3:4] <- c("cruise", "haul")

# Merge hauls and cpue data to get all hauls w/ zero pollock also
cpue2 <- alk %>%
  full_join(hauls, by= c("cruisejoin", "haul", "year", "start_latitude", "start_longitude"))

# Age bin
cpue2$age_bin <- ifelse(cpue2$age>=15, 15, cpue2$age)

# Streamlined data frame
cpue_final2 <- cpue2[,c("year", "age","cruisejoin", "haul", "hauljoin", 
                        "start_latitude", "start_longitude", "age_bin", 
                        "age_cpue_corr", "weight")]

# Make a column adding the cpue for the 15+ ages together
cpue_final3 <- cpue_final2 %>% 
  group_by(year, age_bin, start_latitude, start_longitude) %>%
  mutate(age_cpue_sum = sum(age_cpue_corr)) %>%
  ungroup()

# Filter out duplicate rows w/ multiple 15 age bins
cpue_final <- cpue_final3 %>% 
  group_by(year, age_bin, start_latitude, start_longitude) %>%
  distinct(age_cpue_sum, .keep_all=TRUE) %>%
  ungroup()

# Sanity check that this  is correct and actually does what I want it to
cpue_final3_15 <- subset(cpue_final3, cpue_final3$age_bin == 15)
cpue_final_15b <- cpue_final3_15 %>% 
  group_by(year, age_bin, start_latitude, start_longitude) %>%
  distinct(age_cpue_sum, .keep_all=TRUE) %>%
  ungroup()
#PASSES - SNW: how to tell? 

# Remove hauls w/ biomass measure but no abundance-at-age from dataset
cpue_final$abundance_yn <- ifelse(!is.na(cpue_final$age_cpue_corr), 1, 0)
cpue_final$zero_biomass_yn <- ifelse(cpue_final$weight == 0, 1, 0)
cpue_final$biomass_na_yn <- ifelse(is.na(cpue_final$weight), 1, 0)

cpue_final_edit <- subset(cpue_final, cpue_final$abundance_yn == 1|cpue_final$zero_biomass_yn == 1)
cpue_final_removed <- subset(cpue_final, cpue_final$abundance_yn != 1 & cpue_final$zero_biomass_yn != 1)

# Isolate hauls w/ a biomass measure but no abundance-at-age
missing_abundance <- subset(cpue_final, cpue_final$age_cpue_corr == 0 & cpue_final$weight > 0)
missing_abundance <- subset(missing_abundance, missing_abundance$year != 2021)
missing_abundance <- subset(missing_abundance, missing_abundance$age_cpue_corr != 0)

# save intermediate steps as csv, if wanted
# write.csv(cpue_final_removed, here(wd, "data", "processed", "missing_abundance.csv"),
#           row.names = FALSE)
# write.csv(cpue_final, here(wd, "data", "processed", "cpue_final_no_zero_corrected.csv"),
#           row.names = FALSE)

# save intermediate steps as rds, if wanted
# saveRDS(cpue_final_removed, here(wd, "data", "processed", "missing_abundance.rds"))
# saveRDS(cpue_final, here(wd, "data", "processed", "cpue_final_no_zero_corrected.rds"))

# Remove weird age zeroes
cpue_final_edit2 <- cpue_final_edit %>% 
  filter(age != 0 | is.na(age))

## Add zeroes for hauls missing some age classes
# Sanity check
# Check that there's some zeros
Data=cpue_final_edit2
Data <- as.data.frame(Data)
out <- tapply(Data[,'age_cpue_sum'], INDEX = list(Data[,'year'], Data[,'age_bin']), 
              FUN = function(x){sum(x == 0)})
# FAILS CHECK

# Sanity check
# Check that sample size is equal for all ages in a given year
out2 <- tapply(Data[, 'age_cpue_sum'], INDEX = list(Data[, 'year'], Data[, 'age_bin']), 
               FUN=length)
# FAILS CHECK

# Sanity check
# Check that each hauljoin has 1-15 observations
table(table(Data[, 'hauljoin']))
# PASSES CHECK

# Make column numeric
Data$age_cpue_sum <- as.numeric(Data$age_cpue_sum)

Data_fix <- tapply(Data[,'age_cpue_sum'], INDEX = list("hauljoin"=Data[, 'hauljoin'], 
                                                       "age_bin"= Data[, 'age_bin']), 
                   FUN = sum)
Data_fix <- ifelse(is.na(Data_fix), 0, Data_fix)

# rebuild wide-form data frame
Data_wide <- cbind(Data_fix, Data[match(rownames(Data_fix), Data[, 'hauljoin']), 
                                  c("start_latitude", "start_longitude", "year")])

# build long-form data frame
Data_long <- expand.grid(dimnames(Data_fix))
Data_long <- cbind(Data_long, "age_cpue_sum" = as.vector(Data_fix), 
                   Data[match(Data_long[, 'hauljoin'], Data[, 'hauljoin']),
                        c("start_latitude", "start_longitude", "year")])

# Sanity check
# Check that there's some zeros
out <- tapply(Data_long[,'age_cpue_sum'], INDEX = list(Data_long[, 'year'],
                                                       Data_long[, 'age_bin']), 
              FUN = function(x){sum(x == 0)} )
# PASSES CHECK

# Sanity check
# Check that sample size is equal for all ages in a given year
out2 <- tapply( Data_long[, 'age_cpue_sum'], INDEX=list(Data_long[, 'year'], 
                                                        Data_long[, 'age_bin']), 
                FUN = length )
# PASSES CHECK

glimpse(Data_long)

# Remove age zeros
Data_long <- subset(Data_long, Data_long$age_bin != 0)

# Save final cpue data
write.csv(Data_long, here(wd, "data", "processed", "cpue_final.csv"), row.names = FALSE)
# saveRDS(Data_long, here(wd, "data", "processed", "cpue_final.rds"))


#### SPECIMEN DATA (IE WEIGHT) ------------------------------------------------
# Remove rows where species code == NA
specimen <- specimen %>% filter(species_code == 21740)

# Determine False/True if weight data available
specimen$weight_available <- specimen$weight %>% 
  is.na %>% 
  `!` 
# Determine False/True if age data available
specimen$age_available <- specimen$age%>% 
  is.na %>% 
  `!` 
# Designate categories 1, 2, 3, 4 for level of data available 
specimen$data_available <- case_when(specimen$weight_available == TRUE & specimen$age_available == FALSE ~ 2,
                                     specimen$age_available == TRUE & specimen$weight_available == FALSE  ~ 3,
                                     specimen$age_available == FALSE & specimen$weight_available == FALSE ~ 1,
                                     specimen$age_available == TRUE & specimen$weight_available == TRUE ~ 4)

# Make 15+ bin in new variable called "age_bin"
specimen$age_bin <- ifelse(specimen$age >= 15, 15, specimen$age)

# Merge with haul information

# Merge dataframes and make new dataframe
data2 <- specimen %>% 
  left_join(cruise, by = "cruisejoin")

# Create column w/ EBS and NBS based on stratum
data2$location <- ifelse(data2$stratum == 81 | data2$stratum == 70 |data2$stratum == 71,
                         'northern', 'eastern')

##### Calculate weights from weight-length relationship w=aL^b
# Using parameters from Jim Ianelli's spreadsheet
# 1=male, 2=female
B_m <- 3.038
B_f <- 2.986
B_unknown <- 2.9954
A_m <- 0.000004919
A_f <- 0.000006681
A_unknown <- 0.0000063611

data2$weight_calculated <- case_when(data2$sex == 1 ~ A_m * data2$length ^ B_m,
                                     data2$sex == 2 ~ A_f  *data2$length ^ B_f,
                                     data2$sex == 3 ~ A_unknown * data2$length ^ B_unknown)

# Compare weights, if wanted
plot(data2$weight_calculated ~ data2$weight, 
     ylab = "Calculated weight (w=aL^b)", 
     xlab = "Measured weight")
abline(a = 0,b = 1)
abline(lm(weight_calculated ~ weight, data = data2), col = "blue")

# Save intermediate step as csv, if wanted
# write.csv(data2, file = here(wd, "data", "processed", "specimen_data_edited_complete.csv"), 
#           row.names = FALSE)

# Create streamlined dataframe for individual weight data
specimen_data <- data2[, c('year', 'start_latitude', 'start_longitude', 
                           'region.x','survey_name', 'location','hauljoin', 
                           'cruisejoin', 'sex', 'age_bin', 'age', 'weight', 
                           'weight_calculated','length', 'data_available')]

# Save as csv and rds
write.csv(specimen_data, file = here(wd, "data", "processed", "specimen_data_edited_streamlined.csv"),
          row.names = FALSE)
# saveRDS(specimen_data, file = here(wd, "data", "processed", "specimen_data_edited_streamlined.rds"))


##### Handling extrapolation grid ---------------------------------------------
# # Add ebs and nbs extrapolation grid csv's from Jim Thorson
# ebs_grid <- read.csv("~/Dropbox/Mac/Documents/Documents/Pollock/data/EBSThorsonGrid.csv")
# nbs_grid <- read.csv("~/Dropbox/Mac/Documents/Documents/Pollock/data/NBSThorsonGrid.csv")
# 
# # Add columns for "eastern" and "northern"
# ebs_grid$location <- paste("eastern")
# nbs_grid$location <- paste("northern")
# 
# # Merge grids
# bs_grid <- rbind(ebs_grid, nbs_grid)
# 
# # Convert m^2 to km^2 and label column w/ correct VAST formatting
# bs_grid$Area_km2 <- bs_grid$Shape_Area / 1000000
# 
# # Save extrapolation grid as rds file
# saveRDS(bs_grid, file = "bs_grid.rds")

# Read in extrapolation grid
user_region <- readRDS(here(wd, "data", "bs_grid.rds"))

#### Combine CPUE and specimen (weight) data for model fitting
# Upload CPUE and specimen data, if needed
cpue_final <- read.csv(here(wd, "data", "processed", "cpue_final.csv"))
specimen_data <- read.csv(here(wd, "data", "processed", "specimen_data_edited_streamlined.csv"))

# Combine datasets
cpue_final$start_latitude <- as.numeric(cpue_final$start_latitude)
cpue_final$start_longitude <- as.numeric(cpue_final$start_longitude)

example1 <- specimen_data
example2 <- cpue_final
example1$lat <- example1$start_latitude
example1$lon <- example1$start_longitude
example2$weight <- NULL

# Combine weight_calculated and weight observed
example1$weight_combined <- ifelse(!is.na(example1$weight), example1$weight, example1$weight_calculated)

# bind
example <- bind_rows(example1, example2)

# Remove weird 0 ages
example <- subset(example, example$age_bin > 0)

test <- example1$hauljoin %in% example2$hauljoin

# Limit to before 2019 if wanted (ie exclude 2021 abundance data)
# example <- subset(example, year < 2020)

# Sanity check on data
check <- as.data.frame(example)
check2 <- aggregate(weight_combined ~ year + age_bin, check, FUN = length)
check3 <- check2 %>% spread(key = "age_bin", value = "weight_combined")
ggplot(check2, aes(x = age_bin, y = weight_combined, group = year, color = year)) + 
  geom_line()
check4 <- aggregate(age_cpue_sum ~ year + age_bin, check, FUN = length)
check4 <- check4 %>% spread(key = "age_bin", value = "age_cpue_sum")

# Save combined data
write.csv(example, here(wd, "data", "processed", "data_combined.csv"))
saveRDS(example, here(wd, "data", "processed", "data_combined.rds"))
