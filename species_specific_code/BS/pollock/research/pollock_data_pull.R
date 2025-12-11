#' Script for pulling CPUE by numbers for Pribilof pollock indices project. 
#' Based on the Prepare_Bering_data.R script.

# devtools::install_github(repo = "afsc-gap-products/gapindex@v3.0.2", 
#                          dependencies = TRUE)
library(gapindex)
library(here)
library(dplyr)

# Connect to Oracle
if (file.exists("Z:/Projects/ConnectToOracle.R")) {
  source("Z:/Projects/ConnectToOracle.R")
} else {
  # For those without a ConnectToOracle file
  channel <- gapindex::get_connected(check_access = F)
}

## checks to see if connection has been established
odbcGetInfo(channel)

# Pull catch and effort data --------------------------------------------------
species_code <- c(21740, 21741)

# First, pull data from the standard EBS stations
ebs_standard_data <- get_data(year_set = 1982:as.integer(format(Sys.Date(), "%Y")),
                              survey_set = "EBS",
                              spp_codes = species_code,
                              pull_lengths = FALSE, 
                              haul_type = 3, 
                              abundance_haul = "Y",
                              channel = channel,
                              remove_na_strata = TRUE)

#' Next, pull data from hauls that are not included in the design-based index
#' production (abundance_haul == "N") but are included in VAST. By default, the 
#' gapindex::get_data() function will filter out hauls with negative performance 
#' codes (i.e., poor-performing hauls).
ebs_other_data <- get_data(year_set = c(1994, 2001, 2005, 2006),
                           survey_set = "EBS",
                           spp_codes = species_code,
                           pull_lengths = FALSE, 
                           haul_type = 3, 
                           abundance_haul = "N",
                           channel = channel, 
                           remove_na_strata = TRUE)

# Combine the EBS standard and EBS other data into one list. 
ebs_data <- list(
  survey = ebs_standard_data$survey,
  survey_design = ebs_standard_data$survey_design,
  #' Some cruises are shared between the standard and other EBS cruises, so the 
  #' unique() wrapper is there to remove duplicate cruise records. 
  cruise = unique(rbind(ebs_standard_data$cruise,
                        ebs_other_data$cruise)),
  haul = rbind(ebs_standard_data$haul,
               ebs_other_data$haul),
  catch = rbind(ebs_standard_data$catch,
                ebs_other_data$catch),
  species = ebs_standard_data$species,
  strata = ebs_standard_data$strata)

# Calculate CPUE and export ---------------------------------------------------
ebs_cpue <- calc_cpue(gapdata = ebs_data) %>%
  select("YEAR", "LATITUDE_DD_START",
         "LONGITUDE_DD_START", "CPUE_NOKM2") %>%
  transmute(cpue = CPUE_NOKM2, 
            year = as.integer(YEAR),
            lat = LATITUDE_DD_START,
            lon = LONGITUDE_DD_START)

write.csv(ebs_cpue, 
          here("species_specific_code", "BS", "pollock", "research", "pollock_num.csv"), 
          row.names = FALSE)
