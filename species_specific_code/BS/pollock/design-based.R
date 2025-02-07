# Design-based index of abundance and age composition using the gapindex 
# pacakge, used for comparison to the model-based indices & comps.
# By: Sophia Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.09.25

# devtools::install_github("afsc-gap-products/gapindex")
library(gapindex)
library(dplyr)
library(here)

# Connect to Oracle (will only work on the NOAA network or VPN)
# Get design-based comps from gapindex ----------------------------------------
if (file.exists("Z:/Projects/ConnectToOracle.R")) {
  source("Z:/Projects/ConnectToOracle.R")
} else {
  # For those without a ConnectToOracle file
  channel <- get_connected()
}

# check to see if connection has been established
odbcGetInfo(channel_products)


# Design-based comps from gapindex --------------------------------------------
gapindex_comps <- function(region) {
  gapindex_data_ebs <- gapindex::get_data(
    year_set = c(1982:2024),
    survey_set = region,
    spp_codes = 21740,
    haul_type = 3,
    abundance_haul = "Y",
    pull_lengths = T,
    channel = channel_products
  )
  
  # Calculate and zero-fill CPUE
  cpue <- calc_cpue(gapdata = gapindex_data_ebs)
  
  # Calculate biomass (w/variance), mean/variance CPUE across strata
  biomass_stratum <- calc_biomass_stratum(gapdata = gapindex_data_ebs, 
                                          cpue = cpue)
  
  
  # Calculate size compositions by stratum
  sizecomp_stratum <- calc_sizecomp_stratum(gapdata = gapindex_data_ebs,
                                            cpue = cpue,
                                            abundance_stratum = biomass_stratum,
                                            spatial_level = "stratum",
                                            fill_NA_method = "BS")
  
  # Calculate age-length key 
  alk <- calc_alk(gapdata = gapindex_data_ebs,
                  unsex = "all")
  
  # Calculate age composition by stratum
  age_comp_stratum <- calc_agecomp_stratum(gapdata = gapindex_data_ebs,
                                           alk = alk,
                                           sizecomp_stratum = sizecomp_stratum)
  
  # Calculate aggregated age composition across regions
  age_comp_region <- calc_agecomp_region(gapdata = gapindex_data_ebs,
                                         agecomp_stratum = age_comp_stratum)
  return(age_comp_region)
}
ebs_comps <- gapindex_comps("EBS")
nbs_comps <- gapindex_comps("NBS")

# Calculate proportions
proportions <- function(df, label) {
  props <- df
  props[props$AGE == 0, ]$AGE <- 1  # create lower plus group
  props[props$AGE >= 16, ]$AGE <- 15  # create upper plus group
  props <- props %>%
    group_by(AGE, YEAR) %>%
    summarize(count = sum(POPULATION_COUNT)) %>% # summarize across sex & plus groups
    group_by(YEAR) %>%
    mutate(total = sum(count)) %>%  # total per year
    group_by(AGE, YEAR) %>%
    summarize(proportion = count / total) %>%  # calculate proportion
    arrange(YEAR, AGE)  # sort to match other output
  props$Region <- label
  
  return(props)
}

props_all <- rbind(proportions(ebs_comps, "EBS"), 
                   proportions(nbs_comps, "NBS"),
                   proportions(rbind(ebs_comps, nbs_comps), "Both"))
  
write.csv(props_all, 
          here("species_specific_code", "BS", "pollock", "results", "Comps", "gapindex_comps_2024.csv"),
          row.names = FALSE)


# # Zack's design-based age comps for numbers of fish & both EBS & NBS ----------
# pollock_bs_age_props_db <- RODBC::sqlQuery(
#   channel = channel_products,
#   query = "
# WITH 
# 
# -- AGGREGATE AGECOMPS ACROSS NBS AND EBS WITH THE PLUS GROUP
# AGE_ALL_W_PLUSGROUP AS (
# SELECT YEAR, 
# CASE
#     WHEN AGE >= 15 THEN 15 -- PLUS GROUP
#     ELSE AGE
# END AS AGE,
# SUM(POPULATION_COUNT) as POPULATION_COUNT
# FROM GAP_PRODUCTS.AGECOMP 
# WHERE SPECIES_CODE = 21740 -- WALLEYE POLLOCK
# AND AREA_ID in (99900, -- EBS STANDARD + NW AREA
# 99902                  -- NBS
# )
# AND AGE >= 0
# 
# GROUP BY YEAR, 
# CASE
#     WHEN AGE >= 15 THEN 15 -- PLUS GROUP
#     ELSE AGE
# END
# ),
# 
# -- AGGREGATE AGECOMPS FOR THE NBS AND EBS (SEPARATLEY) WITH THE PLUS GROUP
# AGE_REGION_W_PLUSGROUP AS (
# SELECT YEAR, AREA_ID,
# CASE
#     WHEN AGE >= 15 THEN 15 -- PLUS GROUP
#     ELSE AGE
# END AS AGE,
# SUM(POPULATION_COUNT) as POPULATION_COUNT
# FROM GAP_PRODUCTS.AGECOMP 
# WHERE SPECIES_CODE = 21740 -- WALLEYE POLLOCK
# AND AREA_ID in (99900, -- EBS STANDARD + NW AREA
# 99902                  -- NBS
# )
# AND AGE >= 0
# 
# GROUP BY YEAR, AREA_ID,
# CASE
#     WHEN AGE >= 15 THEN 15 -- PLUS GROUP
#     ELSE AGE
# END
# )
# 
# SELECT YEAR, 'Both' AS AREA_ID, AGE, 
# ROUND(POPULATION_COUNT * 1.0 / SUM(POPULATION_COUNT) OVER (PARTITION BY YEAR), 7) AS PROPORTION
# FROM AGE_ALL_W_PLUSGROUP
# 
# UNION
# 
# SELECT YEAR, CASE
#     WHEN AREA_ID = 99900 THEN 'EBS'
#     ELSE 'NBS'
# END AS AREA_ID, 
# AGE, 
# ROUND(POPULATION_COUNT * 1.0 / SUM(POPULATION_COUNT) OVER (PARTITION BY AREA_ID, YEAR), 7) AS PROPORTION
# FROM AGE_REGION_W_PLUSGROUP
# 
# ORDER BY AREA_ID, YEAR, AGE                                 
# ")
# 
# ## Query years of age data for each region
# ebs_years <- 
#   unique(x = pollock_bs_age_props_db$YEAR[pollock_bs_age_props_db$AREA_ID == "EBS"])
# nbs_years <- 
#   unique(x = pollock_bs_age_props_db$YEAR[pollock_bs_age_props_db$AREA_ID == "NBS"])
# 
# write.csv(pollock_bs_age_props_db, 
#           here("species_specific_code", "BS", "pollock", "results", "Comps", "DB_comps_num_2024.csv"),
#           row.names = FALSE)
