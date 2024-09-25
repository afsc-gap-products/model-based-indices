# Design-based index of abundance and age composition using the gapindex 
# pacakge, used for comparison to the model-based indices & comps.
# By: Sophia Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.09.25

library(gapindex)

# Connect to Oracle (will only work on the NOAA network or VPN)
sql_channel <- gapindex::get_connected()

# Pull data for EBS pollock from Oracle
data <- gapindex::get_data(
  year_set = 1982:lubridate::year(today()),
  survey_set = "EBS",
  spp_codes = 21740,  # only pollock
  pull_lengths = TRUE,
  haul_type = 3,
  abundance_haul = "Y",
  sql_channel = sql_channel
)

# Calculate and zero-fill CPUE
cpue <- calc_cpue(racebase_tables = data)

# Calculate biomass (w/variance), mean/variance CPUE across strata
biomass_stratum <- calc_biomass_stratum(racebase_tables = data, 
                                        cpue = cpue)
# Calculate size compositions by stratum
sizecomp_stratum <- calc_sizecomp_stratum(racebase_tables = data,
                                          racebase_cpue = cpue,
                                          racebase_stratum_popn = biomass_stratum,
                                          spatial_level = "stratum",
                                          fill_NA_method = "BS")

# Calculate age-length key 
alk <- subset(x = calc_alk(racebase_tables = data,
                           unsex = c("all", "unsex")[1],
                           global = FALSE),
              subset = AGE_FRAC > 0)

