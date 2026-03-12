##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## EBS/NBS data pull for ModSquad VAST input via gapindex R package
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Connect to Oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(gapindex)
channel <- gapindex::get_connected(check_access = F)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Species-Specific Constants. Toggle species row
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_info <- data.frame(species_name = c("Pacific_cod"),
                           species_code = c(21720),
                           start_year = 1990,
                           current_year = 2025,
                           plus_group = c(12), 
                           start_year_age = c(1994))[1, ]

## Set constants
start_year <- species_info$start_year
current_year <- species_info$current_year
species_code <- species_info$species_code
species_name <- species_info$species_name
start_year_age <- species_info$start_year_age
plus_group <- species_info$plus_group

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull GOA catch and effort data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_standard_data <- 
  gapindex::get_data(survey_set = "GOA",
                     year_set = start_year:current_year,
                     spp_codes = species_code,
                     pull_lengths = TRUE,
                     taxonomic_source = "GAP_PRODUCTS.TAXONOMIC_CLASSIFICATION",
                     channel = channel)
goa_nonstandard_data <-
  gapindex::get_data(survey_set = "GOA",
                     year_set = start_year:current_year,
                     spp_codes = species_code,
                     abundance_haul = "N", haul_type = 24,
                     pull_lengths = TRUE,
                     taxonomic_source = "GAP_PRODUCTS.TAXONOMIC_CLASSIFICATION",
                     channel = channel)

## Append the haul, catch, length, and specimen data from the non-standard
## hauls with those of the standard hauls
goa_data <- goa_standard_data
goa_data$haul <- rbind(goa_data$haul, goa_nonstandard_data$haul)
goa_data$catch <- rbind(goa_data$catch, goa_nonstandard_data$catch)
goa_data$size <- rbind(goa_data$size, goa_nonstandard_data$size)
goa_data$specimen <- rbind(goa_data$specimen, goa_nonstandard_data$specimen)

goa_cpue <- gapindex::calc_cpue(gapdata = goa_data)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate Total Abundance 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_popn_stratum <- gapindex::calc_biomass_stratum(gapdata = goa_data,
                                                   cpue = goa_cpue)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate numerical CPUE for a given haul/sex/length bin
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sizecomp <- gapindex::calc_sizecomp_stratum(
  gapdata = goa_data,
  cpue = goa_cpue, 
  abundance_stratum = goa_popn_stratum,
  spatial_level = "haul",
  fill_NA_method = "AIGOA")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate globally filled Age-Length Key for EBS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (start_year != start_year_age) { ## i.e., Pcod when age data start 1994
  
  goa_data$survey <- subset(x = goa_data$survey,
                            subset = YEAR >= start_year_age)
  
  goa_data$size <- merge(x = goa_data$size,
                         y = goa_data$cruise[, c("CRUISEJOIN", "CRUISE"), with = F],
                         by = "CRUISEJOIN")
  goa_data$size <- subset(x = goa_data$size, 
                          CRUISE >= start_year_age * 100)
  
  goa_data$specimen <- merge(x = goa_data$specimen,
                             y = goa_data$cruise[, c("CRUISEJOIN", "CRUISE"), with = F],
                             by = "CRUISEJOIN")
  goa_data$specimen <- subset(x = goa_data$specimen, 
                              CRUISE >= start_year_age * 100)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Decompose the numerical CPUE (S_ijklm) for a given haul and sex/length 
##   bin across ages. S_ijklm is the numerical CPUE of the ith station in 
##   stratum j for species k, length bin l, and sex m. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_alk <- gapindex::calc_alk(gapdata = goa_data, 
                              unsex = "unsex", 
                              global = TRUE)
sizecomp <- merge(x = sizecomp,
                  y = goa_alk,
                  by = c("SURVEY", "YEAR", "SPECIES_CODE", 
                         "SEX", "LENGTH_MM"),
                  allow.cartesian = TRUE) |>
  transform(AGE_CPUE_NOKM2 = S_ijklm_NOKM2 * AGE_FRAC)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Aggregate numerical CPUE across lengths for a given age/sex/haul. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
age_cpue <- rbind(
  ## Aggregate ages younger than the `plus_group`
  stats::aggregate(AGE_CPUE_NOKM2 ~ AGE + HAULJOIN + SPECIES_CODE,
                   data = sizecomp,
                   FUN = sum,
                   drop = F,
                   subset = AGE < plus_group),
  ## Aggregate ages at or older than the `plus_group` as one age
  cbind(AGE = plus_group,
        stats::aggregate(AGE_CPUE_NOKM2 ~ HAULJOIN + SPECIES_CODE,
                         data = sizecomp,
                         FUN = sum,
                         drop = F,
                         subset = AGE >= plus_group))
)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Zero-fill CPUEs for missing ages. First we create a grid of all possible
##   HAULJOINs and ages. But, there are hauls with positive numerical catches
##   but no associated size data. These hauls will be removed.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Hauls with zero catch
unique_hauls_cpue_zeros <- 
  sort(x = unique(x = goa_cpue$HAULJOIN[goa_cpue$CPUE_NOKM2 == 0]))
## Hauls with postive count data
unique_hauls_cpue_pos <- 
  sort(x = unique(x = goa_cpue$HAULJOIN[goa_cpue$CPUE_NOKM2 > 0]))
## Hauls with size data
unique_hauls_sizecomp <- 
  sort(x = unique(x = sizecomp$HAULJOIN))

every_combo_of_ages <- 
  expand.grid(HAULJOIN = c(unique_hauls_cpue_zeros,
                           unique_hauls_cpue_pos[unique_hauls_cpue_pos %in% 
                                                   unique_hauls_sizecomp]),
              SPECIES_CODE = species_code,
              AGE = min(age_cpue$AGE, na.rm = TRUE):plus_group)

## Merge the age_cpue onto the `every_combo_of_ages` df using "HAULJOIN", 
## "SPECIES_CODE", and "AGE" as a composite key. all.x = TRUE will reveal
## missing ages denoted by an NA AGE_CPUE_NOKM2 value
every_combo_of_ages <- merge(x = every_combo_of_ages,
                             y = age_cpue,
                             by = c("HAULJOIN", "SPECIES_CODE", "AGE"),
                             all.x = TRUE)

## Turn NA AGE_CPUE_NOKM2 values to zero
every_combo_of_ages$AGE_CPUE_NOKM2[
  is.na(x = every_combo_of_ages$AGE_CPUE_NOKM2)
] <- 0

## Append the latitude and longitude information from the `ebs_nbs_cpue` df
## using "HAULJOIN" as a key. 
age_cpue <- merge(x = every_combo_of_ages,
                  y = goa_cpue[, c("HAULJOIN", "SURVEY", "YEAR",
                                   "LATITUDE_DD_START", "LONGITUDE_DD_START")],
                  by = "HAULJOIN")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format Data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
data_geostat_agecomps <- with(age_cpue, 
                              data.frame(region = SURVEY,
                                         hauljoin = HAULJOIN,
                                         year = as.integer(YEAR),
                                         lon = LONGITUDE_DD_START,
                                         lat = LATITUDE_DD_START,
                                         age = AGE,
                                         cpue_n_km2 = AGE_CPUE_NOKM2)) 

if (start_year != start_year_age) {
  data_geostat_agecomps <- subset(x = data_geostat_agecomps,
                                  subset = year >= start_year_age)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save output
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir_out <- paste0("species_specific_code/BS/", species_name, "/", phase, "/data/")
# if (!dir.exists(paths = dir_out)) dir.create(path = dir_out, recursive = T)
# for (ifile in c("data_geostat_index", 
#                 "data_geostat_agecomps", 
#                 "sizecomp", "alk", "ebs_data", "nbs_data")) 
#   saveRDS(object = get(x = ifile), 
#           file = paste0(dir_out, ifile, ".RDS"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot Age-Length Key
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pdf(file = paste0(dir_out, "/goa_alk.pdf"), width = 8, height = 11, onefile = TRUE)
par(mfrow = c(3, 2), mar = c(4, 3, 3, 2), oma = c(1, 1, 0, 0))
for (iyear in sort(unique(goa_alk$YEAR)) ){
  for (isex in 1:2) {
    temp_alk <- 
      tidyr::spread(
        data = subset(goa_alk,
                      subset = SEX == isex & YEAR == iyear,
                      select = c("LENGTH_MM", "AGE", "AGE_FRAC")),
        key = AGE, value = AGE_FRAC)
    
    temp_alk[is.na(x = temp_alk)] <- 0
    
    image(y = seq(min(as.integer(names(x = temp_alk)[-1])), 
                  max(as.integer(names(x = temp_alk)[-1])) ), 
          x = seq(min(temp_alk$LENGTH_MM), max(temp_alk$LENGTH_MM), by = 10), 
          z = as.matrix(temp_alk[, -1]),
          col = colorRampPalette(colors = c("white", 
                                            hcl.colors(12, "YlOrRd", 
                                                       rev = TRUE),
                                            "black"))(100), 
          axes = F, ann = F)
    box(); mtext(side = 1, text = "Length (mm)", line = 3); 
    mtext(side = 2, text = "Age (yr)", line = 2.5); 
    mtext(side = 3, text = c("Males", "Females", "Unsexed")[isex], font = 2)
    axis(side = 1)
    axis(side = 2, las = 1, at = seq(min(goa_alk$AGE), max(goa_alk$AGE)))
    abline(h = seq(min(goa_alk$AGE), max(goa_alk$AGE)), lwd = 0.5,
           col = "grey", lty = "dashed")
    
    
    legend("topleft", legend = paste(iyear), bty = "n", cex = 1.5)
    
    if (isex == 3)
      legend("right", legend = seq(1, 0, -0.1), 
             fill = rev(c("white",  hcl.colors(9, "YlOrRd", 
                                               rev = TRUE),
                          "black")), bty = "n")
  }
}
# dev.off()
