library(dplyr)
library(VAST)
library(tictoc)

# Set species, model -------------------------------------------------------

which_model <- c("hindcast", "production")[2]
compare <- FALSE # If compare = TRUE, using prior year's alk
species <- 21720
species_name <- "pacific_cod"

workDir <- paste0(getwd(),"/species_specific_code/BS/", 
                  species_name, "/", which_model, "/")
if(!dir.exists(workDir))
  dir.create(path = workDir, recursive = TRUE)
if(!dir.exists(paste0(workDir, "results_age/")))
  dir.create(path = paste0(workDir, "results_age/"), recursive = TRUE)

# Record sessionInfo -------------------------------------------------------
sink(file = paste0(workDir, "results_age/session_info.txt"), 
     type = "output")
sessionInfo()
sink()

# Make sure package versions are correct for current year ------------------
current_year <- 2024
prev_year <- current_year-1
VAST_cpp_version <- "VAST_v14_0_1"
pck_version <- c("VAST" = "3.10.0",
                 "FishStatsUtils" = "2.12.0",
                 "Matrix" = "1.5-3",
                 "TMB" = "1.9.2",
                 "DHARMa" = "0.4.6")

for (pck in 1:length(pck_version)) {
  temp_version <- packageVersion(pkg = names(pck_version)[pck])
  
  if(temp_version == pck_version[pck])
    message(paste0("The version of the '", names(pck_version)[pck], 
                   "' package (", temp_version, ") is consistent",
                   " with the ", current_year, " TOR."))
  
  if(!temp_version == pck_version[pck])
    message(paste0("WARNING: ", 
                   "The version of the '", names(pck_version)[pck], 
                   "' package (", temp_version, ") is NOT consistent",
                   " with the ", current_year, " TOR. Please update the '", 
                   names(pck_version)[pck], "' package to ", 
                   pck_version[pck]))
  
  rm(pck, temp_version)
}

# Import data -------------------------------------------------------------
Data_Geostat <- readRDS(file = paste0(workDir, 
                                      "data/data_geostat_agecomps.RDS"))

# Settings ----------------------------------------------------------------
Region <- c("eastern_bering_sea","northern_bering_sea")
Method <- "Mesh"
knot_method <- "grid"
grid_size_km <- 25
n_x <- 50 # Specify number of "knots"
FieldConfig <- c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig <- c("Beta1"=0, "Beta2"=0, "Epsilon1"=4, "Epsilon2"=4)
OverdispersionConfig <- c("Eta1"=0, "Eta2"=0)
ObsModel <- c(2,4)
Options <- c("Calculate_Range" = FALSE, "Calculate_effective_area" = FALSE, "treat_nonencounter_as_zero" = TRUE )
Aniso <- FALSE
Npool <- 100 #20
fine_scale <- TRUE
BiasCorr <- TRUE
max_cells <- 2000

strata.limits <- data.frame('STRATA'="All_areas")

# Make settings 
settings <- make_settings( 
                      n_x = n_x,
                      Region = Region,
                      purpose = "index2",
                      fine_scale = fine_scale,
                      strata.limits = strata.limits,
                      ObsModel = ObsModel,
                      FieldConfig = FieldConfig,
                      RhoConfig = RhoConfig,
                      OverdispersionConfig = OverdispersionConfig,
                      Options = Options,
                      use_anisotropy = Aniso,
                      Version = VAST_cpp_version,
                      max_cells = max_cells,
                      knot_method = knot_method,
                      bias.correct = BiasCorr
)

strata_names = c("Both","EBS","NBS")

# Fit model ------------------------------------------------------------

tic("running model")
fit = fit_model( "settings"=settings, 
                 "Lat_i"=Data_Geostat[,'Lat'], 
                 "Lon_i"=Data_Geostat[,'Lon'], 
                 "t_i"=Data_Geostat[,'Year'],  # "t_i"=rep(2019,nrow(Data_Geostat)),
                 "c_i"=Data_Geostat[,'Age'], 
                 "b_i"=Data_Geostat[,'Catch_N'], 
                 "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                 "v_i"=Data_Geostat[,'Vessel'],
                 Npool = Npool, 
                 test_fit=FALSE, # set to FALSE if want to avoid interruption due to marginal fit issues
                 create_strata_per_region=TRUE,
                 "working_dir" = paste0(workDir,"/results_age/"),
                 "CompileDir" = paste0(workDir,"/results_age/")
)
toc()  

# Save results
saveRDS(fit, file = paste0(workDir,"/results_age/",species_name,"_VASTfit.RDS"))

## Parameter estimates
saveRDS(object = fit$ParHat, 
        file = paste0(workDir, "results_age/starting_parameters.RDS"))

## General output plots, DHARMa residuals
results <- FishStatsUtils::plot_results( 
  fit = fit, 
  working_dir = paste0(workDir, "results_age/"),
  plot_set = NULL,
  strata_names = strata_names, 
  check_residuals = TRUE)

saveRDS(object = results, 
        file = paste0(workDir, "results_age/VASTresults.RDS"))

## Mapping information
map_list = FishStatsUtils::make_map_info( 
  "Region" = settings$Region, 
  "spatial_list" = fit$spatial_list, 
  "Extrapolation_List" = fit$extrapolation_list)

## Predicted density maps for each age group across years
n_ages <- length(unique(Data_Geostat$Age))
Year_Set <- fit$year_labels

if(!dir.exists(paste0(workDir, "results_age/predicted_density/"))){
  dir.create(paste0(workDir, "results_age/predicted_density/"))
}

FishStatsUtils::plot_maps( 
  fit = fit, 
  plot_set = 3,
  n_cells = 200,
  category_names = c(paste0("Age ", 1:(n_ages-1)), 
                     paste0("Age ", n_ages, "+")),
  Obj = fit$tmb_list$Obj, 
  year_labels = Year_Set,
  PlotDF = map_list[["PlotDF"]],
  working_dir = paste0(workDir, "results_age/predicted_density/")) 

## Predicted Spatial Fields
if(!dir.exists(paste0(workDir, "results_age/spatial_effects/"))){
  dir.create(paste0(workDir, "results_age/spatial_effects/"))
}
FishStatsUtils::plot_maps( 
  fit = fit, 
  plot_set = 16:17,
  n_cells = 200,
  category_names = c(paste0("Age ", 1:(n_ages-1)), 
                     paste0("Age ", n_ages, "+")),
  year_labels = Year_Set,
  Obj = fit$tmb_list$Obj, 
  PlotDF = map_list[["PlotDF"]],
  working_dir = paste0(workDir, "results_age/spatial_effects/")) 

## Predicted Spatiotemporal Fields
if(!dir.exists(paste0(workDir, "results_age/spatiotemporal_effects/"))){
  dir.create(paste0(workDir, "results_age/spatiotemporal_effects/"))
}
FishStatsUtils::plot_maps( 
  fit = fit, 
  plot_set = 6:7,
  n_cells = 200,
  category_names = c(paste0("Age ", 1:(n_ages-1)), 
                     paste0("Age ", n_ages, "+")),
  year_labels = Year_Set,
  Obj = fit$tmb_list$Obj, 
  PlotDF = map_list[["PlotDF"]],
  working_dir = paste0(workDir, "results_age/spatiotemporal_effects/")) 

## Predicted Proportions
if(!dir.exists(paste0(workDir, "results_age/proportions/"))){
  dir.create(paste0(workDir, "results_age/proportions/"))
}

proportions <- FishStatsUtils::calculate_proportion( 
  TmbData = fit$data_list, 
  Index = results$Index, 
  year_labels = Year_Set, 
  years_to_plot = which(fit$year_labels != 2020),
  strata_names = strata_names, 
  DirName = paste0(workDir, "results_age/proportions/"))

prop <- t(data.frame(proportions$Prop_ctl))
colnames(prop) <- c(paste0("age_", seq(from = 1, length = ncol(prop) - 1 )),
                    paste0("age_", ncol(prop), "+"))
rownames(prop) <- as.vector(sapply(X = strata_names, 
                                   FUN = function(x) paste0(x, "_", Year_Set)))

saveRDS(object = proportions, 
        file = paste0(workDir, "results_age/proportionsVAST_proportions.RDS"))
write.csv(x = prop, 
          file = paste0(workDir, "results_age/proportions/clean_proportions.csv"))


# Plot design-based estimate vs VAST estimate ------------------------------

## load model-based index ----

# If starting from Index.csv ...
# mb <- mb %>% filter(Stratum != "Both") %>%
#   mutate(age = Category - 1) %>%
#   select(stratum = Stratum,
#          year = Time,
#          age,
#          estimate_n = Estimate) %>%
#   mutate(estimator = "mb") %>%
#   group_by(stratum, year) %>%
#   mutate(N = sum(estimate_n)) %>%
#   group_by(age) %>%
#   mutate(proportion = estimate_n / N) %>%
#   ungroup()

mb <- read.csv(paste0(workDir, "results_age/proportions/clean_proportions.csv")) 

# elongate dataset to fields YEAR (year), recode vast ages to start at 0,
# AREA_ID (Both, EBS, or NBS), AGE (age), PROPORTION (proportion)
names(x = mb)[-1] <- 
  na.omit(as.numeric(gsub("[^0-9]", "", names(x = mb)[-1])) - 1)
mb$YEAR <- as.numeric(sub(".*_", "", mb$X))
mb$AREA_ID <- sub("_.*", "", mb$X)

mb <- 
  tidyr::pivot_longer(data = subset(x = mb, select = -X), 
               cols = c(paste(0:12)), 
               names_to = "AGE", 
               values_to = "PROPORTION")
mb$AGE <- as.numeric(x = mb$AGE)
head(mb)

## load design-based (db) estimated proportions from GAP_PRODUCTS.AGECOMP ----
chl <- gapindex::get_connected() # enter credentials in pop-out window

## Or from raw data, old way... 
# db_dat <- RODBC::sqlQuery(channel = chl, 
#                        query = "
# WITH FILTERED_STRATA AS (
# SELECT 
# AREA_ID,
# DESCRIPTION 
# FROM GAP_PRODUCTS.AKFIN_AREA
# WHERE AREA_TYPE = 'REGION'
# AND SURVEY_DEFINITION_ID IN (143, 98))
# 
# -- Select columns for output data
# SELECT 
# AGECOMP.AGE, 
# AGECOMP.YEAR,
# AGECOMP.POPULATION_COUNT, 
# AGECOMP.SEX, 
# DESCRIPTION
# 
# -- Identify what tables to pull data from
# FROM GAP_PRODUCTS.AKFIN_AGECOMP AGECOMP
# JOIN FILTERED_STRATA STRATA 
# ON STRATA.AREA_ID = AGECOMP.AREA_ID
# 
# -- Filter data results
# WHERE SPECIES_CODE = 21720")
# 
# db <- db_dat %>% 
#   mutate(stratum = recode(DESCRIPTION, 
#                           `EBS Standard Plus NW Region: All Strata` = "EBS", 
#                           `NBS Region: All Strata` = "NBS"),
#          age = AGE) %>% 
#   filter(stratum != "EBS Standard Region: All Strata",
#          age != -9,
#          YEAR > 1993) %>%
#   select(stratum,
#          year = YEAR, 
#          age, 
#          sex = SEX,
#          estimate_n = POPULATION_COUNT) %>%
#   mutate(estimator = "db") %>%
#   group_by(stratum, year, age) %>%
#   mutate(estimate_n = sum(estimate_n)) %>%
#   group_by(stratum, year) %>%
#   mutate(N = sum(estimate_n)) %>%
#   group_by(age) %>%
#   mutate(proportion = estimate_n / N) %>%
#   ungroup()

db <- RODBC::sqlQuery(channel = chl,
                                        query = "
WITH 

-- AGGREGATE AGECOMPS ACROSS NBS AND EBS WITH THE PLUS GROUP
AGE_ALL_W_PLUSGROUP AS (
SELECT YEAR, 
CASE
    WHEN AGE >= 12 THEN 12 -- PLUS GROUP
    ELSE AGE
END AS AGE,
SUM(POPULATION_COUNT) as POPULATION_COUNT
FROM GAP_PRODUCTS.AGECOMP 
WHERE SPECIES_CODE = 21720 -- PACIFIC COD
AND AREA_ID in (99900, -- EBS STANDARD + NW AREA
99902                  -- NBS
)
AND AGE >= 0

GROUP BY YEAR, 
CASE
    WHEN AGE >= 12 THEN 12 -- PLUS GROUP
    ELSE AGE
END
),

-- AGGREGATE AGECOMPS FOR THE NBS AND EBS (SEPARATELY) WITH THE PLUS GROUP
AGE_REGION_W_PLUSGROUP AS (
SELECT YEAR, AREA_ID,
CASE
    WHEN AGE >= 12 THEN 12 -- PLUS GROUP
    ELSE AGE
END AS AGE,
SUM(POPULATION_COUNT) as POPULATION_COUNT
FROM GAP_PRODUCTS.AGECOMP 
WHERE SPECIES_CODE = 21720 -- PACIFIC COD
AND AREA_ID in (99900, -- EBS STANDARD + NW AREA
99902                  -- NBS
)
AND AGE >= 0

GROUP BY YEAR, AREA_ID,
CASE
    WHEN AGE >= 12 THEN 12 -- PLUS GROUP
    ELSE AGE
END
)

SELECT YEAR, 'Both' AS AREA_ID, AGE, 
ROUND(POPULATION_COUNT * 1.0 / SUM(POPULATION_COUNT) OVER (PARTITION BY YEAR), 7) AS PROPORTION
FROM AGE_ALL_W_PLUSGROUP

UNION

SELECT YEAR, CASE
    WHEN AREA_ID = 99900 THEN 'EBS'
    ELSE 'NBS'
END AS AREA_ID, 
AGE, 
ROUND(POPULATION_COUNT * 1.0 / SUM(POPULATION_COUNT) OVER (PARTITION BY AREA_ID, YEAR), 7) AS PROPORTION
FROM AGE_REGION_W_PLUSGROUP

ORDER BY AREA_ID, YEAR, AGE                                 
")

# filter to years > 1994 for cod specifically
db <- filter(db, YEAR >= 1994)

## Query years of age data for each region
ebs_years <- 
  unique(x = db$YEAR[db$AREA_ID == "EBS"])
nbs_years <- 
  unique(x = db$YEAR[db$AREA_ID == "NBS"])

##   Merge design- and model-based age comps ----
agecomp_merge <- 
  merge(x = db,
        y = mb, all.y = TRUE,
        by = c("YEAR", "AREA_ID", "AGE"), 
        suffixes = c("_DBE", "_MBE"), sort = TRUE)

# ## combine indices and plot old way
# d <- bind_rows(db, mb) %>% filter(year != 2020)
# saveRDS(d, file=paste0(workDir,"results_age/db_mb_index.RDS"))
# 
# library(ggplot2)
# ggplot(filter(d, year %in% 2021:2023, stratum == "EBS"), 
#        aes(age, proportion, group = estimator, color = estimator)) +
#   geom_point(position = position_dodge(width = 0.5)) + 
#   facet_grid(year ~ . ) +
#   theme_bw()
# ggsave(file=paste0(workDir,"results_age/db_mb_index_comparison_ebs.png"), 
#        height = 8, width = 5, units = c("in"))
# 
# ggplot(filter(d, year %in% 2021:2023, stratum == "NBS"), 
#        aes(age, proportion, group = estimator, color = estimator)) +
#   geom_point(position = position_dodge(width = 0.5)) + 
#   facet_grid(year ~ . ) +
#   theme_bw()
# ggsave(file=paste0(workDir,"results_age/db_mb_index_comparison_nbs.png"), 
#        height = 8, width = 5, units = c("in"))


##   Plot age composition comparisons for each year with survey data ----
## EBS
pdf(file = paste0(workDir, "results_age/ebsagecomp_comps.pdf"),
    onefile = TRUE, width = 6, height = 6)
par(mfrow = c(3, 3), mar = c(2, 2, 1, 1), oma = c(2, 3, 0, 0))
for (iyear in ebs_years) {
  temp_df <- subset(agecomp_merge, YEAR == iyear & AREA_ID == "EBS")
  temp_df <- temp_df[order(temp_df$AGE), ]
  
  plot(PROPORTION_DBE ~ AGE, ann = F,
       data = temp_df,
       las = 1, pch = 16, col = "black", ylim = c(0, 0.57))
  legend("topright", legend = c("DB", "MB"), col = c("black", "red"),
         lty = 1, pch = 16)
  legend("topleft", legend = iyear, bty = 'n')
  lines(PROPORTION_DBE ~ AGE,
        data = temp_df, col = "black")
  points(PROPORTION_MBE ~ AGE,
         data = temp_df,
         pch = 16, col = "red")
  lines(PROPORTION_MBE ~ AGE,
        data = temp_df, col = "red")
  
  if (iyear %in% ebs_years[c(1, 10, 19, 28)]) {
    mtext(side = 1, text = "Age (years)", outer = TRUE, line = 0.5, font = 2)
    mtext(side = 2, text = "Proportion", outer = TRUE, line = 1, font = 2)
  }
}
dev.off()

## NBS
pdf(file = paste0(workDir, "results_age/nbsagecomp_comps.pdf"),
    onefile = TRUE, width = 5, height = 5)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1), oma = c(2, 3, 0, 0))
for (iyear in nbs_years) {
  temp_df <- subset(agecomp_merge, YEAR == iyear & AREA_ID == "NBS")
  temp_df <- temp_df[order(temp_df$AGE), ]
  
  plot(PROPORTION_DBE ~ AGE, ann = F,
       data = temp_df, main = iyear, las = 1, 
       pch = 16, col = "black", ylim = c(0, 0.57))
  legend("topright", legend = c("DB", "MB"), col = c("black", "red"),
         lty = 1, pch = 16)
  legend("topleft", legend = iyear, bty = 'n')
  lines(PROPORTION_DBE ~ AGE,
        data = temp_df, col = "black")
  points(PROPORTION_MBE ~ AGE,
         data = temp_df,
         pch = 16, col = "red")
  lines(PROPORTION_MBE ~ AGE,
        data = temp_df, col = "red")
}
mtext(side = 1, text = "Age (years)", outer = TRUE, line = 0.5, font = 2)
mtext(side = 2, text = "Proportion", outer = TRUE, line = 1, font = 2)
dev.off()

## Both Regions
pdf(file = paste0(workDir, "results_age/bsagecomp_comps.pdf"),
    onefile = TRUE, width = 5, height = 5)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1), oma = c(2, 3, 0, 0))
for (iyear in nbs_years) {
  temp_df <- subset(agecomp_merge, YEAR == iyear & AREA_ID == "Both")
  temp_df <- temp_df[order(temp_df$AGE), ]
  
  plot(PROPORTION_DBE ~ AGE, ann = F,
       data = temp_df, main = iyear, las = 1, 
       pch = 16, col = "black", ylim = c(0, 0.57))
  legend("topright", legend = c("DB", "MB"), col = c("black", "red"),
         lty = 1, pch = 16)
  legend("topleft", legend = iyear, bty = 'n')
  lines(PROPORTION_DBE ~ AGE,
        data = temp_df, col = "black")
  points(PROPORTION_MBE ~ AGE,
         data = temp_df,
         pch = 16, col = "red")
  lines(PROPORTION_MBE ~ AGE,
        data = temp_df, col = "red")
}
mtext(side = 1, text = "Age (years)", outer = TRUE, line = 0.5, font = 2)
mtext(side = 2, text = "Proportion", outer = TRUE, line = 1, font = 2)
dev.off()

# with(with(RODBC::sqlQuery(channel = chl,
#                           query = "
# select * from gap_products.agecomp
# where species_code = 21720
# and year = 2010
# and area_id = 99900
# and age > 0"),
#           
#           aggregate(POPULATION_COUNT ~ AGE,
#                     FUN = sum)),
#      POPULATION_COUNT / sum(POPULATION_COUNT))