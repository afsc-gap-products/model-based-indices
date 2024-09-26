##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Age Composition VAST comparisons with design-based 
##                estiamtes of age composition
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Bering Sea yellowfin sole, 2024 production run
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import libraries, connect to Oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(gapindex)
library(tidyr)
chl <- gapindex::get_connected(check_access = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import VAST Age Compositions, elongate dataset to fields YEAR (year), 
##   AREA_ID (Both, EBS, or NBS), AGE (age), PROPORTION (proportion)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
yfs_bs_age_props_mbe <- 
  read.csv(file = "species_specific_code/BS/yellowfin_sole/clean_proportions.csv")
names(x = yfs_bs_age_props_mbe)[-1] <- 
  gsub("[^0-9]", "", names(x = yfs_bs_age_props_mbe)[-1])
yfs_bs_age_props_mbe$YEAR <- as.numeric(sub(".*_", "", yfs_bs_age_props_mbe$X))
yfs_bs_age_props_mbe$AREA_ID <- sub("_.*", "", yfs_bs_age_props_mbe$X)

yfs_bs_age_props_mbe <- 
  pivot_longer(data = subset(x = yfs_bs_age_props_mbe, select = -X), 
               cols = c(paste(1:20)), 
               names_to = "AGE", 
               values_to = "PROPORTION")
yfs_bs_age_props_mbe$AGE <- as.numeric(x = yfs_bs_age_props_mbe$AGE)
head(yfs_bs_age_props_mbe)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate design-based age proportions from GAP_PRODUCTS.AGECOMP 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
yfs_bs_age_props_dbe <- RODBC::sqlQuery(channel = chl,
                                        query = "
WITH 

-- AGGREGATE AGECOMPS ACROSS NBS AND EBS WITH THE PLUS GROUP
AGE_ALL_W_PLUSGROUP AS (
SELECT YEAR, 
CASE
    WHEN AGE >= 20 THEN 20 -- PLUS GROUP
    ELSE AGE
END AS AGE,
SUM(POPULATION_COUNT) as POPULATION_COUNT
FROM GAP_PRODUCTS.AGECOMP 
WHERE SPECIES_CODE = 10210 -- YELLOWFIN SOLE
AND AREA_ID in (99900, -- EBS STANDARD + NW AREA
99902                  -- NBS
)
AND AGE >= 0

GROUP BY YEAR, 
CASE
    WHEN AGE >= 20 THEN 20 -- PLUS GROUP
    ELSE AGE
END
),

-- AGGREGATE AGECOMPS FOR THE NBS AND EBS (SEPARATLEY) WITH THE PLUS GROUP
AGE_REGION_W_PLUSGROUP AS (
SELECT YEAR, AREA_ID,
CASE
    WHEN AGE >= 20 THEN 20 -- PLUS GROUP
    ELSE AGE
END AS AGE,
SUM(POPULATION_COUNT) as POPULATION_COUNT
FROM GAP_PRODUCTS.AGECOMP 
WHERE SPECIES_CODE = 10210 -- YELLOWFIN SOLE
AND AREA_ID in (99900, -- EBS STANDARD + NW AREA
99902                  -- NBS
)
AND AGE >= 0

GROUP BY YEAR, AREA_ID,
CASE
    WHEN AGE >= 20 THEN 20 -- PLUS GROUP
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

## Query years of age data for each region
ebs_years <- 
  unique(x = yfs_bs_age_props_dbe$YEAR[yfs_bs_age_props_dbe$AREA_ID == "EBS"])
nbs_years <- 
  unique(x = yfs_bs_age_props_dbe$YEAR[yfs_bs_age_props_dbe$AREA_ID == "NBS"])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge design- and model-based age comps
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
agecomp_merge <- 
  merge(x = yfs_bs_age_props_dbe,
        y = yfs_bs_age_props_mbe, all.y = TRUE,
        by = c("YEAR", "AREA_ID", "AGE"), 
        suffixes = c("_DBE", "_MBE"), sort = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot age composition comparisons for each year with survey data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## EBS
pdf(file = paste0("species_specific_code/BS/yellowfin_sole/",
                  "production/ebsagecomp_comps.pdf"),
    onefile = TRUE, width = 6, height = 6)
par(mfrow = c(3, 3), mar = c(2, 2, 1, 1), oma = c(2, 3, 0, 0))
for (iyear in ebs_years) {
  temp_df <- subset(agecomp_merge, YEAR == iyear & AREA_ID == "EBS")
  temp_df <- temp_df[order(temp_df$AGE), ]
  
  plot(PROPORTION_DBE ~ AGE, ann = F,
       data = temp_df,
       las = 1, pch = 16, col = "black", ylim = c(0, 0.4))
  legend("topright", legend = c("DBE", "VAST"), col = c("black", "red"),
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
    mtext(side = 2, text = "Propotion", outer = TRUE, line = 1, font = 2)
  }
}
dev.off()


## NBS
pdf(file = paste0("species_specific_code/BS/yellowfin_sole/",
                  "production/nbsagecomp_comps.pdf"),
    onefile = TRUE, width = 5, height = 5)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1), oma = c(2, 3, 0, 0))
for (iyear in nbs_years) {
  temp_df <- subset(agecomp_merge, YEAR == iyear & AREA_ID == "NBS")
  temp_df <- temp_df[order(temp_df$AGE), ]
  
  plot(PROPORTION_DBE ~ AGE, ann = F,
       data = temp_df, main = iyear, las = 1, 
       pch = 16, col = "black", ylim = c(0, 0.4))
  legend("topright", legend = c("DBE", "VAST"), col = c("black", "red"),
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
mtext(side = 2, text = "Propotion", outer = TRUE, line = 1, font = 2)
dev.off()

## Both Regions
pdf(file = paste0("species_specific_code/BS/yellowfin_sole/",
                  "production/bsagecomp_comps.pdf"),
    onefile = TRUE, width = 5, height = 5)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1), oma = c(2, 3, 0, 0))
for (iyear in nbs_years) {
  temp_df <- subset(agecomp_merge, YEAR == iyear & AREA_ID == "Both")
  temp_df <- temp_df[order(temp_df$AGE), ]
  
  plot(PROPORTION_DBE ~ AGE, ann = F,
       data = temp_df, main = iyear, las = 1, 
       pch = 16, col = "black", ylim = c(0, 0.4))
  legend("topright", legend = c("DBE", "VAST"), col = c("black", "red"),
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
mtext(side = 2, text = "Propotion", outer = TRUE, line = 1, font = 2)
dev.off()

with(with(RODBC::sqlQuery(channel = chl,
                query = "
select * from gap_products.agecomp
where species_code = 10210
and year = 2010
and area_id = 99900
and age > 0"),
     
     aggregate(POPULATION_COUNT ~ AGE,
               FUN = sum)),
     POPULATION_COUNT / sum(POPULATION_COUNT))
