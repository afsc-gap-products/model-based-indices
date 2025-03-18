##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       2023 GOA ESP VAST Hindcast compariosn (Pcod and pollock)
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Compare center of gravity (rotated) and effective area 
##                occupied between the 2022 production and 2023 hindcast 
##                VAST runs.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set up directories for a given species.
##   Download 2022 production run and create a folder on the 
##   Desktop for storage.
##   Set plot layout.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfcol = c(3, 2), mar = c(2, 5, 1, 1), oma = c(2.5, 0, 2, 0))
species_list <- c("Pacific cod" = "Gadus_macrocephalus",
                  "walleye pollock" = "Gadus_chalcogrammus")

for (ispp in 1:length(species_list)) {
  species <- species_list[ispp]
  dir_2022 <- paste0("C:/Users/zack.oyafuso/Desktop/Production_2022/", 
                     names(species), " rotated axis/")
  dir_2023 <- paste0("species_specific_code/GOA/Pcod_WPlk_ESP/", 
                     species, "_w140/")
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Pull index so we can query which years had data.
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cod_2022_index <- subset(x = read.csv(file = paste0(dir_2022, 
                                                      "Table_for_SS3.csv")), 
                           subset = Estimate_metric_tons > 0)
  
  years <- cod_2022_index$Year
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Compare Center of Gravity
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cog_2022 <- subset(x = read.csv(file = paste0(dir_2022, "COG.csv")),
                     subset = m %in% 3:4 & Year %in% years)
  cog_2023 <- subset(x = read.csv(file = paste0(dir_2023,
                                                "results/output/COG.csv")),
                     subset = m %in% 3:4 & Year %in% years)
  
  ## Plot COG on the rotated x-axis
  plot(COG_hat ~ Year, subset = m == 3, data = cog_2022,
       ylim = range(c((cog_2022$COG_hat-cog_2022$SE)[cog_2022$m == 3],
                      (cog_2022$COG_hat+cog_2022$SE)[cog_2022$m == 3])),
       ylab = "Rotated Eastings (km)", pch = 16, col = "black", las = 1)
  lines(COG_hat ~ Year, subset = m == 3, data = cog_2022)
  with(subset(x = cog_2022, subset = m == 3),
       segments(x0 = Year, x1 = Year, y0 = COG_hat - SE, y1 = COG_hat + SE))
  
  
  points(COG_hat ~ I(Year+0.5), subset = m == 3, data = cog_2023, pch = 16, col = "red")
  lines(COG_hat ~ I(Year+0.5), subset = m == 3, data = cog_2023, col = "red")
  with(subset(x = cog_2023, subset = m == 3),
       segments(x0 = Year + 0.5, x1 = Year + 0.5, 
                y0 = COG_hat - SE, y1 = COG_hat + SE,
                col = "red"))
  
  mtext(side = 3, text = names(species), font = 2, line = 1)
  
  ## Plot COG on rotated y-axis
  plot(COG_hat ~ Year, subset = m == 4, data = cog_2022,
       ylab = "Rotated Northings (km)", 
       ylim = range(c((cog_2022$COG_hat-cog_2022$SE)[cog_2022$m == 4],
                      (cog_2022$COG_hat+cog_2022$SE)[cog_2022$m == 4])),
       pch = 16, col = "black", las = 1)
  lines(COG_hat ~ Year, subset = m == 4, data = cog_2022)
  with(subset(x = cog_2022, subset = m == 4),
       segments(x0 = Year, x1 = Year, y0 = COG_hat - SE, y1 = COG_hat + SE))
  
  points(COG_hat ~ I(Year + 0.5), 
         subset = m == 4, data = cog_2023, 
         pch = 16, col = "red")
  lines(COG_hat ~ I(Year + 0.5), 
        subset = m == 4, data = cog_2023, 
        col = "red")
  with(subset(x = cog_2023, subset = m == 4),
       segments(x0 = Year + 0.5, x1 = Year + 0.5, 
                y0 = COG_hat - SE, y1 = COG_hat + SE,
                col = "red"))
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Effective Area Occupied
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eao_2022 <- subset(x = read.csv(file = paste0(dir_2022, 
                                                "ln_effective_area.csv")),
                     subset = year %in% years)
  eao_2023 <- subset(x = read.csv(file = paste0(dir_2023, "results/output/",
                                                "ln_effective_area.csv")),
                     subset = year %in% years)
  
  plot(Estimate ~ year, data = eao_2022, 
       ylim = range(c(eao_2022$Estimate - eao_2022$Std..Error,
                      eao_2022$Estimate + eao_2022$Std..Error, 
                      eao_2023$Estimate - eao_2023$Std..Error,
                      eao_2023$Estimate + eao_2023$Std..Error)),
       ylab = "ln(Area)", pch = 16, col = "black", las = 1)
  lines(Estimate ~ year, data = eao_2022)
  segments(x0 = eao_2022$year, x1 = eao_2022$year,
           y0 = eao_2022$Estimate - eao_2022$Std..Error,
           y1 = eao_2022$Estimate + eao_2022$Std..Error)
  
  points(Estimate ~ I(year+0.5), data = eao_2023, pch = 16, col = "red")
  lines(Estimate ~ I(year+0.5), data = eao_2023, col = "red")
  segments(x0 = eao_2023$year + 0.5, x1 = eao_2023$year + 0.5,
           y0 = eao_2023$Estimate - eao_2023$Std..Error,
           y1 = eao_2023$Estimate + eao_2023$Std..Error, 
           col = "red")
  
  if (species == "Gadus_macrocephalus")
    legend("bottomleft", 
           legend = c("2022 Production (+/- SE)", "2023 Hindcast (+/- SE)"), 
           pch = 16, lty = 1, col = c("black", "red"), bty = "n")
}

mtext(side = 1, text = "Year", outer = TRUE, line = 1, font = 2)
