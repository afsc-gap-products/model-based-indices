# Compare COG and EAO from 2023 (VAST w/ bias correction) and 2025 (sdmTMB
# without bias correction). Make sure to donwload the 2023 COGs and EAO from
# the google drive folder.
# https://drive.google.com/drive/folders/1pzSxThFPRcxYDyV1Ih-mX6L77sXZKkMQ
# Author: Zack Oyafuso

species <- c("Gadus_chalcogrammus", "Gadus_macrocephalus", "Atheresthes_stomias")

for (species in c("Gadus_chalcogrammus", "Gadus_macrocephalus", "Atheresthes_stomias")) {
  dir <- paste0("species_specific_code/GOA/", species, "/production/")
  
  ## Import the COG and EAO from 2025 (using sdmTMB) and 2023 (using VAST)
  cog_sdmtmb <- read.csv(file = paste0(dir, "cog.csv"))
  cog_vast <- read.csv(file = paste0(dir, "COG_vast.csv"))
  
  eao_sdmtmb <- read.csv(file = paste0(dir, "area_occupied.csv"))
  eao_vast <- read.csv(file = paste0(dir, "ln_effective_area_vast.csv"))
  names(x = eao_vast) <- c("Estimate", "SE", "year")
  
  pdf(file = paste0(dir, "cog_eao_bridging.pdf"), width = 6, height = 8)
  ## Plot COG (Eastings)
  par(mfrow = c(3, 1), mar = c(4, 4, 1, 1), oma = c(1, 1, 2, 0))
  plot(est_x ~ year, data = cog_sdmtmb, type = "n", las = 1, 
       xlab = "Year", ylab = "COG Eastings (km)",
       ylim = c(min(est_x-2*se_x), 1.1 * max(est_x+2*se_x)))
  polygon(x = c(cog_sdmtmb$year, rev(cog_sdmtmb$year)),
          y = c(cog_sdmtmb$est_x - 2*cog_sdmtmb$se_x, 
                rev(cog_sdmtmb$est_x + 2*cog_sdmtmb$se_x)),
          col = rgb(255,160,122, maxColorValue=255, alpha = 100), border = F)
  points(est_x ~ year, data = cog_sdmtmb, pch = 16, cex = 0.5,
         col = rgb(250,128,114, maxColorValue=255))
  lines(est_x ~ year, data = cog_sdmtmb,
        col = rgb(250,128,114, maxColorValue=255))
  legend("topright", legend = c("VAST (as of 2023)", "sdmTMB (as of 2025)"),
         fill = c(rgb(64,224,208, maxColorValue=255, alpha = 100),
                  rgb(255, 160, 122, maxColorValue = 255, alpha = 100)),
         cex = 1.2)
  
  with(subset(x = cog_vast,
              subset = m == 1 & Year %in% cog_sdmtmb$year),
       polygon(x = c(Year, rev(Year)),
               y = c(COG_hat - 2*SE, 
                     rev(COG_hat + 2*SE)),
               col = rgb(64,224,208, maxColorValue=255, alpha = 100), 
               border = F)
  )
  
  
  points(COG_hat ~ Year, pch = 16, cex = 0.5,
         data = cog_vast,
         subset = m == 1 & Year %in% cog_sdmtmb$year,
         col = rgb(0,128,128, maxColorValue=255))
  lines(COG_hat ~ Year, pch = 16, cex = 0.5,
        data = cog_vast,
        subset = m == 1 & Year %in% cog_sdmtmb$year,
        col = rgb(0,128,128, maxColorValue=255))
  
  ## Plot COG (Northings)
  plot(est_y ~ year, data = cog_sdmtmb, type = "n", las = 1, 
       xlab = "Year", ylab = "COG Northings (km)",
       ylim = c(min(est_y-2*se_y), max(est_y+2*se_y)))
  polygon(x = c(cog_sdmtmb$year, rev(cog_sdmtmb$year)),
          y = c(cog_sdmtmb$est_y - 2*cog_sdmtmb$se_y, 
                rev(cog_sdmtmb$est_y + 2*cog_sdmtmb$se_y)),
          col = rgb(255, 160, 122, maxColorValue = 255, alpha = 100), 
          border = F)
  points(est_y ~ year, data = cog_sdmtmb, pch = 16, cex = 0.5,
         col = rgb(250, 128, 114, maxColorValue=255))
  lines(est_y ~ year, data = cog_sdmtmb,
        col = rgb(250, 128, 114, maxColorValue=255))
  
  with(subset(x = cog_vast,
              subset = m == 2 & Year %in% cog_sdmtmb$year),
       polygon(x = c(Year, rev(Year)),
               y = c(COG_hat - 2*SE, 
                     rev(COG_hat + 2*SE)),
               col = rgb(64,224,208, maxColorValue=255, alpha = 100), border = F)
  )
  
  
  points(COG_hat ~ Year, pch = 16, cex = 0.5,
         data = cog_vast,
         subset = m == 2 & Year %in% cog_sdmtmb$year,
         col = rgb(0, 128, 128, maxColorValue=255))
  lines(COG_hat ~ Year, pch = 16, cex = 0.5,
        data = cog_vast,
        subset = m == 2 & Year %in% cog_sdmtmb$year,
        col = rgb(0, 128, 128, maxColorValue=255))
  
  ## Plot Effective Area Occupied
  plot(log_est ~ year, data = eao_sdmtmb, type = "n", las = 1, 
       xlab = "Year", ylab = "Effective Area Occupied (log-km2)",
       ylim = c(min(log_est-2*se), max(log_est+2*se)))
  polygon(x = c(eao_sdmtmb$year, rev(eao_sdmtmb$year)),
          y = c(eao_sdmtmb$log_est - 2*eao_sdmtmb$se, 
                rev(eao_sdmtmb$log_est + 2*eao_sdmtmb$se)),
          col = rgb(255,160,122, maxColorValue=255, alpha = 100), 
          border = F)
  points(log_est ~ year, data = eao_sdmtmb, pch = 16, cex = 0.5,
         col = rgb(250,128,114, maxColorValue=255))
  lines(log_est ~ year, data = eao_sdmtmb,
        col = rgb(250,128,114, maxColorValue=255))
  
  with(eao_vast,
       polygon(x = c(year, rev(year)),
               y = c(Estimate  - 2*SE, 
                     rev(Estimate  + 2*SE)),
               col = rgb(64,224,208, maxColorValue=255, alpha = 100), border = F)
  )
  
  
  points(Estimate ~ year, pch = 16, cex = 0.5,
         data = eao_vast,
         col = rgb(0,128,128, maxColorValue=255))
  lines(Estimate ~ year, pch = 16, cex = 0.5,
        data = eao_vast,
        col = rgb(0,128,128, maxColorValue=255))
  
  mtext(side = 3, outer = TRUE, 
        text = gsub(x = species, 
                    pattern = "_", 
                    replacement = " "))
  dev.off()
}
