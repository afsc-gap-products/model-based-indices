index_comp <- 
  read.csv(file = "species_specific_code/BS/yellowfin_sole/dbe_vs_vast_comp.csv")
vast_results <- 
  readRDS(file = "species_specific_code/BS/yellowfin_sole/VASTresults.RDS")
vast_fit <- 
  readRDS(file = "species_specific_code/BS/yellowfin_sole/yellowfin_sole_VASTfit.RDS")

{
  png(file = "species_specific_code/BS/yellowfin_sole/production/bsindex_comps.png",
      width = 8, height = 6, units = "in", res = 500)
  
  par(mfrow = c(2, 2), mar = c(2.5, 5, 3, 0.5), oma = c(2, 0, 0, 0))
  plot(BIOMASS_MT_VAST/1e6 ~ YEAR, las = 1, 
       xlab = "", ylab = "Biomass (million mt)",
       subset = REGION == "EBS" & YEAR != 2020, 
       xlim = c(1980, max(index_comp$YEAR)), ylim = c(1, 4), 
       data = index_comp, main = "EBS Biomass Index",
       col = "red", pch = 16)
  lines(BIOMASS_MT_VAST/1e6 ~ YEAR,
        subset = REGION == "EBS" & YEAR != 2020,
        data = index_comp,
        col = "red")
  
  points(BIOMASS_MT_DBE/1e6 ~ YEAR,
         subset = REGION == "EBS" & YEAR != 2020,
         data = index_comp, pch = 16,
         col = "black")
  lines(BIOMASS_MT_DBE/1e6 ~ YEAR,
        subset = REGION == "EBS" & YEAR != 2020,
        data = index_comp,
        col = "black")
  
  plot(BIOMASS_CV_VAST ~ YEAR, las = 1, 
       main = "Coefficient of Variation\nof EBS Biomass Index",
       subset = REGION == "EBS" & YEAR != 2020, 
       xlim = c(1980, max(index_comp$YEAR)), ylim = c(0, 0.15),
       data = index_comp, xlab = "", ylab = "Coefficient of Variation",
       col = "red", pch = 16)
  lines(BIOMASS_CV_VAST ~ YEAR,
        subset = REGION == "EBS" & YEAR != 2020,
        data = index_comp,
        col = "red")
  
  points(BIOMASS_CV_DBE ~ YEAR,
         subset = REGION == "EBS" & YEAR != 2020,
         data = index_comp, pch = 16,
         col = "black")
  lines(BIOMASS_CV_DBE ~ YEAR,
        subset = REGION == "EBS" & YEAR != 2020,
        data = index_comp,
        col = "black")
  
  plot(BIOMASS_MT_VAST/1e6 ~ YEAR, las = 1, 
       xlab = "", ylab = "Biomass (million mt)",
       subset = REGION == "NBS" & YEAR != 2020, 
       xlim = c(1980, max(index_comp$YEAR)), ylim = c(0, 1),
       data = index_comp, main = "NBS Biomass Index",
       col = "red", pch = 16)
  lines(BIOMASS_MT_VAST/1e6 ~ YEAR,
        subset = REGION == "NBS" & YEAR != 2020,
        data = index_comp,
        col = "red")
  
  points(BIOMASS_MT_DBE/1e6 ~ YEAR,
         subset = REGION == "NBS" & YEAR != 2020,
         data = index_comp, pch = 16,
         col = "black")
  lines(BIOMASS_MT_DBE/1e6 ~ YEAR,
        subset = REGION == "NBS" & YEAR != 2020,
        data = index_comp,
        col = "black")
  
  plot(BIOMASS_CV_VAST ~ YEAR, las = 1,
       main = "Coefficient of Variation\nof NBS Biomass Index",
       subset = REGION == "NBS" & YEAR != 2020, 
       xlim = c(1980, max(index_comp$YEAR)), ylim = c(0, 0.7),
       data = index_comp, xlab = "", ylab = "Coefficient of Variation",
       col = "red", pch = 16)
  lines(BIOMASS_CV_VAST ~ YEAR,
        subset = REGION == "NBS" & YEAR != 2020,
        data = index_comp,
        col = "red")
  
  points(BIOMASS_CV_DBE ~ YEAR,
         subset = REGION == "NBS" & YEAR != 2020,
         data = index_comp, pch = 16,
         col = "black")
  lines(BIOMASS_CV_DBE ~ YEAR,
        subset = REGION == "NBS" & YEAR != 2020,
        data = index_comp,
        col = "black")
  
  mtext(side = 1, text = "Year", outer = TRUE, line = 0.5, font = 2)
  
  dev.off()
  }
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Center of Gravity
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cog <- as.data.frame(vast_results$Range$COG_Table)
cog$Year <- as.numeric(cog$Year)
cog$m <- as.numeric(cog$m)
cog$COG_hat <- as.numeric(cog$COG_hat)
cog$SE <- as.numeric(cog$SE)

png(file = "species_specific_code/BS/yellowfin_sole/production/cog.png",
    width = 8, height = 3, res = 500, units = "in")
par(mfrow = c(1, 2), mar = c(2, 4, 1, 1), oma = c(1, 0, 0, 0))
for (m_axis in 1:2) {
  se_boundaries <- c((cog$COG_hat + cog$SE)[cog$m == m_axis],
                     rev((cog$COG_hat - cog$SE)[cog$m == m_axis]))
  plot(COG_hat ~ Year, 
       data = cog, 
       ylab = c("Eastings (km)", "Northings (km)")[m_axis], xlab = "",
       xlim = c(1980, max(cog$Year)),
       ylim = range(se_boundaries), 
       subset = m == m_axis, 
       las = 1, type = "b", pch = 16, col = "red", cex = 0.5)
  lines(COG_hat ~ Year, data = cog, subset = m == m_axis, col = "red")
  polygon(x = c(unique(cog$Year), rev(unique(cog$Year))),
          y = se_boundaries,
          col = rgb(1,0,0,0.2),
          border = F)
}
mtext(side = 1, text = "Year", outer = TRUE)
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

png(file = "species_specific_code/BS/yellowfin_sole/production/effective_area.png",
    width = 5, height = 5, res = 500, units = "in")
effective_area <- vast_results$Range$SD_effective_area_ctl[1, , , ]
par(mfrow = c(1, 1), mar = c(4, 4, 1, 1))
matplot(x = sort(x = unique(cog$Year)),
        y = effective_area[, , "Estimate"] / 1000,
        xlab = "Year",
        ylab = "Effective Area Occupied (thousand square km)",
        xlim = c(1980, max(cog$Year)),
        ylim = c(0, max(effective_area[, , "Estimate"] +
                          effective_area[, , "Std. Error"])) / 1000 ,
        col = c("red", "darkgreen", "blue"), pch = 1, las = 1, cex = 0.75)
matlines(x = sort(x = unique(cog$Year)),
         y = effective_area[, , "Estimate"] / 1000,
         col = c("red", "darkgreen", "blue"), pch = 16, lty = 1)
segments(x0 = sort(x = unique(cog$Year)),
         y0 = (effective_area[, 1, "Estimate"] - 
                 effective_area[, 1, "Std. Error"]) / 1000,
         y1 = (effective_area[, 1, "Estimate"] + 
                 effective_area[, 1, "Std. Error"]) / 1000, 
         col = "red")
segments(x0 = sort(x = unique(cog$Year)),
         y0 = (effective_area[, 2, "Estimate"] - 
                 effective_area[, 2, "Std. Error"])/ 1000,
         y1 = (effective_area[, 2, "Estimate"] + 
                 effective_area[, 2, "Std. Error"])/ 1000, col = "darkgreen")
segments(x0 = sort(x = unique(cog$Year)),
         y0 = (effective_area[, 3, "Estimate"] - 
                 effective_area[, 3, "Std. Error"])/ 1000,
         y1 = (effective_area[, 3, "Estimate"] + 
                 effective_area[, 3, "Std. Error"])/ 1000, col = "blue")
legend("bottomright", legend = c("EBS + NBS", "EBS", "NBS"), 
       col = c("red", "darkgreen", "blue"), pch = 1, lty = 1, ncol = 1, bty = "n")
dev.off()
