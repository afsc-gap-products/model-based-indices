index_comp <- 
  read.csv(file = "species_specific_code/BS/yellowfin_sole/dbe_vs_vast_comp.csv")

pdf(file = "species_specific_code/BS/yellowfin_sole/production/bsindex_comps.pdf",
    width = 8, height = 6)
par(mfrow = c(2, 2))
plot(BIOMASS_MT_VAST/1e6 ~ YEAR, las = 1, ylab = "Total Biomass (million mt)",
     subset = REGION == "EBS" & YEAR != 2020, yli = c(1, 4),
     data = index_comp, main = "EBS YFS Total Biomass",
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

plot(BIOMASS_CV_VAST ~ YEAR, las = 1, main = "CV of EBS YFS Total Biomass Estimate",
     subset = REGION == "EBS" & YEAR != 2020, yli = c(0, 0.15),
     data = index_comp, ylab = "CV",
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

plot(BIOMASS_MT_VAST/1e6 ~ YEAR, las = 1, ylab = "Total Biomass (million mt)",
     subset = REGION == "NBS" & YEAR != 2020, yli = c(0, 1),
     data = index_comp, main = "NBS YFS Total Biomass",
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

plot(BIOMASS_CV_VAST ~ YEAR, las = 1, main = "CV of NBS YFS Total Biomass Estimate",
     subset = REGION == "NBS" & YEAR != 2020, ylim = c(0, 0.7),
     data = index_comp, ylab = "CV",
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
dev.off()
