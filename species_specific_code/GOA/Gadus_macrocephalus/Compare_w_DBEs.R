library(RODBC)

# MBE <- read.csv("C:/Users/zack.oyafuso/Downloads/Index.csv")
# MBE <- subset(x = MBE, Estimate > 0)
# names(MBE) <- c("Category", "Time", "Stratum", "Units", 
#                 "Estimate", "SE", "CV")

source("R/get_connected.R")

DBE <- RODBC::sqlQuery(
  channel = channel, 
  query = paste("SELECT YEAR, SPECIES_CODE, TOTAL_BIOMASS, BIOMASS_VAR",
                "FROM GOA.BIOMASS_TOTAL WHERE",
                "SPECIES_CODE = 21720 AND",
                "YEAR >= 1990"))

DBE <- DBE[order(DBE$YEAR), ]
DBE$BIOMASS_SD <- sqrt(DBE$BIOMASS_VAR)
DBE$BIOMASS_CV <- DBE$BIOMASS_SD / DBE$TOTAL_BIOMASS

plot(TOTAL_BIOMASS/1000 ~ YEAR, data = DBE, type = "p", pch = 16,
     ylim = c(0, 1000), las = 1)
lines(TOTAL_BIOMASS/1000 ~ YEAR, data = DBE)
segments(x0 = DBE$YEAR,
         x1 = DBE$YEAR,
         y0 = DBE$TOTAL_BIOMASS/1000 - DBE$BIOMASS_SD/1000,
         y1 = DBE$TOTAL_BIOMASS/1000 + DBE$BIOMASS_SD/1000)

# points(Estimate/1e6 ~ Time,
#        data = MBE,
#        pch = 16,
#        col = "red")
# lines(Estimate/1e6 ~ Time,
#        data = MBE,
#        pch = 16,
#        col = "red")
# segments(x0 = MBE$Time,
#          x1 = MBE$Time,
#          y0 = MBE$Estimate/1e6 - MBE$SE/1e6,
#          y1 = MBE$Estimate/1e6 + MBE$SE/1e6,
#          col = "red")
