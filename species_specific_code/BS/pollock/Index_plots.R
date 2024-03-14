# Custom plots of VAST index of abundance results for Bering sea pollock
# By: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2023.10.17
# Date updated: 2024.03.14

library(here)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(reshape2)
library(gapindex)

# Set ggplot theme
devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
theme_set(theme_sleek())

this_year <- 2023

### VAST index ----------------------------------------------------------------
index <- read.csv("Index.csv")

# Calculate 2023 difference for both areas
diff_23 <- ((index[42, 5] - index[41, 5]) / index[41, 5]) * 100
mean_23 <- (index[42, 5] / mean((index %>% filter(Stratum == "EBS"))[, 5])) * 100

# Calculate 2023 difference for EBS
ebs_23_diff <- ((index[84, 5] - index[83, 5]) / index[83, 5]) * 100
ebs_23_mean <- (index[84, 5] / mean((index %>% filter(Stratum == "EBS"))[, 5])) * 100

# Caclulate percent NBS for 2023
nbs_23 <- (index[126, 5] / index[42, 5]) * 100

colnames(index)[6] <- "error"  # easier column name for plotting

index_all_areas <- ggplot(index %>% filter(Time != 2020), 
                          aes(x = Time, y = (Estimate / 1000000000))) +
  geom_line(alpha = 0.4) +
  geom_pointrange(aes(ymin = (Estimate / 1000000000) - (error / 1000000000),
                      ymax = (Estimate / 1000000000) + (error / 1000000000)), alpha = 0.8) +
  ylim(0, NA) +
  xlab("Year") + ylab("Index (Mt)") +
  facet_wrap(~ Stratum, ncol = 1)
index_all_areas

ggsave(index_all_areas, filename = here("species_specific_code", "BS", "pollock", "plots", "index_all_areas.png"),
       width=130, height=160, units="mm", dpi=300)

# Plot just EBS
index_ebs <- ggplot(index %>% filter(Stratum == "EBS" & Time != 2020), aes(x = Time, y = (Estimate / 1000000000))) +
  geom_line(alpha = 0.4) +
  geom_pointrange(aes(ymin = (Estimate / 1000000000) - (error / 1000000000),
                      ymax = (Estimate / 1000000000) + (error / 1000000000)), alpha = 0.8) +
  ylim(0, NA) +
  xlab("Year") + ylab("Index (Mt)") 
index_ebs


### Cold pool extent covariate ------------------------------------------------
# Cold pool covariabe
cold_pool <- read.csv(here("output", "cold_pool_scaled_formatted.csv")) %>%
  mutate(extent = case_when(area_lte2_km2 >= 0 ~ "greater",
                            area_lte2_km2 < 0 ~ "less"))

cold_pool_plot <- ggplot() +
  geom_bar(data = cold_pool, aes(x = Year, y = area_lte2_km2, fill = extent), stat = "identity") +
  scale_fill_manual(values = c("cornflowerblue", "darkred")) +
  ylab("Cold pool covariate")
cold_pool_plot

ggsave(cold_pool, filename = here("output", "cold_pool_covariate.png"),
       width = 120, height = 100, unit = "mm", dpi = 300)


# Cold pool vs. index value ---------------------------------------------------
cold_pollock_cor <- cor(index[index$Stratum == "EBS", ]$Estimate, cold_pool$area_lte2_km2)

# Index relative to mean
mean_index <- mean(index[index$Stratum == "EBS", ]$Estimate)
relative_index <- cbind.data.frame(Year = index[index$Stratum == "EBS", ]$Time, 
                                   Index = (index[index$Stratum == "EBS", ]$Estimate - mean_index) /
                                     mean_index)
relative_index[relative_index$Year == 2020, ]$Index <- 0

cold_pollock_plot <- cold_pool_plot +
  geom_point(data = relative_index, aes(x = Year, y = Index), size = 2) +
  geom_line(data = relative_index, aes(x = Year, y = Index), size = 1)
cold_pollock_plot
