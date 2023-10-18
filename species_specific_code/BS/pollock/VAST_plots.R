# Custom plots of VAST results for pollock
# By: Sophia Wassermann
# By: Caitlin I. Allen Akselrud
# Contact: sophia.wassermann@noaa.gov
# Date created: 2023.10.17

library(here)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

### VAST index ----------------------------------------------------------------
index <- read.csv("Index.csv")
colnames(index)[6] <- "error"  # easier column name for plotting

index_all_areas <- ggplot(index, aes(x = Time, y = (Estimate / 1000000000))) +
  geom_line(linetype = "dashed") +
  geom_pointrange(aes(ymin = (Estimate / 1000000000) - (error / 1000000000),
                      ymax = (Estimate / 1000000000) + (error / 1000000000)), alpha = 0.8) +
  scale_shape(solid = FALSE) +
  xlab("Year") + ylab("Index (Mt)") +
  facet_wrap(~ Stratum, ncol = 1)
index_all_areas

ggsave(index_all_areas, filename = here("species_specific_code", "BS", "pollock", "plots", "index_all_areas.png"),
       width=130, height=160, units="mm", dpi=300)

# Plot just EBS
index_ebs <- ggplot(index %>% filter(Stratum == "EBS"), aes(x = Time, y = (Estimate / 1000000000))) +
  geom_line(linetype = "dashed") +
  geom_pointrange(aes(ymin = (Estimate / 1000000000) - (error / 1000000000),
                      ymax = (Estimate / 1000000000) + (error / 1000000000)), alpha = 0.8) +
  scale_shape(solid = FALSE) +
  xlab("Year") + ylab("Index (Mt)") 
index_ebs