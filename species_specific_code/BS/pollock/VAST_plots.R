# Custom plots of VAST results for pollock
# By: Sophia Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2023.10.17

library(here)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggsidekick)
# Set ggplot theme
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


### Compare VAST comps to hindcast --------------------------------------------
# Plot Difference between hindcast and production run -------------------------
new_props <- read.csv("proportions.csv")[, -1]
old_props <- read.csv("proportions_2022.csv")[, -1]

new_props <- subset(new_props, new_props$Year < 2023)

check_props <- round(new_props[,1:15] - old_props[,1:15], 4)
check_props_tab <- cbind(check_props, new_props[,16:17])
check_props_abs <- round(abs(new_props[,1:15] - old_props[,1:15]), 4)
check_props_abs_tab <-  cbind(check_props_abs, new_props[,16:17])
options(scipen=999)
write.csv(check_props_tab, here("VAST_results", "bridge_props_2023.csv"))
write.csv(check_props_abs_tab, here("VAST_results", "bridge_props_abs_2023.csv"))

colnames(check_props_tab)[1:15] <- 1:15
props_plot <- melt(check_props_tab, id.vars = c("Year", "Region"), 
                   variable.name = "Age", value.name = "Proportion") %>%
  # add column for coloring the bars in the plot based on positive/negative
  mutate(sign = case_when(Proportion >= 0 ~ "positive",
                          Proportion < 0 ~ "negative"))

# Plot both regions together and without 2020
comp_diff <- ggplot(props_plot %>% filter(Region == "EBS" & Year != 2020), 
                    aes(x = Age, y = Proportion, fill = sign)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = c("cornflowerblue", "darkred")) +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  ylab("Difference between 2023 production and 2022 hindcast") +
  facet_wrap(~ Year, ncol = 8) 
comp_diff

ggsave(comp_diff, filename = here("VAST_results", "2023_age_comp_diff.png"),
       width=200, height=130, units="mm", dpi=300)

# Plot hindcast and production runs together ----------------------------------
colnames(new_props)[1:15] <- 1:15
colnames(old_props)[1:15] <- 1:15
new_props_long <- melt(new_props, id.vars = c("Year", "Region"),
                       variable.name = "Age", value.name = "Proportion")
new_props_long$version <- "2023"
old_props_long <- melt(old_props, id.vars = c("Year", "Region"),
                       variable.name = "Age", value.name = "Proportion")
old_props_long$version <- "2022 hindcast"

all_props <- rbind.data.frame(new_props_long, old_props_long) %>%
  filter(Year != 2020 & Region == "EBS") %>%  # remove 2020
  ggplot(., aes(x = Age, y = Proportion, fill = version)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("gray", "darkslateblue")) +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  ylab("Proportion-at-age") +
  facet_wrap(~ Year, ncol = 8) 
all_props

ggsave(all_props, filename = here("VAST_results", "2023_age_comp_compare.png"),
       width=200, height=130, units="mm", dpi=300)


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

### ESP plots -----------------------------------------------------------------
# Center of gravity -----------------------------------------------------------
cog <- read.csv(here("VAST_results", "COG.csv")) 
cog$m[cog$m == 1] <- "Eastings"
cog$m[cog$m == 2] <- "Northings"

cog_plot <- ggplot(cog, aes(x = Year, y = COG_hat)) +
  geom_line(alpha = 0.4) +
  geom_pointrange(aes(ymin = (COG_hat - SE), ymax = (COG_hat + SE))) +
  ylab("Center of Gravity") +
  facet_wrap(~ m, ncol = 1, scales = "free_y")
cog_plot

ggsave(cog_plot, filename = here("VAST_results", "2023_pollock_COG.png"),
       width = 150, height = 180, unit = "mm", dpi = 300)

# Area occupied ---------------------------------------------------------------
options(scipen = 999)
area <- read.csv(here("VAST_results", "ln_effective_area.csv")) 
area$Region <- c(rep("Both", length(1982:this_year) - 1), 
                 rep("EBS", length(1982:this_year) - 1),
                 rep("NBS", length(1982:this_year) - 1))
colnames(area)[2] <- "error"
area$Estimate <- exp(area$Estimate)
area$error <- area$error

area_plot <- ggplot(area, aes(x = Year, y = Estimate, color = Region)) +
  geom_line(alpha = 0.4) +
  geom_pointrange(aes(ymin = (Estimate - (Estimate * error)), 
                      ymax = (Estimate + (Estimate * error))),
                  alpha = 0.8) +
  scale_y_continuous(labels = scales::comma, limits = c(0, NA)) +
  ylab("Effective area occupied (km^2)")
area_plot

ggsave(area_plot, filename = here("VAST_results", "2023_pollock_area.png"),
       width = 150, height = 100, unit = "mm", dpi = 300)
