# Custom plots of VAST results for pollock
# By: Sophia Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2023.10.17

library(here)
library(tidyverse)
library(ggplot2)
library(viridis)
library(dplyr)
library(reshape2)
library(gapindex)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

this_year <- 2023

workDir <- here("species_specific_code","BS","pollock")  # define working directory

### VAST index ----------------------------------------------------------------
index <- read.csv(here(workDir, "results", "VAST Index", "Index.csv"))

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
  xlab("Year") + ylab("Index of Abundance (Mt)") +
  facet_wrap(~ Stratum, ncol = 1)
index_all_areas

ggsave(index_all_areas, filename = here(workDir, "plots", "index_all_areas.png"),
       width=130, height=160, units="mm", dpi=300)

# Plot just EBS
ebs <- index %>% filter(Stratum == "EBS" & Time != 2020)
colnames(ebs)[6] <- "error"
index_ebs <- ggplot(ebs, aes(x = Time, y = (Estimate / 1000000000))) +
  geom_line(alpha = 0.4) +
  geom_pointrange(aes(ymin = (Estimate / 1000000000) - (error / 1000000000),
                      ymax = (Estimate / 1000000000) + (error / 1000000000)), alpha = 0.8) +
  ylim(0, NA) +
  xlab("Year") + ylab("Index of Abundance (Mt)") 
index_ebs

ggsave(index_ebs, filename = here(workDir, "plots", "index_ebs.png"),
       width=130, height=90, units="mm", dpi=300)

# Compare to design-based index from gapindex ---------------------------------
sql_channel <- gapindex::get_connected()  # connect to Oracle

# Pull data
gapindex_data <- gapindex::get_data(
  year_set = c(1982:this_year),
  survey_set = "EBS",
  spp_codes = c(21740),   
  haul_type = 3,
  abundance_haul = "Y",
  pull_lengths = F,
  sql_channel = sql_channel)

# Fill zeros and calculate CPUE
cpue <- gapindex::calc_cpue(racebase_tables = gapindex_data)
# warning: no catch records found for species code 21741

# Calculate stratum-level biomass, abundance, mean CPUE, and variances
biomass_stratum <- gapindex::calc_biomass_stratum(
  racebase_tables = gapindex_data,
  cpue = cpue)
# Warning: EBS + NW output only includes years 1987-present

# Calculate aggregated biomass and abundance 
biomass_subareas <- gapindex::calc_biomass_subarea(
  racebase_tables = gapindex_data,
  biomass_strata = biomass_stratum)
# Warning: EBS + NW output only includes years 1987-present

# Plot design-based index
# Filter index to standard EBS & NW area 
db_index <- biomass_subareas %>%
  filter(AREA_ID == 99900)

ggplot(db_index, aes(x = YEAR, y = (BIOMASS_MT))) +
  geom_line(alpha = 0.4) +
  geom_point() +
  # geom_pointrange(aes(ymin = (BIOMASS_MT) - (BIOMASS_VAR),
  #                     ymax = (BIOMASS_MT) + (BIOMASS_VAR)), alpha = 0.8) +
  ylim(0, NA) +
  xlab("Year") + ylab("Index (Mt)") 

# Combine indices together and plot
db_mb <- rbind.data.frame(cbind.data.frame(Year = db_index$YEAR,
                                           Biomass = db_index$BIOMASS_MT / 1000000,  
                                           # Error = db_index$BIOMASS_VAR / 1000000,  # what's going on here???
                                           Error = 0, 
                                           Source = "gapindex"),
                          cbind.data.frame(Year = ebs$Time,
                                           Biomass = ebs$Estimate / 1000000000,  # to Mt
                                           Error = ebs$error / 1000000000,  # to Mt
                                           Source = "VAST")) %>%
  mutate(Source = factor(Source, levels = c("VAST", "gapindex"))) %>%
  ggplot(., aes(x = Year, y = Biomass, color = Source, shape = Source)) +
  geom_line(alpha = 0.4) +
  geom_pointrange(aes(ymin = (Biomass) - (Error),
                      ymax = (Biomass) + (Error))) +
  scale_color_manual(values = c("black", "darkslateblue")) +
  ylim(0, NA) +
  xlab("Year") + ylab("EBS Index (Mt)")
db_mb

ggsave(db_mb, filename = here(workDir, "plots", "2023_db_vast.png"),
       width=170, height=90, units="mm", dpi=300)

### Plot VAST age comps -------------------------------------------------------
proportions <- read.csv(here(workDir, "results", "Comps", "proportions.csv"))[, -1]
colnames(proportions)[1:15] <- 1:15
props <- melt(proportions, id.vars = c("Year", "Region"),
              variable.name = "Age", value.name = "Proportion")
props_ebs <- props %>% 
  filter(Region == "EBS") %>%
  arrange(Year, Age)

# Set up shifting colors for each year to track cohorts
colors <- rep(1:16, length(1982:this_year))
props_ebs$color <- colors[1:nrow(props_ebs)]

prop_plot <- ggplot(props_ebs, aes(x = Age, y = Proportion, fill = color)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis(option = "turbo") +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  scale_y_continuous(limits = c(0, 0.5), breaks = c(0, 0.2, 0.4)) +
  ylab("Proportion") + 
  guides(fill = "none") +
  facet_wrap(~ Year, ncol = 4, dir = "v") +
  theme(strip.text.x = element_blank()) +
  geom_text(x = 13, y = 0.45, aes(label = Year), color = "grey30", size = 2.8)
prop_plot

ggsave(prop_plot, filename = here(workDir, "plots", "2023_age_comp.png"),
       width=120, height=180, units="mm", dpi=300)

library(gganimate)
# Animated version, just for fun.
prop_gif <- ggplot(props_ebs, aes(x = Age, y = Proportion, fill = color)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis(option = "turbo") +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  scale_y_continuous(limits = c(0, 0.51), breaks = c(0, 0.2, 0.4)) +
  ylab("Proportion") + 
  guides(fill = "none") +
  labs(title = "Year: {next_state}") +
  transition_states(Year, transition_length = 1, state_length = 2, wrap = TRUE) +
  ease_aes("linear")
# animate(prop_gif)
anim_save(here(workDir, "plots", "2023_age_comp.gif"), prop_gif, fps = 2)


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

ggsave(cold_pool_plot, filename = here("output", "cold_pool_covariate.png"),
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
