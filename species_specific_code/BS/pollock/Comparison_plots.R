# Script for comparing two indices of abundance, either the hindcast and 
# production runs, or between models.
# By: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.03.14

library(here)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)

# Set ggplot theme
# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
theme_set(theme_sleek())

# Compare Indices of Abundance ------------------------------------------------
# Read in indices & make sure columns for year = Time, Estimate, Error are named correctly
index1 <- read.csv("Index.csv")
colnames(index1)[6] <- "error"
index2 <- read.csv("Index_2022.csv")
colnames(index2)[6] <- "error"

# Function combining & plotting any number of indices. 
compare_index <- function(indices, names){
  df <- data.frame()
  for(i in 1:length(indices)) {
    index <- indices[[i]]
    index$version <- names[i]
    df <- rbind.data.frame(df, index)
  }
  
  plot <- ggplot(df %>% filter(Time != 2020), 
                 aes(x = Time, y = (Estimate / 1000000000), 
                     color = version, shape = version)) +
    geom_line(alpha = 0.3) +
    geom_pointrange(aes(ymin = (Estimate / 1000000000) - (error / 1000000000),
                        ymax = (Estimate / 1000000000) + (error / 1000000000)), alpha = 0.8) +
    ylim(0, NA) +
    xlab("Year") + ylab("Index (Mt)") +
    scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    facet_wrap(~ Stratum, scales = "free_y", ncol = 1)
  
  return(plot)
}

compare_index(indices = list(index1, index2), 
              names = c("2023 production", "2022 production"))


# Compare Age Compositions ----------------------------------------------------
# Read in age comp model results (and remove rownames column)
new_props <- read.csv("proportions.csv")[, -1]
old_props <- read.csv("proportions_2022.csv")[, -1]

## Combine age comp models into one plot --------------------------------------
compare_props <- function(props, names, last_year) {
  df <- data.frame()
  for(i in 1:length(props)) {
    prop <- props[[i]]
    colnames(prop)[1:15] <- 1:15
    prop <- melt(prop, id.vars = c("Year", "Region"),
                 variable.name = "Age", value.name = "Proportion")
    prop$version <- names[i]
    df <- rbind.data.frame(df, prop)
  }
  
  plot <- ggplot(df %>% filter(Year != 2020 & Region == "EBS"),
                 aes(x = Age, y = Proportion, fill = version)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_x_discrete(breaks = c(1, 5, 10, 15)) +
    ylab("Proportion-at-age") +
    facet_wrap(~ Year, ncol = 8) 
  return(plot)
}

all_props <- compare_props(props = list(new_props, old_props),
                           names = c("2023 production", "2023 hindcast"))
all_props

# Plot difference between two models ------------------------------------------
new_props <- subset(new_props, new_props$Year < 2023)

comp_difference <- function(new, old, names, save_results = FALSE) {
  # Get difference between new and old props
  check_props <- round(new[,1:15] - old[,1:15], 4)
  check_props_tab <- cbind(check_props, new[,16:17])
  check_props_abs <- round(abs(new[,1:15] - old[,1:15]), 4)
  check_props_abs_tab <-  cbind(check_props_abs, new[,16:17])
  
  if(save_results == TRUE) {  # save to drive, if you want. Check file paths.
    options(scipen=999)
    write.csv(check_props_tab, here("VAST_results", "bridge_props_2023.csv"))
    write.csv(check_props_abs_tab, here("VAST_results", "bridge_props_abs_2023.csv"))
  }
  
  colnames(check_props_tab)[1:15] <- 1:15
  props_plot <- melt(check_props_tab, id.vars = c("Year", "Region"), 
                     variable.name = "Age", value.name = "Proportion") %>%
    # add column for coloring the bars in the plot based on positive/negative
    mutate(sign = case_when(Proportion >= 0 ~ "positive",
                            Proportion < 0 ~ "negative"))
  
  # Plot both regions together and without 2020
  label = paste0("Difference between ", names[1], " and ", names[2])
  plot <- ggplot(props_plot %>% filter(Region == "EBS" & Year != 2020), 
                      aes(x = Age, y = Proportion, fill = sign)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = c("cornflowerblue", "darkred")) +
    scale_x_discrete(breaks = c(1, 5, 10, 15)) +
    ylab(label) +
    facet_wrap(~ Year, ncol = 8) 
  
  return(plot)
}

comp_diff <- comp_difference(new = new_props, old = old_props,
                             names = c("2023 production", "2023 hindcast"))
comp_diff

# Save plots ------------------------------------------------------------------
ggsave(all_props, filename = here("VAST_results", "2023_age_comp_compare.png"),
       width=200, height=130, units="mm", dpi=300)
ggsave(comp_diff, filename = here("VAST_results", "2023_age_comp_diff.png"),
       width=200, height=130, units="mm", dpi=300)

