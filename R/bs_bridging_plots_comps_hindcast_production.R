# Script for comparing comps from hindcast to production
# By: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.03.14
# modified by Lewis Barnett on 2025.09.25

library(here)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)

# Set ggplot theme
if (!requireNamespace("ggsidekick", quietly = TRUE)) {
  devtools::install_github("seananderson/ggsidekick")
}
library(ggsidekick)
theme_set(theme_sleek())

# Set up ----------------------------------------------------------------------
phase <- c("hindcast", "production")[2] # specify analysis phase

sp <- 1 # specify species from species vector
species <- c("yellowfin_sole", "pollock", "pacific_cod")[sp]
age_classes <- c(20, 15, 13)[sp]
start_age <- c(1 , 1, 0)[sp]
max_age <- c(20, 15, 12)[sp]

# Set year
this_year <- as.numeric(format(Sys.Date(), "%Y"))
if(phase == "hindcast") {this_year <- this_year - 1}  

# Set working directory specific to species & phase
workDir <- here("species_specific_code", "BS", species, phase)

# Compare Age Compositions ----------------------------------------------------
new <- cbind(read.csv(here(workDir, "results_age/tinyVAST_props.csv")), version = "production")
old <- cbind(read.csv(here("species_specific_code", "BS", species, "hindcast", 
                      "results_age", "tinyVAST_props_2024.csv")), version = "hindcast")

new <- dplyr::bind_rows(new) |>
        filter(year < this_year) |> # drop current year from production output
        select(-region)

old <- bind_cols(year = old$year, old[, 1:age_classes], version = old$version) # unify old format

compare_props <- function(props, names) {
  df <- data.frame()
  for(i in 1:length(props)) {
    prop <- props[[i]]
    prop <- prop[, names(prop) != "version"]
    colnames(prop)[2:(age_classes + 1)] <- start_age:max_age
    prop <- melt(prop, id.vars = "year",
                 variable.name = "age", value.name = "proportion")
    prop$version <- names[i]
    df <- rbind.data.frame(df, prop)
  }
  
  barplot <- ggplot(df, aes(x = age, y = proportion, fill = version)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_x_discrete(breaks = pretty(start_age:max_age)) +
    ylab("proportion-at-age") +
    facet_wrap(~ year, ncol = 6, dir = "v")
  
  boxplot <- ggplot(df, aes(x = age, y = proportion, color = version, fill = version)) +
    geom_boxplot(alpha = 0.5) +
    scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_x_discrete(breaks = pretty(start_age:max_age)) +
    ylab("proportion-at-age") 
  
  return(list(barplot = barplot, boxplot = boxplot))
}


sum_props <- compare_props(props = list(new, old),
                           names = c("production", "hindcast"))
sum_props$barplot

ggsave(sum_props$barplot, filename = here(workDir, "results_age", "bridge_comp_by_year.png"),
       width=200, height=180, units="mm", dpi=300)