# Script for comparing comps
# By: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.03.14
# modified by Lewis Barnett on 2025.03.05

library(here)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)

# Set ggplot theme
# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
theme_set(theme_sleek())

# Directories
workDir <- here("species_specific_code", "BS", "pacific_cod")
this_year <- 2024

# save_dir <- paste0(this_year, " Production")
save_dir <- paste0("hindcast/results_age")

# Compare Age Compositions ----------------------------------------------------
new_props <- cbind(read.csv("species_specific_code/BS/pacific_cod/hindcast/results_age/tinyVAST_props_2024.csv"), distribution = "tinyVAST")
old_props <- read.csv(here(workDir, "production", "results_age", "proportions", "clean_proportions.csv"))[c(1:26,28:31),-1]

old_props$year <- new_props$year
old_props$distribution <- "VAST"
names(old_props) <- names(new_props)

## Combine age comp models into one plot --------------------------------------
compare_props <- function(props, names, last_year) {
  df <- data.frame()
  for(i in 1:length(props)) {
    prop <- props[[i]]
    colnames(prop)[1:13] <- 0:12
    prop <- melt(prop, id.vars = c("year", "distribution"),
                 variable.name = "Age", value.name = "Proportion")
    prop$version <- names[i]
    df <- rbind.data.frame(df, prop)
  }
  
  barplot <- ggplot(df, aes(x = Age, y = Proportion, fill = version)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_x_discrete(breaks = c(1, 5, 10, 15)) +
    ylab("Proportion-at-age") +
    facet_wrap(~ year, ncol = 6, dir = "v")
  
  boxplot <- ggplot(df, aes(x = Age, y = Proportion, color = version, fill = version)) +
    geom_boxplot(alpha = 0.5) +
    scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_x_discrete(breaks = c(1, 5, 10, 15)) +
    ylab("Proportion-at-age") 

  return(list(barplot = barplot, boxplot = boxplot))
}

sum_props_sub <- compare_props(props = list(new_props, old_props),
                               names = c("tinyVAST", "VAST"))
sum_props_sub$barplot
sum_props_sub$boxplot

# Plot difference between two models ------------------------------------------
comp_difference <- function(new, old, names, save_results = FALSE) {
  # new <- subset(new, new$year < this_year)  # make sure new dataset is the same length as the old
  
  # Get difference between new and old props
  check_props <- round(new[,1:13] - old[,1:13], 4)
  check_props_tab <- cbind(check_props, new[,14:15])
  check_props_abs <- abs(check_props)
  check_props_abs_tab <-  cbind(check_props_abs, new[,14:15])
  
  if(save_results == TRUE) {  # save to drive, if you want. Check file paths.
    options(scipen=999)
    write.csv(check_props_tab, here(results_dir, save_dir, "bridge_props.csv"))
    write.csv(check_props_abs_tab, here(results_dir, save_dir, "bridge_props_abs.csv"))
  }
  
  colnames(check_props_tab)[1:13] <- 0:12
  props_plot <- melt(check_props_tab, id.vars = "year", 
                     variable.name = "Age", value.name = "Proportion") %>%
    # add column for coloring the bars in the plot based on positive/negative
    mutate(sign = case_when(Proportion >= 0 ~ "positive",
                            Proportion < 0 ~ "negative"))
  
  # Plot both regions together and without 2020
  label <- paste0("Difference between ", names[1], " and ", names[2])
  plot <- ggplot(props_plot, aes(x = Age, y = Proportion, fill = sign)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = c("cornflowerblue", "darkred")) +
    scale_x_discrete(breaks = c(1, 5, 10, 15)) +
    ylab(label) +
    facet_wrap(~ year, ncol = 6, dir = "v") 
  
  return(plot)
}

comp_diff <- comp_difference(new = new_props, old = old_props,
                             names = c("tinyVAST", "VAST"),
                             save_results = FALSE)
comp_diff

# Percent difference
comp_percent_diff <- function(new, old, names, save_results = FALSE) {
  # new <- subset(new, new$year < this_year)  # make sure new dataset is the same length as the old
  
  # Get difference between new and old props
  check_props <- round((((new[,1:13] - old[,1:13]) / old[,1:13]) *100), 4)
  check_props_tab <- cbind(check_props, new[,14:15])
  check_props_abs <- abs(check_props)
  check_props_abs_tab <-  cbind(check_props_abs, new[,14:15])
  
  if(save_results == TRUE) {  # save to drive, if you want. Check file paths.
    options(scipen=999)
    write.csv(check_props_tab, here(results_dir, save_dir, "bridge_props.csv"))
    write.csv(check_props_abs_tab, here(results_dir, save_dir, "bridge_props_abs.csv"))
  }
  
  colnames(check_props_tab)[1:13] <- 0:12
  props_plot <- melt(check_props_tab, id.vars = "year", 
                     variable.name = "Age", value.name = "Proportion") %>%
    # add column for coloring the bars in the plot based on positive/negative
    mutate(sign = case_when(Proportion >= 0 ~ "positive",
                            Proportion < 0 ~ "negative"))
  
  # Plot both regions together and without 2020
  label <- paste0("Percent difference between ", names[1], " and ", names[2])
  plot <- ggplot(props_plot, aes(x = Age, y = Proportion, fill = sign)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = c("cornflowerblue", "darkred")) +
    scale_x_discrete(breaks = c(1, 5, 10, 15)) +
    ylab(label) +
    facet_wrap(~ year, ncol = 6) 
  
  return(plot)
}

per_diff <- comp_percent_diff(new = new_props, old = old_props,
                              names = c("tinyVAST", "VAST"),
                              save_results = FALSE)
per_diff

## Check trends in difference between models by age & year
comp_trends <- function(new, old, names) {
  # new <- subset(new, new$year < this_year)  # make sure new dataset is the same length as the old
  
  # Get difference between new and old props
  check_props <- round(new[,1:13] - old[,1:13], 4)
  check_props_tab <- cbind(check_props, new[,14:15])
  check_props_abs <- round(abs(new[,1:13] - old[,1:13]), 4)
  check_props_abs_tab <-  cbind(check_props_abs, new[,14:15])
  
  colnames(check_props_tab)[1:13] <- 0:12
  props <- melt(check_props_tab, id.vars = "year", 
                     variable.name = "Age", value.name = "Proportion") %>%
    # add column for coloring the bars in the plot based on positive/negative
    mutate(sign = case_when(Proportion >= 0 ~ "positive",
                            Proportion < 0 ~ "negative"))
  
  # Plot difference over ages
  label <- paste0("Difference between ", names[1], " and ", names[2])
  plot_age <- ggplot(props, aes(x = Age, y = Proportion, color = sign)) +
    # geom_violin() +
    geom_jitter(height = 0, width = 0.1, alpha = 0.6, show.legend = FALSE) +
    scale_color_manual(values = c("cornflowerblue", "darkred")) +
    scale_x_discrete(breaks = c(1, 5, 10, 15)) +
    ylab(label) 
  
  # Plot difference over years
  plot_year <- ggplot(props, aes(x = year, y = Proportion, color = sign)) +
    # geom_violin() +
    geom_jitter(height = 0, width = 0.1, alpha = 0.6, show.legend = FALSE) +
    scale_color_manual(values = c("cornflowerblue", "darkred")) +
    ylab("")

  plot_both <- ggpubr::ggarrange(plot_age, plot_year)  # combine plots
  return(plot_both)
}

comp_trends <- comp_trends(new = new_props, old = old_props,
                           names = c("tinyVAST", "VAST"))
comp_trends

# tinyVAST plots save ---------------------------------------------------------
# ggsave(comp_diff, filename = here(workDir, "results", save_dir, "comp_diff.png"),
#        width=200, height=200, units="mm", dpi=300)
# ggsave(per_diff, filename = here(workDir, "results", save_dir, "comp_per_diff.png"),
#        width=200, height=200, units="mm", dpi=300)
# ggsave(comp_trends, filename = here(workDir, "results", save_dir, "comp_trends.png"),
#        width=260, height=120, units="mm", dpi=300)
# ggsave(sum_props_sub$boxplot, filename = here(workDir, "results", save_dir, "tinyVAST_summary.png"),
#        width=200, height=120, units="mm", dpi=300)
# ggsave(sum_props_sub$barplot, filename = here(workDir, "results", save_dir, "tinyVAST_by_year.png"),
#        width=200, height=120, units="mm", dpi=300)