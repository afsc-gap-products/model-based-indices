# Script for comparing two indices of abundance, either the hindcast and 
# production runs, or between models.
# By: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.03.14

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
workDir <- here("species_specific_code", "BS", "pollock")
this_year <- 2024

# save_dir <- paste0(this_year, " Production")
save_dir <- "tinyVAST"


# Compare Indices of Abundance ------------------------------------------------
# Read in indices & make sure columns for year = Time, Estimate, Error are named correctly
index1 <- read.csv(here(workDir, "results", "VAST Index", "Index.csv"))
colnames(index1)[6] <- "error"
index1$Estimate <- index1$Estimate / 1000000000
index1$error <- index1$error / 1000000000

index2 <- read.csv(here(workDir, "results", "2023 Production", "VAST Index", "Index.csv"))
colnames(index2)[6] <- "error"
index2$Estimate <- index2$Estimate / 1000000000
index2$error <- index2$error / 1000000000
 
# index3 <- read.csv(here(workDir, "results", "2023 Production", "VAST Index", "Index.csv"))
# colnames(index3)[6] <- "error"
# index3$Estimate <- index3$Estimate / 1000000000
# index3$error <- index3$error / 1000000000

# # When needed, sum across ages
# index2 <- index2 %>% group_by(Time, Stratum) %>%
#   summarize(Estimate = sum(Estimate),
#             error = mean(error))
# 
# # Get design-based index from density-dependent correction --------------------
# ddc_orig_biomass <- read.csv(here(workDir, "results", "design-based", "biomass_densdep_corrected_EBSonly_2024.csv"))
# ddc_index <- ddc_orig_biomass %>%
#   group_by(year) %>%
#   summarize(Estimate = sum(biomass_MT_ha) / 1000000)
# ddc_index$error <- 0
# ddc_index$Stratum <- "EBS"
# ddc_index <- ddc_index[, c(1, 4, 2, 3)]
# colnames(ddc_index)[1] <- "Time"
  
# Combine and plot indices ----------------------------------------------------
# Set names for old and new index
names_index <- c("DDC DB 2024", "MB 2024")

compare_index <- function(indices, names, ebs_only = FALSE) {
  df <- data.frame()
  for(i in 1:length(indices)) {
    index <- indices[[i]]
    index <- index[c("Time", "Stratum", "Estimate", "error")]
    # index$Estimate <- index$Estimate / 1000000000  # convert to million tons
    # index$error <- index$error / 1000000000  # convert to million tons
    index$version <- names[i]
    df <- rbind.data.frame(df, index)
  }
  
  df <- df %>% filter(Time !=2020)
  
  if(ebs_only == TRUE) {
    df <- df %>% filter(Time !=2020 & Stratum == "EBS") 
  }
  
  plot <- ggplot(df, 
                 aes(x = Time, y = Estimate, 
                     color = version, shape = version)) +
    geom_line(alpha = 0.3) +
    geom_pointrange(aes(ymin = Estimate - error, ymax = Estimate + error), 
                    alpha = 0.8) +
    ylim(0, NA) +
    xlab("Year") + ylab("Index of Abundance (Mt)") +
    scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.9) 
  
  if(ebs_only == FALSE) {
    plot <- plot + facet_wrap(~ Stratum, scales = "free_y", ncol = 1)
  }
  
  return(plot)
}

index_comp <- compare_index(indices = list(ddc_index, index1), 
                            names = c(names_index[1], names_index[2]),
                            ebs_only = TRUE)
index_comp

# Difference between two indices
index_difference <- function(new, old, names, save_results = FALSE) {
  # Only run for EBS estimate (and no 2020)
  new <- new %>% filter(Stratum == "EBS" & Time != 2020)
  old <- old %>% filter(Stratum == "EBS" & Time != 2020)
  # make sure new dataset is the same length as the old
  new <- new %>% filter(Time %in% min(old$Time):max(old$Time)) 
  df <- cbind.data.frame(Year = new$Time,
                         Stratum = new$Stratum,
                         Difference = ((new$Estimate - old$Estimate) / old$Estimate) * 100)

  if(save_results == TRUE) {  # save to drive, if you want. Check file paths.
    write.csv(df, here(results_dir, save_dir, "index_difference.csv"))
  }
  
  df <- df %>% mutate(sign = case_when(Difference >= 0 ~ "Positive",
                                       Difference < 0 ~ "Negative"))
  
  label <- paste0("Percent difference between ", names[1], " and ", names[2])
  plot <- ggplot(df, aes(x = Year, y = Difference, fill = sign)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = c("cornflowerblue", "darkred")) +
    ylab(label) + labs(color = "Difference")
  
  return(plot)
}

index_diff <- index_difference(new = index1, old = index2,
                               names = c(names_index[1], names_index[2]))
index_diff

# Compare Age Compositions ----------------------------------------------------
# Read in age comp model results (and remove rownames column)
# old_props <- read.csv(here(workDir, "results", "2023 Production", "Comps", "proportions.csv"))
new_props <- cbind(read.csv(here(workDir, "results", "Comps", "proportions.csv")), distribution = "VAST")[, -1]
# tiny_tweedie <- cbind(read.csv(here(workDir, "results", "tinyVAST", "tinyVAST_props.csv")), distribution = "Tweedie")
# tiny_dg <- cbind(read.csv(here(workDir, "results", "tinyVAST", "tinyVAST_props_dg.csv")), distribution = "Delta Gamma")
# tiny_mesh <- cbind(read.csv(here(workDir, "results", "tinyVAST", "tinyVAST_props_mesh.csv")), distribution = "Tweedie")
# tiny_mesh_dg <- cbind(read.csv(here(workDir, "results", "tinyVAST", "tinyVAST_props_mesh_dg.csv")), distribution = "Delta Gamma")
# tiny_logn <- cbind(read.csv(here(workDir, "results", "tinyVAST", "tinyVAST_props_logn.csv")), distribution = "Lognormal")
# tiny_bc <- cbind(read.csv(here(workDir, "results", "tinyVAST", "tinyVAST_props_bc.csv")), distribution = "Tweedie")
# tiny_dg_bc <- cbind(read.csv(here(workDir, "results", "tinyVAST", "tinyVAST_props_dg_bc.csv")), distribution = "Delta Gamma")
tiny_mesh_dglink <- cbind(read.csv(here(workDir, "results", "tinyVAST", "tinyVAST_props_mesh_dglink.csv")), distribution = "DG (poisson-link)")[, -17]
sdm_props <- dcast(cbind(read.csv(here(workDir, "results", "sdmTMB_age_prop.csv"))[, 2:4]), formula = year ~ Age)
# tiny_coarse <- cbind(read.csv(here(workDir, "results", "tinyVAST_props_coarse_2024.csv")), distribution = "DG (poisson-link)")
gapindex_comps_raw <- read.csv(here(workDir, "results", "Comps", "gapindex_comps_2024.csv"))

# Reshape sdmTMB comps 
colnames(sdm_props)[1] <- "Year"
sdm_props$distribution <- "sdmTMB"
sdm_props <- sdm_props[, c(2:16, 1, 17)]

# Reshape gapindex comps
gapindex_comps <- gapindex_comps_raw %>% 
  filter(Region == "Both")
gapindex_comps <- gapindex_comps[, -4] %>%
  dcast(formula = YEAR ~ AGE) 
colnames(gapindex_comps)[1] <- "Year"
gapindex_comps <- gapindex_comps[, c(2:16, 1)]
gapindex_comps$distribution <- "gapindex"

# Update old_props to match tinyVAST test output - only EBS
tiny_years <- c(1980:2019, 2021:this_year)
new_props <- new_props %>% filter(Year %in% tiny_years & Region == "Both")
new_props <- new_props[, -17]

# Set names for old and new comps
# names_comps <- c("original", "original", "original", "VAST mesh", "VAST mesh", "original", "bias correction", "bias correction")

## Combine age comp models into one plot --------------------------------------
compare_props <- function(props, names, last_year) {
  df <- data.frame()
  for(i in 1:length(props)) {
    prop <- props[[i]]
    colnames(prop)[1:15] <- 1:15
    prop <- melt(prop, id.vars = c("Year", "distribution"),
                 variable.name = "Age", value.name = "Proportion")
    prop$version <- names[i]
    df <- rbind.data.frame(df, prop)
  }
  
  barplot <- ggplot(df, aes(x = Age, y = Proportion, fill = version)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_x_discrete(breaks = c(1, 5, 10, 15)) +
    ylab("Proportion-at-age") +
    facet_wrap(~ Year, ncol = 6, dir = "v")
  
  boxplot <- ggplot(df, aes(x = Age, y = Proportion, color = version, fill = version)) +
    geom_boxplot(alpha = 0.5) +
    scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_x_discrete(breaks = c(1, 5, 10, 15)) +
    ylab("Proportion-at-age") 

  return(list(barplot = barplot, boxplot = boxplot))
}

# comp_all_plots <- compare_props(props = list(new_props, tiny_tweedie, tiny_dg, 
#                                              tiny_mesh, tiny_mesh_dg, tiny_logn, 
#                                              tiny_bc, tiny_dg_bc, sdm_props),
#                                 names = c("VAST", "1 - original", "1 - original", 
#                                           "3 - VAST mesh", "3 - VAST mesh", "1 - original", 
#                                           "2 - bias correction", "2 - bias correction", "sdmTMB"))

# summary_props_all <- comp_all_plots$boxplot + 
#   facet_wrap(~ distribution)
# summary_props_all

sum_props_sub <- compare_props(props = list(tiny_mesh_dglink, sdm_props, gapindex_comps),
                               names = c("tinyVAST", "sdmTMB", "gapindex"))
sum_props_sub$barplot
sum_props_sub$boxplot

# Plot difference between two models ------------------------------------------
comp_difference <- function(new, old, names, save_results = FALSE) {
  # new <- subset(new, new$Year < this_year)  # make sure new dataset is the same length as the old
  
  # Get difference between new and old props
  check_props <- round(new[,1:15] - old[,1:15], 4)
  check_props_tab <- cbind(check_props, new[,16:17])
  check_props_abs <- abs(check_props)
  check_props_abs_tab <-  cbind(check_props_abs, new[,16:17])
  
  if(save_results == TRUE) {  # save to drive, if you want. Check file paths.
    options(scipen=999)
    write.csv(check_props_tab, here(results_dir, save_dir, "bridge_props.csv"))
    write.csv(check_props_abs_tab, here(results_dir, save_dir, "bridge_props_abs.csv"))
  }
  
  colnames(check_props_tab)[1:15] <- 1:15
  props_plot <- melt(check_props_tab, id.vars = c("Year", "Region"), 
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
    facet_wrap(~ Year, ncol = 6, dir = "v") 
  
  return(plot)
}

comp_diff <- comp_difference(new = tiny_mesh, old = new_props,
                             names = c("Tweedie/old mesh", "VAST"),
                             save_results = FALSE)
comp_diff

# Percent difference
comp_percent_diff <- function(new, old, names, save_results = FALSE) {
  # new <- subset(new, new$Year < this_year)  # make sure new dataset is the same length as the old
  
  # Get difference between new and old props
  check_props <- round((((new[,1:15] - old[,1:15]) / old[,1:15]) *100), 4)
  check_props_tab <- cbind(check_props, new[,16:17])
  check_props_abs <- abs(check_props)
  check_props_abs_tab <-  cbind(check_props_abs, new[,16:17])
  
  if(save_results == TRUE) {  # save to drive, if you want. Check file paths.
    options(scipen=999)
    write.csv(check_props_tab, here(results_dir, save_dir, "bridge_props.csv"))
    write.csv(check_props_abs_tab, here(results_dir, save_dir, "bridge_props_abs.csv"))
  }
  
  colnames(check_props_tab)[1:15] <- 1:15
  props_plot <- melt(check_props_tab, id.vars = c("Year", "Region"), 
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
    facet_wrap(~ Year, ncol = 6) 
  
  return(plot)
}

per_diff <- comp_percent_diff(new = tiny_mesh, old = new_props,
                              names = c("Tweedie/old mesh", "VAST"),
                              save_results = FALSE)
per_diff

## Check trends in difference between models by age & year
comp_trends <- function(new, old, names) {
  # new <- subset(new, new$Year < this_year)  # make sure new dataset is the same length as the old
  
  # Get difference between new and old props
  check_props <- round(new[,1:15] - old[,1:15], 4)
  check_props_tab <- cbind(check_props, new[,16:17])
  check_props_abs <- round(abs(new[,1:15] - old[,1:15]), 4)
  check_props_abs_tab <-  cbind(check_props_abs, new[,16:17])
  
  colnames(check_props_tab)[1:15] <- 1:15
  props <- melt(check_props_tab, id.vars = c("Year", "Region"), 
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
  plot_year <- ggplot(props, aes(x = Year, y = Proportion, color = sign)) +
    # geom_violin() +
    geom_jitter(height = 0, width = 0.1, alpha = 0.6, show.legend = FALSE) +
    scale_color_manual(values = c("cornflowerblue", "darkred")) +
    ylab("")

  plot_both <- ggpubr::ggarrange(plot_age, plot_year)  # combine plots
  return(plot_both)
}

comp_trends <- comp_trends(new = tiny_mesh, old = new_props,
                           names = c("Tweedie/old mesh", "VAST"))
comp_trends

# tinyVAST plots save ---------------------------------------------------------
ggsave(summary_props_all, filename = here(workDir, "results", save_dir, "tiny_summary.png"),
       width=180, height=120, units="mm", dpi=300)
ggsave(sum_tweedie_all, filename = here(workDir, "results", save_dir, "tweedie_compare.png"),
       width=200, height=180, units="mm", dpi=300)
ggsave(sum_tweedie_sum, filename = here(workDir, "results", save_dir, "tweedie_summary.png"),
       width=200, height=120, units="mm", dpi=300)
ggsave(comp_diff, filename = here(workDir, "results", save_dir, "tweedie_mesh_diff.png"),
       width=200, height=200, units="mm", dpi=300)
ggsave(per_diff, filename = here(workDir, "results", save_dir, "tweedie_mesh_per_diff.png"),
       width=200, height=200, units="mm", dpi=300)
ggsave(comp_trends, filename = here(workDir, "results", save_dir, "tweedie_mesh_trends.png"),
       width=260, height=120, units="mm", dpi=300)

sum_props_sub$boxplot
ggsave(sum_props_sub$boxplot, filename = here(workDir, "results", save_dir, "tinyVAST_summary.png"),
       width=200, height=120, units="mm", dpi=300)

# Save plots ------------------------------------------------------------------
# ggsave(index_comp, filename = here(workDir, "results", save_dir, "index_comparison.png"),
#        width=170, height=120, units="mm", dpi=300)
# ggsave(index_diff, filename = here(workDir, "results", save_dir, "index_difference.png"),
#        width=170, height=120, units="mm", dpi=300)
# ggsave(all_props, filename = here(workDir, "results", save_dir, "age_comp_compare.png"),
#        width=200, height=180, units="mm", dpi=300)
# ggsave(summary_props, filename = here(workDir, "results", save_dir, "age_comp_summary.png"),
#        width=200, height=120, units="mm", dpi=300)
# ggsave(comp_diff, filename = here(workDir, "results", save_dir, "age_comp_diff.png"),
#        width=200, height=200, units="mm", dpi=300)
# ggsave(per_diff, filename = here(workDir, "results", save_dir, "age_comp_per_diff.png"),
#        width=200, height=200, units="mm", dpi=300)
# ggsave(comp_trends, filename = here(workDir, "results", save_dir, "age_comp_trends.png"),
#        width=260, height=120, units="mm", dpi=300)
