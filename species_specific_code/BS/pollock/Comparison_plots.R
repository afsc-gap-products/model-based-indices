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

save_dir <- paste0(this_year, " production")

# Compare Indices of Abundance ------------------------------------------------
# Read in indices & make sure columns for year = Time, Estimate, Error are named correctly
index1 <- read.csv(here(workDir, "results", "VAST Index", "Index.csv"))
colnames(index1)[6] <- "error"
index1$Estimate <- index1$Estimate / 1000000000
index1$error <- index1$error / 1000000000

# index2 <- read.csv(here(workDir, "results", "2024 Hindcast", "VAST Index", "Index.csv"))
# colnames(index2)[6] <- "error"
# index2$Estimate <- index2$Estimate / 1000000000
# index2$error <- index2$error / 1000000000
# 
# index3 <- read.csv(here(workDir, "results", "2023 Production", "VAST Index", "Index.csv"))
# colnames(index3)[6] <- "error"
# index3$Estimate <- index3$Estimate / 1000000000
# index3$error <- index3$error / 1000000000

# Get design-based index from GAP_PRODUCTS - first connect to Oracle
if (file.exists("Z:/Projects/ConnectToOracle.R")) {
  source("Z:/Projects/ConnectToOracle.R")
} else {
  # For those without a ConnectToOracle file
  channel <- odbcConnect(dsn = "AFSC", 
                         uid = rstudioapi::showPrompt(title = "Username", 
                                                      message = "Oracle Username", 
                                                      default = ""), 
                         pwd = rstudioapi::askForPassword("Enter Password"),
                         believeNRows = FALSE)
}

# check to see if connection has been established
odbcGetInfo(channel)

# Pull index 
wp_db_index <- RODBC::sqlQuery(channel = channel,
                               query = "SELECT 
                                    YEAR, 
                                    CASE
                                     WHEN AREA_ID = 99900 THEN 'EBS'
                                     WHEN AREA_ID = 99902 THEN 'NBS'
                                    END AS REGION, 
                                    BIOMASS_MT,
                                    SQRT(BIOMASS_VAR)/BIOMASS_MT AS BIOMASS_CV
                                    FROM GAP_PRODUCTS.BIOMASS
                                    WHERE SPECIES_CODE = 21740
                                    AND AREA_ID in (99900, 99902)")

colnames(wp_db_index) <- c("Time", "Stratum", "Estimate", "error")
wp_db_index$Estimate <- wp_db_index$Estimate / 1000000
wp_db_index$error <- wp_db_index$error / 1000000


# # When needed, sum across ages
# index2 <- index2 %>% group_by(Time, Stratum) %>%
#   summarize(Estimate = sum(Estimate),
#             error = mean(error)) 

# Set names for old and new index
names_index <- list(old = "2024 Design-based", new = "2024 Model-based")

# Combine & plot any number of indices. 
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

index_comp <- compare_index(indices = list(wp_db_index, index1), 
                            names = c(names_index$old, names_index$new),
                            ebs_only = TRUE)
index_comp
ggsave(index_comp, filename = here(workDir, "plots", "index_dbmb_2024.png"),
       width=130, height=90, units="mm", dpi=300)

# Difference between two indices
index_difference <- function(new, old, names, save_results = FALSE) {
  new <- new <- subset(new, new$Time < this_year)  # make sure new dataset is the same length as the old
  df <- cbind.data.frame(Year = new$Time,
                         Stratum = new$Stratum,
                         Difference = ((new$Estimate - old$Estimate) / old$Estimate) * 100)

  if(save_results == TRUE) {  # save to drive, if you want. Check file paths.
    write.csv(df, here(results_dir, save_dir, "index_difference.csv"))
  }
  
  df <- df %>% 
    filter(Stratum == "EBS") %>% # only Eastern Bering Sea for simplicity
    # add column for coloring the bars in the plot based on positive/negative
    mutate(sign = case_when(Difference >= 0 ~ "Positive",
                            Difference < 0 ~ "Negative")) %>%
    filter(Year != 2020)
  
  label <- paste0("Percent difference between ", names[1], " and ", names[2])
  plot <- ggplot(df, aes(x = Year, y = Difference, fill = sign)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = c("cornflowerblue", "darkred")) +
    ylab(label) + labs(color = "Difference")
  
  return(plot)
}

index_diff <- index_difference(new = index1, old = wp_db_index,
                               names = c(names_index$new, names_index$old))
index_diff
ggsave(index_diff, filename = here(workDir, "plots", "index_dbmb_diff_2024.png"),
       width=130, height=110, units="mm", dpi=300)

# Compare Age Compositions ----------------------------------------------------
# Read in age comp model results (and remove rownames column)
old_props <- read.csv(here(workDir, "results", "2024 Hindcast", "Comps", "proportions.csv"))[, -1]
new_props <- read.csv(here(workDir, "results", "Comps", "proportions.csv"))
# new_props_dg <- read.csv(here(workDir, "results", "tinyVAST", "tinyVAST_props_dg.csv"))

# # Update old_props to match tinyVAST test output
# tiny_years <- c(1980:2019, 2021:2023)
# old_props <- old_props %>% filter(Year %in% tiny_years & Region == "EBS")

# Set names for old and new comps
names_comps <- list(old = "2024 Hindcast", new = "2024 Production")

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
  
  plot <- ggplot(df, aes(x = Age, y = Proportion, fill = version)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_x_discrete(breaks = c(1, 5, 10, 15)) +
    ylab("Proportion-at-age") +
    facet_wrap(~ Year, ncol = 6) 
  return(plot)
}

all_props <- compare_props(props = list(new_props, old_props),
                           names = c(names_comps$new, names_comps$old))
all_props

# Plot difference between two models ------------------------------------------
comp_difference <- function(new, old, names, save_results = FALSE) {
  new <- subset(new, new$Year < this_year)  # make sure new dataset is the same length as the old
  
  # Get difference between new and old props
  check_props <- round(new[,1:15] - old[,1:15], 4)
  check_props_tab <- cbind(check_props, new[,16:17])
  check_props_abs <- round(abs(new[,1:15] - old[,1:15]), 4)
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
    facet_wrap(~ Year, ncol = 6) 
  
  return(plot)
}

comp_diff <- comp_difference(new = new_props, old = old_props,
                             names = c(names$new, names$old),
                             save_results = FALSE)
comp_diff

## Check trends in difference between models by age & year
comp_trends <- function(new, old, names) {
  new <- subset(new, new$Year < this_year)  # make sure new dataset is the same length as the old
  
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

comp_trends <- comp_trends(new = new_props, old = old_props,
                           names = c(names$new, names$old))
comp_trends

 # Save plots ------------------------------------------------------------------
ggsave(index_comp, filename = here(workDir, "results", save_dir, "index_comparison.png"),
       width=130, height=160, units="mm", dpi=300)
ggsave(index_diff, filename = here(workDir, "results", save_dir, "index_difference.png"),
       width=130, height=180, units="mm", dpi=300)
ggsave(all_props, filename = here(workDir, "results", save_dir, "age_comp_compare.png"),
       width=200, height=180, units="mm", dpi=300)
ggsave(comp_diff, filename = here(workDir, "results", save_dir, "age_comp_diff.png"),
       width=200, height=200, units="mm", dpi=300)
ggsave(comp_trends, filename = here(workDir, "results", save_dir, "age_comp_trends.png"),
       width=260, height=120, units="mm", dpi=300)
