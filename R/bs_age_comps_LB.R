# Run all Bering Sea age compositions in tinyVAST
# Author: Sophia Wassermann
# Date: 26-3-2025

library(tinyVAST)
library(fmesher)
library(sf)
library(here)
library(dplyr)
library(ggplot2)

# Set ggplot theme
if (!requireNamespace("ggsidekick", quietly = TRUE)) {
  devtools::install_github("seananderson/ggsidekick")
}
library(ggsidekick)
theme_set(theme_sleek())

# Set up ----------------------------------------------------------------------
phase <- c("hindcast", "production")[1] # specify analysis phase

sp <- 1 # specify species from species vector
species <- c("yellowfin_sole", "pollock", "pacific_cod")[sp]

# Set year
this_year <- as.numeric(format(Sys.Date(), "%Y"))
if(phase == "hindcast") {this_year <- this_year - 1}  

# Set region (one of "EBS", "NBS", or "both")
region <- "both"

# Set working directory specific to species & phase
workDir <- here("species_specific_code", "BS", species, phase)

# Read in and format data -----------------------------------------------------
if(species == "pollock"){
  dat <- read.csv(here(workDir, "data", paste0("VAST_ddc_alk_", this_year, ".csv")))  
  dat <- transmute(dat,
                   cpue = CPUE_num * 100, # converts cpue from kg/ha to kg/km^2
                   year = as.integer(Year),
                   lat = Lat,
                   lon = Lon,
                   age = Age
  )
}

if(species %in% c("yellowfin_sole", "pacific_cod")){
  dat <- readRDS(here(workDir, "data", "data_geostat_agecomps.RDS"))  
  dat <- rename(dat, cpue = cpue_n_km2)
}

dat <- dat[!is.na(dat$cpue),]

# Add Year-Age interaction
ages <- unique(dat$age)
dat$age_f <- factor(paste0("age_", dat$age))
dat$year_age <- interaction(dat$year, dat$age_f)

# Project data to UTM
dat <- st_as_sf(dat, coords = c('lon','lat'), crs = st_crs(4326))
dat <- st_transform(dat, crs = st_crs("+proj=utm +zone=2 +units=km"))
# Add UTM coordinates as columns X & Y
dat <- cbind(st_drop_geometry(dat), st_coordinates(dat))

# For yellowfin sole and Pacific cod, drop data for levels with all zeros 
if(species %in% c("yellowfin_sole", "pacific_cod")) {
  N_z <- tapply(dat$cpue, 
                INDEX = dat$year_age, 
                FUN = \(x) sum(x > 0))
  year_age_to_drop <- names(which(N_z == 0))
  dat <- subset(dat, !(year_age %in% year_age_to_drop))
  dat <- droplevels(dat)
  
  # Check results
  tapply(dat$cpue,
         INDEX = list(dat$age, dat$year), 
         FUN = \(x) sum(x > 0))
}

if(species == "pollock") {
  year_age_to_drop <- NULL
}

# Check for ages with < 100 positive hauls ever (VAST would pool for hyperparams)
N_a = tapply(dat$cpue,
             INDEX = dat$age,
             FUN=\(x) sum(x>0) )
which(N_a < 100)

# Inputs to tinyVAST ----------------------------------------------------------
sem <- "\n  "
for(i in min(ages):max(ages)) {
  sem <- paste0(sem, "age_",i, " <-> age_", i, ", sd", i, "\n  ")
}

# Constant AR1 spatio-temporal term across ages & different variances for each age
dsem <- "\n  "
for(i in min(ages):max(ages)) {
  dsem <- paste0(dsem,  "age_", i, " -> age_", i, ", 1, lag1\n")
}

# # TinyVAST mesh
# mesh <- fm_mesh_2d(loc = dat[,c("X","Y")],
#                    cutoff = 50)

# Format mesh for tinyVAST (using a function from sdmTMB; same method as sdmTMB index bridging)
old_mesh <- sdmTMB::make_mesh(dat, 
                              xy_cols = c("X", "Y"), 
                              mesh = fmesher::fm_as_fm(readRDS(file = here("meshes/bs_vast_mesh_50_knots.RDS"))),
                              fmesher_func = fm_mesh_2d()) 

# Fit model -------------------------------------------------------------------
fit <- tinyVAST(
  formula = cpue ~ 0 + year_age,
  data = dat,
  space_term = sem,
  spacetime_term = dsem,
  family = setNames(
    lapply(ages, function(x) delta_gamma(type = "poisson-link")), 
    paste0("age_", ages)
    ),
  space_columns = c("X", "Y"),
  spatial_domain = old_mesh$mesh,
  time_column = "year",
  variable_column = "age_f",
  distribution_column = "age_f",
  delta_options = list(
    formula = ~ 0 + year_age,
    space_term = sem,
    spacetime_term = dsem
    ),
  control = tinyVASTcontrol(
    getsd = TRUE,
    silent = FALSE
    # , profile = c("alpha_j", "alpha2_j"), # for experimentation
    # , newton_loops = 1, # add newton loop(s) as needed to improve convergence
    # , tmb_par = fit$parameter_estimates # restart at prior best parameters
  )
)
fit$run_time
sanity(fit)

# Save fit object (create directory for results first, if it doesn't exist)
if (!dir.exists(here(workDir, "results_age"))) {
  dir.create(here(workDir, "results_age"))
}

saveRDS(fit, here(workDir, "results_age", "tinyVAST_fit.RDS"))

# Age composition expansion ---------------------------------------------------
# Load fit object if needed
if(!exists("fit")) {
  fit <- readRDS(here(workDir, "updated tinyVAST", "tinyVAST_fit.RDS"))
}

start <- Sys.time()
get_abundance <- function(region) {
  # Read in coarsened extrapolation grid
  if(region == "EBS") {grid <- read.csv(here("extrapolation_grids", "ebs_coarse_grid_5nm.csv"))}
  if(region == "NBS") {grid <- read.csv(here("extrapolation_grids", "nbs_coarse_grid_5nm.csv"))}
  if(region == "both") {grid <-  read.csv(here("extrapolation_grids", "bering_coarse_grid_5nm.csv"))}
  
  N_jz <- expand.grid(age_f = fit$internal$variables, year = sort(unique(dat$year)))
  N_jz$year_age <- interaction(N_jz$year, N_jz$age)
  N_jz <- cbind(N_jz, "abundance" = NA, "SE" = NA)
  
  areas <- newdata <- NULL
  
  for(j in which(!(N_jz$year_age %in% year_age_to_drop))) {
    # Make inputs, including 'block' column
      tmp <- data.frame(
        X = grid$X,
        Y = grid$Y,
        year = N_jz[j, "year"], 
        age_f = N_jz[j, "age_f"]
      )
      tmp$year_age <- paste(tmp$year, tmp$age_f, sep=".")
      newdata <- rbind(newdata, cbind(tmp, block = j))
      areas <- c(areas, grid$area_km2)
  }
  
  # Area-expansion (split into chunks to avoid memory limitation)
  x <- 1:nrow(newdata)
  chunk_size <- max(x) / 4
  chunk <- split(x, ceiling(seq_along(x) / chunk_size))
  
    gc()
    index1 <- integrate_output(
      fit,
      area = areas[chunk[[1]]],
      block = newdata$block[chunk[[1]]],
      newdata = newdata[chunk[[1]],],
      apply.epsilon = TRUE,
      bias.correct = FALSE,
      intern = TRUE,
      getsd = FALSE
    )
    gc()
    index2 <- integrate_output(
      fit,
      area = areas[chunk[[2]]],
      block = newdata$block[chunk[[2]]],
      newdata = newdata[chunk[[2]],],
      apply.epsilon = TRUE,
      bias.correct = FALSE,
      intern = TRUE,
      getsd = FALSE
    )
    gc()
    index3 <- integrate_output(
      fit,
      area = areas[chunk[[3]]],
      block = newdata$block[chunk[[3]]],
      newdata = newdata[chunk[[3]],],
      apply.epsilon = TRUE,
      bias.correct = FALSE,
      intern = TRUE,
      getsd = FALSE
    )
    gc()
    index4 <- integrate_output(
      fit,
      area = areas[chunk[[4]]],
      block = newdata$block[chunk[[4]]],
      newdata = newdata[chunk[[4]],],
      apply.epsilon = TRUE,
      bias.correct = FALSE,
      intern = TRUE,
      getsd = FALSE
    )
  index <- rbind(index1, index2, index3, index4)
    
  N_jz[, "abundance"] <- index[3] / 1e9 # scale to billions of individuals
  N_jz[is.na(N_jz)] <- 0 # replace NAs for combinations with 0 encounters
  
  N_ct <- array(N_jz$abundance, 
                dim = c(length(fit$internal$variables), length(unique(dat$year))),
                dimnames = list(fit$internal$variables, sort(unique(dat$year))))
  return(N_ct)
}

abundance <- get_abundance(region = region)
end <- Sys.time()
expand_time <- difftime(end, start, units = "hours")
expand_time

# Save abundance
write.csv(abundance, here(workDir, "results_age", "tinyVAST_abundance.csv"), row.names = FALSE)

# Run expansion for other regions
# ebs <- get_abundance(region = "EBS")
# gc()
# nbs <- get_abundance(region = "NBS")
# gc()
# both <- get_abundance(region = "both")

# Calculate proportions & plot ------------------------------------------------
calc_props <- function(df, area = region) {
  prop <- df / outer(rep(1, nrow(df)), colSums(df))
  prop <- tibble::rownames_to_column(data.frame(t(prop)), "year")
  prop$region <- area
  prop$year <- as.integer(prop$year)
  return(prop)
}

props <- calc_props(abundance)

# props <- rbind.data.frame(calc_props(ebs, "EBS"), 
#                           calc_props(nbs, "NBS"), 
#                           calc_props(both, "both"))

# Save proportions
write.csv(props, here(workDir, "results_age", "tinyVAST_props.csv"), row.names = FALSE)

# Plot proportions with colors to track cohort strength
plot_proportions <- function(area = region) {
  props <- props %>% filter(region == area)  # in case there are proportions calculated for multiple areas
  colors <- rep(1:(length(ages) + 1), length(min(props$year):this_year))  # color ID for the plot
  colnames(props) <- c("year", min(ages):max(ages), "region")
  plot <- reshape2::melt(props, 
                              id.vars = c("year", "region"),
                              variable.name = "age", 
                              value.name = "proportion") %>%
    arrange(year, age) %>%  # reformat to make colors work
    mutate(color = colors[1:nrow(.)]) %>%  # colors tracking cohorts
    ggplot(., aes(x = age, y = proportion, fill = color)) +
    geom_bar(stat = "identity", position = "dodge") +
    viridis::scale_fill_viridis(option = "turbo") +  # colorblind-friendly rainbow palette
    scale_x_discrete(breaks = pretty(ages)) +  # better breakpoints
    # scale_y_continuous(limits = c(0, 0.5), breaks = c(0, 0.2, 0.4)) +
    # ylab(paste0("proportion (", area, ")")) +   # define the region in the y-axis legend
    guides(fill = "none") +  # no legend 
    facet_wrap(~ year, ncol = 4, dir = "v") +  # years fill column-wise for cohort tracking
    theme(strip.text.x = element_blank()) +  # remove year label from top of plot & move to inside boxes
    geom_text(aes(x = max(ages) * 0.89, y = max(proportion) - 0.05, label = year), 
              color = "grey30", size = 2.8)
  
  return(plot)
}

props_plot <- plot_proportions()
props_plot

ggsave(props_plot, filename = here(workDir, "results_age", paste0("Age_comp_", region, ".png")),
       width=120, height=180, units="mm", dpi=300)

# Model diagnostics -----------------------------------------------------------
# working off of this vignette: https://vast-lib.github.io/tinyVAST/articles/mgcv.html
sim <- replicate(n = 100, expr = fit$obj$simulate()$y_i)

res <- DHARMa::createDHARMa(simulatedResponse = sim, 
                            observedResponse = dat$cpue, 
                            fittedPredictedResponse = fitted(fit))

# q-q plot
png(file = here(workDir, "results_age", "qq.png"),
    width = 130, height = 130, units = "mm", res = 300)
DHARMa::plotQQunif(simulationOutput = res, 
                   testOutliers = FALSE, 
                   testUniformity = FALSE,
                   testDispersion = FALSE)
dev.off()

# Spatial residuals
map_list <- list()
if (!dir.exists(here(workDir, "results_age", "spatial_residuals"))) {
  dir.create(here(workDir, "results_age", "spatial_residuals"))
}

if(species == "pacific_cod"){ 
    for(i in min(ages):max(ages)) {
    df <- cbind.data.frame(dat, residuals = res$scaledResiduals) %>%
      filter(age == i) 
    map <- ggplot() +
      geom_point(data = df, aes(x = X, y = Y, color = residuals), shape = 15, size = 0.9) +
      scale_color_gradient2(low = "darkred", mid = "white", high = "darkblue", midpoint = 0.5) +
      xlab("eastings") + ylab("northings") + ggtitle(paste0("age ", i)) +
      facet_wrap(~year)
    map_list[[i+1]] <- map # accounting for age 0s in cod for indexing
    ggsave(map, filename = here(workDir, "results_age", "spatial_residuals", paste0("age_", i, ".png")),
           width = 300, height = 300, units = "mm", dpi = 300)
    }
  }else{
    for(i in min(ages):max(ages)) {
      df <- cbind.data.frame(dat, residuals = res$scaledResiduals) %>%
        filter(age == i) 
      map <- ggplot() +
        geom_point(data = df, aes(x = X, y = Y, color = residuals), shape = 15, size = 0.9) +
        scale_color_gradient2(low = "darkred", mid = "white", high = "darkblue", midpoint = 0.5) +
        xlab("eastings") + ylab("northings") + ggtitle(paste0("age ", i)) +
        facet_wrap(~year)
      map_list[[i]] <- map
      ggsave(map, filename = here(workDir, "results_age", "spatial_residuals", paste0("age_", i, ".png")),
             width = 300, height = 300, units = "mm", dpi = 300)
    }
    }

# Compare proportions to previous model ---------------------------------------
# TODO: change this path to the previous run
previous_props <- read.csv(here("species_specific_code", "BS", species, 
                                "archive",
                                "2025",
                                "hindcast", 
                                "results_age", 
                                "tinyVAST_props.csv"))
previous_name <- "2025 hindcast tinyVAST"  # TODO: name the previous run

# # Reshape VAST output to match tinyVAST output
# tiny_years <- c(1980:2019, 2021:this_year)
# old_props <- old_props %>% filter(Year %in% tiny_years & Region == "Both")
# old_props <- old_props[, c(16, 1:15, 17)]
# colnames(old_props)[c(1, 16, 17)] <- c("year", "age_15", "region")

compare_props <- function(dfs, names) {
  # Reshape dataframes for plotting summaries of the two models
  df <- data.frame()
  for(i in 1:length(dfs)) {
    prop <- dfs[[i]]
    colnames(prop) <- c("year", ages, "region")
    prop <- prop %>% select(-region)  # Be careful if the regions are different.
    prop <-  reshape2::melt(prop, id.vars = "year",
                 variable.name = "age", value.name = "proportion")
    prop$version <- names[i]
    df <- bind_rows(df, prop)
  }
  
  # side-by-side barplot of proportion-at-age for each year for each model
  barplot <- ggplot(df, aes(x = age, y = proportion, fill = version)) +
    geom_bar(stat = "identity", position = "dodge") +
    viridis::scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_x_discrete(breaks = pretty(ages)) + 
    xlab("Age") + ylab("Proportion-at-Age") +
    theme(legend.title = element_blank()) +
    facet_wrap(~ year, ncol = 6, dir = "v")
  barplot
  
  # Boxplot summarizing proportion-at-age across years for each model
  boxplot <- ggplot(df, aes(x = age, y = proportion, color = version, fill = version)) +
    geom_boxplot(alpha = 0.5) +
    viridis::scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    viridis::scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
    scale_x_discrete(breaks = pretty(ages)) +  
    xlab("Age") + ylab("Proportion-at-Age") +
    theme(legend.title = element_blank()) 

  # Reshape data for percent difference between the models
  new <- dfs[[1]] %>% select(-region)
  old <- dfs[[2]] %>% select(-region)
  new <- new %>% filter(year <= max(old$year))  # make sure new dataset is the same length as the old
  
  # Get percent difference
  per <- round((((new[, -1] - old[, -1]) / old[, -1]) *100), 4)
  per_tab <- cbind(per, year = new[, 1])
  
  # Plot the percent difference between the models
  colnames(per_tab)[min(ages):length(ages)] <- min(ages):max(ages)
  diff <-  reshape2::melt(per_tab, id.vars = "year", 
               variable.name = "age", value.name = "proportion") %>%
    # add column for coloring the bars in the plot based on positive/negative
    mutate(sign = case_when(proportion >= 0 ~ "positive",
                            proportion < 0 ~ "negative"))
  
  # Plot both regions together and without 2020
  label <- paste0("Percent difference between ", names[1], " and ", names[2])
  diff_plot <- ggplot(diff, aes(x = age, y = proportion, fill = sign)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = c("cornflowerblue", "darkred")) +
    scale_x_discrete(breaks = pretty(ages)) + 
    xlab("Age") + ylab(label) +
    facet_wrap(~ year, ncol = 6, dir = "v") 
  
  return(list(barplot = barplot, boxplot = boxplot, diff_plot = diff_plot))
}

prop_diff <- compare_props(dfs = list(props, previous_props),
                           names = c(phase, previous_name))
prop_diff$barplot
prop_diff$boxplot
prop_diff$diff_plot

# Save comparison
ggsave(prop_diff$barplot, filename = here(workDir, "results_age", "comps_by_year.png"),
       width=220, height=200, units="mm", dpi=300)
ggsave(prop_diff$boxplot, filename = here(workDir, "results_age", "comps_summary.png"),
       width=200, height=120, units="mm", dpi=300)
ggsave(prop_diff$diff_plot, filename = here(workDir, "results_age", "comp_per_diff.png"),
       width=200, height=200, units="mm", dpi=300)
