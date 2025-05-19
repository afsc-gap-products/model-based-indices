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

sp <- 2 # specify species from species vector
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
dat <- st_as_sf(dat,
                 coords = c('lon','lat'),
                 crs = st_crs(4326))
dat <- st_transform(dat,
                     crs = st_crs("+proj=utm +zone=2 +units=km"))
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

# Inputs to tinyVAST ----------------------------------------------------------
# Constant AR1 spatio-temporal term across ages & different variances for each age
dsem <- "\n  "
for(i in min(ages):max(ages)) {
  dsem <- paste0(dsem,  "age_", i, " -> age_", i, ", 1, lag1\n")
}

# # TinyVAST mesh
# mesh <- fm_mesh_2d(loc = Data[,c("X","Y")],
#                    cutoff = 50)
# 
# Format mesh for tinyVAST (using a function from sdmTMB; same method as sdmTMB index bridging)
old_mesh <- sdmTMB::make_mesh(dat, 
                              xy_cols = c("X", "Y"), 
                              mesh = readRDS(file = here("meshes/bs_vast_mesh_50_knots.RDS")),
                              fmesher_func = fm_mesh_2d()) 

# Fit model -------------------------------------------------------------------
fit <- tinyVAST(
  data = dat,
  formula = cpue ~ 0 + year_age,  
  sem = "",
  dsem = dsem,
  family = setNames(lapply(ages, function(x) delta_gamma(type = "poisson-link")), 
                    paste0("age_", ages)),
  delta_options = list(delta_formula = ~ 0 + year_age), # 2nd linear predictor
  space_column = c("X", "Y"), 
  variable_column = "age_f",
  time_column = "year",
  distribution_column = "age_f",
  spatial_graph = old_mesh,
  control = tinyVASTcontrol(getsd = FALSE,
                            profile = c("alpha_j"),
                            trace = 0)
)

# Save fit object (create directory for results first, if it doesn't exist)
if (!dir.exists(here(workDir, "results_age"))) {
  dir.create(here(workDir, "results_age"))
}

saveRDS(fit, here(workDir, "results_age", "tinyVAST_fit.RDS"))

# Age composition expansion ---------------------------------------------------
# Load fit object if needed
if(!exists("fit")) {
  fit <- readRDS(here(workDir, "results_age", "tinyVAST_fit.RDS"))
}

get_abundance <- function(region) {
  # Read in coarsened extrapolation grid
  if(region == "EBS") {grid <- read.csv(here("extrapolation_grids", "ebs_coarse_grid.csv"))}
  if(region == "NBS") {grid <- read.csv(here("extrapolation_grids", "nbs_coarse_grid.csv"))}
  if(region == "both") {grid <-  read.csv(here("extrapolation_grids", "bering_coarse_grid.csv"))}
  
  N_jz <- expand.grid(age_f = fit$internal$variables, year = sort(unique(dat$year)))
  N_jz <- cbind(N_jz, "abundance" = NA, "SE" = NA)
  for(j in seq_len(nrow(N_jz))) {
    if (N_jz[j, "age_f"] == 1) {
      message("Integrating ", N_jz[j, "year"], " ", N_jz[j, "age_f"], ": ", Sys.time())
    }
    if(is.na(N_jz[j, "abundance"])) {
      newdata <- data.frame(grid, 
                            year = N_jz[j, "year"], 
                            age_f = N_jz[j, "age_f"])
      newdata$year_age <- paste(newdata$year, newdata$age_f, sep = ".")
      # Area-expansion
      index1 <- integrate_output(fit,
                                 area = grid$area_km2,
                                 newdata = newdata,
                                 apply.epsilon = TRUE,
                                 bias.correct = TRUE,
                                 intern = TRUE)
      N_jz[j, "abundance"] <- index1[3] / 1e9
    }
  }
  
  N_ct <- array(N_jz$abundance, 
                dim = c(length(fit$internal$variables), length(unique(dat$year))),
                dimnames = list(fit$internal$variables,sort(unique(dat$year))))
  return(N_ct)
}

abundance <- get_abundance(region = region)

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
  return(prop)
}

props <- calc_props(abundance)

# props <- rbind.data.frame(calc_props(ebs, "EBS"), 
#                           calc_props(nbs, "NBS"), 
#                           calc_props(both, "both"))

# Save proportions
write.csv(props, here(workDir, "results_age", "tinyVAST_props.csv"), row.names = FALSE)

# Plot proportions with colors to track cohort strength
plot_proportions <- function(area = region, save_plot = TRUE) {
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
    geom_text(x = 13, y = 0.45, aes(label = year), color = "grey30", size = 2.8)
  
  return(plot)
  
  if(save_plot == TRUE) {
    ggsave(plot, filename = here(workDir, "results_age", paste0("Age_comp_", area, ".png")),
           width=120, height=180, units="mm", dpi=300)
  }
}

plot_proportions()

# Model diagnostics -----------------------------------------------------------
# working off of this vignette: https://vast-lib.github.io/tinyVAST/articles/mgcv.html
sim <- replicate(n = 100, expr = fit$obj$simulate()$y_i)

res <- DHARMa::createDHARMa(simulatedResponse = sim, 
                            observedResponse = dat$cpue, 
                            fittedPredictedResponse = fitted(fit))

res_data <- residuals(res)

# q-q plot
if (!requireNamespace("qqplotr", quietly = TRUE)) {
  install.packages("qqplotr")
  install.packages("twosamples")  # needed to install this separately on VM
}
library(qqplotr)

qqplot <- ggplot(data.frame(resid = res_data), aes(sample = resid)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  xlim(0, 1) + ylim(0, 1)
# qqplot

ggsave(qqplot, filename = here(workDir, "results_age", "qq.png"),
       width = 130, height = 130, units = "mm", dpi = 300)

# Spatial residuals
map_list <- list()
if (!dir.exists(here(workDir, "results_age", "spatial_residuals"))) {
  dir.create(here(workDir, "results_age", "spatial_residuals"))
}
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
