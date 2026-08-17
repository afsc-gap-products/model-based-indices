#' Test of an age composition model using sdmTMB, following a vignette for 
#' multispecies models: https://github.com/pbs-assess/sdmTMB/blob/main/vignettes/web_only/multispecies.Rmd
#' and a population composition example: https://github.com/pbs-assess/sdmTMB/blob/comps/scratch/composition-example.R
# By: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.12.09

library(here)
library(sdmTMB)
library(dplyr)
library(sp)
library(ggplot2)
library(viridis)
# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
theme_set(theme_sleek())

# Set species -----------------------------------------------------------------
species <- 21740
# this_year <- lubridate::year(Sys.Date())
this_year <- 2024  # different year for debugging
Species <- "pollock"
phase <- c("hindcast","production")[1]
speciesName <- paste0("Walleye_Pollock_age_", as.character(this_year), "_EBS-NBS")
workDir <- here::here("species_specific_code", "BS", Species, phase)
Data <- read.csv(here("species_specific_code", "BS", Species, phase, "data", 
                      paste0("VAST_ddc_alk_", this_year, ".csv")))

# Set up sdmTMB model ---------------------------------------------------------
# Re-format data for sdmTMB
dat <-  transmute(Data,
                  cpue = CPUE_num, # converts cpue from kg/ha to kg/km^2
                  year = as.integer(Year),
                  age = as.integer(Age),
                  Y = Lat,
                  X = Lon) %>% 
  as.data.frame() %>%  # ensure not a tibble
  filter(year != 2020)  # drop dummy 2020 data

dat$year_f <- as.factor(dat$year)
dat$age_f <- as.factor(dat$age)

# Transform lat & lon to UTM
coordinates(dat) <- ~ X + Y
proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")
dat <- as.data.frame(spTransform(dat, CRS("+proj=utm +zone=2")))
# scale to km so values don't get too large
dat$X <- dat$coords.x1 / 1000
dat$Y <- dat$coords.x2 / 1000

# use prior model mesh from VAST
mesh <- make_mesh(dat, xy_cols = c("X", "Y"), 
                  mesh = readRDS(file = here("meshes", "bs_vast_mesh_50_knots.RDS")))

# Fit sdmTMB model ------------------------------------------------------------
fit_sdmTMB <- sdmTMB( 
  cpue ~ 0 + year_f * age_f,
  data = dat, 
  mesh = mesh,
  family = delta_gamma(type = "poisson-link"), 
  time = "year", 
  spatial = "on",
  spatiotemporal = "ar1",
  extra_time = 2020L, 
  silent = FALSE,
  anisotropy = FALSE,
  #control = sdmTMBcontrol(nlminb_loops = 1L, newton_loops = 2L),
  do_fit = TRUE
)

# apply additional optimization loops as needed to reduce gradients
# fit_sdmTMB <- run_extra_optimization(fit_sdmTMB, 
#                                      nlminb_loops = 0, 
#                                      newton_loops = 1)

saveRDS(fit_sdmTMB, file = here("species_specific_code", "BS", Species, 
                                phase, "results_age", "fit_sdmTMB_age.RDS"))

# Diagnostic plots ------------------------------------------------------------
sanity(fit_sdmTMB)

# q-q plot
pdf(file = here("species_specific_code", "BS", Species, phase, "results_age", 
                "qq.pdf"), width = 5, height = 5)
sims <- simulate(fit_sdmTMB, nsim = 500, type = "mle-mvn") 
sims |> dharma_residuals(fit_sdmTMB, test_uniformity = FALSE)
dev.off()


# Get abundance density indices for each year-age combination -----------------
fit_sdmTMB <- readRDS(here("species_specific_code", "BS", Species, phase, 
                           "results_age",  "fit_sdmTMB_age.RDS"))

# Load in grid for EBS, NBS, or entire Bering Sea
load_grid <- function(region) {
  if(region == "EBS") {
    grid <- read.csv(here("extrapolation_grids", "ebs_coarse_grid.csv"))
  }
  if(region == "NBS") {
    grid <- read.csv(here("extrapolation_grids", "nbs_coarse_grid.csv"))
  }
  if(region == "BS") {
    grid <- read.csv(here("extrapolation_grids", "bering_coarse_grid.csv"))
  } 
  
  # replicate extrapolation grids for each year in data
  pred_grid <- replicate_df(grid, "year_f", unique(dat$year_f))
  pred_grid$year <- as.integer(as.character(factor(pred_grid$year_f)))
  
  return(pred_grid)
}

pred_grid <- load_grid(region = "BS")

# Loop over ages - getting an index for each age
ages <- unique(dat$age_f)
ind_list <- lapply(ages, \(a) {
  print(a)
  pred_grid$age_f <- a
  pred <- predict(fit_sdmTMB, newdata = pred_grid, return_tmb_object = TRUE)
  ind <- get_index(obj = pred, area = pred_grid$area_km2, bias_correct = TRUE)
  data.frame(ind, Age = a)
})
ind <- do.call(rbind, ind_list)

# Plot age-structured abundance (index)
ggplot(ind, aes(year, est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray") +
  geom_line() +
  facet_wrap(~Age, scales = "fixed")

write.csv(ind, here(workDir, "results_age",  "sdmTMB_age_abundance.csv"))

# Calculate proportion-at-age -------------------------------------------------
prop <- ind %>% 
  group_by(year) %>%
  mutate(total = sum(est)) %>%
  ungroup() %>%
  group_by(year, Age) %>%
  summarize(proportion = est / total)

# Add column for cohort colors 
colors <- rep(1:16, length(1982:2024))
prop$color <- colors[1:nrow(prop)]

ggplot(prop, aes(x = Age, y = proportion, fill = color)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis(option = "turbo") +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  scale_y_continuous(limits = c(0, 0.5), breaks = c(0, 0.2, 0.4)) +
  ylab("Proportion") + 
  guides(fill = "none") +
  facet_wrap(~ year, ncol = 4, dir = "v") +
  theme(strip.text.x = element_blank()) +
  geom_text(x = 13, y = 0.45, aes(label = year), color = "grey30", size = 2.8)

write.csv(prop, here(workDir, "results_age", "sdmTMB_age_prop.csv"))

# Compare to production VAST model --------------------------------------------
# Read in and reshape VAST comps
vast_props <- read.csv(here(workDir, "results_age", "proportions.csv"))
colnames(vast_props)[1:15] <- 1:15
vast_props <- reshape2::melt(vast_props, id.vars = c("Year", "Region"),
                             variable.name = "Age", value.name = "Proportion") %>%
  filter(Region == "Both") 

# Combine with sdmTMB comps
all_props <- rbind(cbind(year = vast_props$Year, 
                         age = vast_props$Age, 
                         proportion = vast_props$Proportion, 
                         model = "VAST"),
                   cbind(year = prop$year, 
                         age = prop$Age, 
                         proportion = prop$proportion, 
                         model = "sdmTMB")) 
all_props <- data.frame(all_props) %>%
  mutate(year = as.integer(year),
         age = factor(as.integer(age)),
         proportion = as.numeric(proportion),
         model = factor(model))

props_compare <- ggplot(all_props, aes(x = age, y = proportion, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  ylab("Proportion-at-age") +
  facet_wrap(~ year, ncol = 6, dir = "v")
props_compare

props_summary <- ggplot(all_props, aes(x = age, y = proportion, color = model, fill = model)) +
  geom_boxplot(alpha = 0.5) +
  scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
  scale_fill_viridis(discrete = TRUE, option = "plasma", end = 0.9) +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  ylab("Proportion-at-age") 
props_summary

ggsave(props_compare, filename = here(workDir, "results_age", "sdmTMB_compare.png"),
       width=200, height=180, units="mm", dpi=300)
ggsave(props_summary, filename = here(workDir, "results_age", "sdmTMB_summary.png"),
       width=200, height=120, units="mm", dpi=300)
