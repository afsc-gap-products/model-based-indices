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
speciesName <- paste0("Walleye_Pollock_age_", as.character(this_year), "_EBS-NBS")
workDir <- here::here("species_specific_code", "BS", Species)
Data <- read.csv(here("species_specific_code", "BS", Species, "data", 
                      paste0("VAST_ddc_alk_", this_year, ".csv")))

# Set up sdmTMB model ---------------------------------------------------------
# Re-format data for sdmTMB
dat <-  transmute(Data,
                  cpue = CPUE_num, # converts cpue from kg/ha to kg/km^2
                  year = as.integer(Year),
                  age = as.integer(Age),
                  vessel = "missing",
                  effort = 1, # area swept is 1 when using CPUE instead of observed weight
                  Y = Lat,
                  X = Lon,
                  pass = 0) %>% 
  as.data.frame() %>%  # ensure not a tibble
  filter(year != 2020)  # drop dummy 2020 data

dat$year_f <- as.factor(dat$year)

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
  cpue ~ year_f * age,
  spatial_varying = ~ 0 + factor(age),
  data = dat, 
  mesh = mesh,
  family = delta_gamma(type = "poisson-link"), 
  time = "year", 
  spatial = "off",
  spatiotemporal = "iid",
  extra_time = 2020L, #  omit if dummy 2020 included in data
  silent = FALSE,
  anisotropy = TRUE,
  do_fit = TRUE
)
fit_sdmTMB
saveRDS(fit_sdmTMB, file = here("species_specific_code", "BS", Species, "results", "fit_sdmTMB_age.RDS"))

# Diagnostic plots ------------------------------------------------------------
sanity(fit_sdmTMB)

pdf(file = here("species_specific_code", "BS", Species, "qq.pdf"),
    width = 5, height = 5)
sims <- simulate(fit_sdmTMB, nsim = 500, type = "mle-mvn") 
sims |> dharma_residuals(fit_sdmTMB, test_uniformity = FALSE)
dev.off()


# Get abundance density indices for each year-age combination -----------------
fit_sdmTMB <- readRDS(here("species_specific_code", "BS", Species, "results", "fit_sdmTMB_age.RDS"))

# Load in grids
load(here("extrapolation_grids", "eastern_bering_sea_grid.rda"))
load(here("extrapolation_grids", "northern_bering_sea_grid.rda"))

# EBS grid
grid_ll_ebs <- as.data.frame(eastern_bering_sea_grid)
names(grid_ll_ebs) <- tolower(names(grid_ll_ebs))
grid_ll_ebs <- grid_ll_ebs %>% 
  rename(X = lon, Y = lat)
coordinates(grid_ll_ebs) <- ~ X + Y
proj4string(grid_ll_ebs) <- CRS("+proj=longlat +datum=WGS84")
grid_ebs <- as.data.frame(spTransform(grid_ll_ebs, CRS("+proj=utm +zone=2")))
grid_ebs$X <- grid_ebs$coords.x1 / 1000 # scale to km to work with smaller numbers
grid_ebs$Y <- grid_ebs$coords.x2 / 1000 

# NBS grid
grid_ll_nbs <- as.data.frame(northern_bering_sea_grid)
names(grid_ll_nbs) <- tolower(names(grid_ll_nbs))
grid_ll_nbs <- grid_ll_nbs %>% 
  rename(X = lon, Y = lat)
coordinates(grid_ll_nbs) <- ~ X + Y
proj4string(grid_ll_nbs) <- CRS("+proj=longlat +datum=WGS84")
grid_nbs <- as.data.frame(spTransform(grid_ll_nbs, CRS("+proj=utm +zone=2")))
grid_nbs$X <- grid_nbs$coords.x1 / 1000 # scale to km to work with smaller numbers
grid_nbs$Y <- grid_nbs$coords.x2 / 1000 

# Combined grid
grid <- bind_rows(grid_nbs, grid_ebs)

# replicate extrapolation grids for each year in data
pred_grid_ebs <- replicate_df(grid_ebs, "year_f", unique(dat$year_f))
pred_grid_nbs <- replicate_df(grid_nbs, "year_f", unique(dat$year_f))
pred_grid <- replicate_df(grid, "year_f", unique(dat$year_f))
pred_grid_ebs$year <- as.integer(as.character(factor(pred_grid_ebs$year_f)))
pred_grid_nbs$year <- as.integer(as.character(factor(pred_grid_nbs$year_f)))
pred_grid$year <- as.integer(as.character(factor(pred_grid$year_f)))

# Loop over ages - getting an index for each age
ages <- unique(dat$age)
ind_list <- lapply(ages, \(a) {
  print(a)
  pred_grid_ebs$age <- a
  pred <- predict(fit_sdmTMB, newdata = pred_grid_ebs, return_tmb_object = TRUE)
  ind <- get_index(obj = pred, area = pred_grid_ebs$area_in_survey_km2, bias_correct = TRUE)
  data.frame(ind, Age = a)
})
ind <- do.call(rbind, ind_list)

# Plot age-structured abundance (index)
ggplot(ind, aes(year, est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray") +
  geom_line() +
  facet_wrap(~Age, scales = "fixed")

write.csv(ind, here(workDir, "results", "sdmTMB_age_abundance.csv"))

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

write.csv(prop, here(workDir, "results", "sdmTMB_age_prop.csv"))

# Compare to production VAST model --------------------------------------------
# Read in and reshape VAST comps
vast_props <- read.csv(here(workDir, "results", "Comps", "proportions.csv"))
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

ggsave(props_compare, filename = here(workDir, "results", "sdmTMB_compare.png"),
       width=200, height=180, units="mm", dpi=300)
ggsave(props_summary, filename = here(workDir, "results", "sdmTMB_summary.png"),
       width=200, height=120, units="mm", dpi=300)
