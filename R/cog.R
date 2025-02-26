## Calculate Empirical Center of Gravity From Observed Data
## Zack Oyafuso

## Connect to Oracle
library(gapindex)

if (file.exists("Z:/Projects/ConnectToOracle.R")) {
  source("Z:/Projects/ConnectToOracle.R")
} else {
  # For those without a ConnectToOracle file
  channel <- gapindex::get_connected(check_access = FALSE)
}

## Define species and species groupings
rf_groups <- data.frame(
  GROUP_CODE   = c(30060, 30050, 30050, 30050, 30576, 30420, 30152, 30020),
  SPECIES_CODE = c(30060, 30050:30052, 30576, 30420, 30152, 30020)
)

## Pull data
gp_data <- 
  gapindex::get_data(year_set = c(seq(from = 1990, to = 1999, by = 3),
                                  seq(from = 2003, to = 2023, by = 2)),
                     survey_set = "GOA",
                     spp_codes = rf_groups,
                     channel = channel
  )

## Calculate cpue
gp_cpue <- gapindex::calc_cpue(gapdata = gp_data) |> as.data.frame()


calc_weighted_mean <- function(x, w, lwr_p = 0.025, upr_p = 0.975) {
  ## Count the number of non-zero weight records
  n <- length(x = w[w > 0 & !is.na(x = w)])
  
  ## Calculate weighted mean and standard error
  weighted_mean <- weighted.mean(x = x, 
                                 w = w, 
                                 na.rm = TRUE)
  
  weighted_se <- sqrt(x = weighted.mean(x = (x-weighted_mean)^2, 
                                        w = w, 
                                        na.rm = TRUE ) / n)
  
  ## Calculate confidence interval of the weighted mean estimate
  ci <- qnorm(p = c(lwr_p, upr_p), 
              mean = weighted_mean, 
              sd = weighted_se)
  
  return(data.frame(est = weighted_mean, 
                    se = weighted_se, 
                    lwr = ci[1], 
                    upr = ci[2])
         )
}

## Loop over metrics and calculate weighted means, SEs, and CIs 
## for each species and year
cogs <- data.frame()
for (imetric in c("DEPTH_M", "BOTTOM_TEMPERATURE_C", 
                  "LATITUDE_DD_START", "LONGITUDE_DD_START")){
  cogs <- 
    rbind(cogs,
          data.frame(
            metric = imetric,
            do.call(
              what = rbind,
              args = lapply(
                X = split(x = gp_cpue, 
                          f = list(gp_cpue$SPECIES_CODE, gp_cpue$YEAR)),
                FUN = function(df){
                  data.frame(cbind(species_code = unique(df$SPECIES_CODE),
                                   year = unique(df$YEAR),
                                   calc_weighted_mean(x = df[, imetric],
                                                      w = df$CPUE_KGKM2)))
                  
                }
              )
            )
          )
    ) 
}

write.csv(cogs, here::here("output", "rf_cogs.csv"), row.names = FALSE)

## Plotting
library(dplyr)
library(ggplot2)
library(viridis)
library(ggsidekick)
theme_set(theme_sleek())

cogs <- read.csv(here::here("output", "rf_cogs.csv"))

cogs_plot <- cogs %>%
  mutate(species_code = factor(species_code)) %>%
  mutate(species_code = case_when(
    species_code == 30020 ~ "Shortspine Thornyhead",
    species_code == 30050 ~ "Rougheye & Blackspotted",
    species_code == 30060 ~ "Pacific Ocean Perch",
    species_code == 30152 ~ "Dusky Rockfish",
    species_code == 30420 ~ "Northern Rockfish",
    species_code == 30576 ~ "Shortraker Rockfish",
  )) %>%
  mutate(metric = case_when(
    metric == "BOTTOM_TEMPERATURE_C" ~ "Bottom Temp (°C)",
    metric == "DEPTH_M" ~ "Depth (m)",
    metric == "LATITUDE_DD_START" ~ "Latitude",
    metric == "LONGITUDE_DD_START" ~ "Longitude"
  ))

## Plot all variables as a time-series
pal <- nationalparkcolors::park_palette("Saguaro")
ts_plot <- ggplot(cogs_plot, aes(x = year, y = est)) +
  geom_line(aes(color = species_code)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = species_code), alpha = 0.4) +
  xlab("Year") + ylab("Weighted Mean") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  facet_wrap(~ metric, scales = "free_y")
ts_plot

## Plot latitude & longitude together in a 'sparkleplot'
cog_lat <- cogs_plot[cogs_plot$metric == "Latitude", c(2:4, 6:7)]
colnames(cog_lat)[3:5] <- c("est_lat", "lwr_lat", "upr_lat")
cog_lon <- cogs_plot[cogs_plot$metric == "Longitude", c(2:4, 6:7)]
colnames(cog_lon)[3:5] <- c("est_lon", "lwr_lon", "upr_lon")

cog_sparkle <- cog_lat %>% left_join(cog_lon, by = c("species_code", "year"))

sparkle <- ggplot(data = cog_sparkle, aes(x = est_lon, y = est_lat, color = year)) +
  geom_point() +
  # With arrow
  # geom_segment(data = cog_error2 %>% filter(Year >= this_year - 10), 
  #              aes(x = X, y = Y, xend = X2, yend = Y2), 
  #              alpha = 0.8, arrow = arrow(length = unit(0.03, "npc"))) +
  # Without arrow
  # geom_segment(data = cog_error2 %>% filter(Year >= this_year - 10), 
  #              aes(x = X, y = Y, xend = X2, yend = Y2), 
  #              alpha = 0.8) +
  geom_errorbar(aes(ymin = lwr_lat, ymax = upr_lat, color = year), alpha = 0.4) +
  geom_errorbarh(aes(xmin = lwr_lon, xmax = upr_lon, color = year), alpha = 0.4) +
  scale_color_viridis(option = "plasma", discrete = FALSE, end = 0.9) +
  xlab("Longitude (°W)") + ylab("Latitude (°N)") +
  facet_wrap(~species_code)
sparkle

## Plot points on a map
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
sf::sf_use_s2(FALSE)  # turn off spherical geometry
map <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = cog_sparkle, aes(x = est_lon, y = est_lat, color = year), size = 1) +
  # geom_errorbar(aes(ymin = lwr_lat, ymax = upr_lat, color = year), alpha = 0.4) +
  # geom_errorbarh(aes(xmin = lwr_lon, xmax = upr_lon, color = year), alpha = 0.4) +
  coord_sf(xlim = c(-162, -135), ylim = c(54, 60), expand = FALSE) +
  scale_color_viridis(option = "plasma", discrete = FALSE, end = 0.9) +
  scale_x_continuous(breaks = c(-160, -140)) +
  scale_y_continuous(breaks = c(55, 60)) +
  labs(x = NULL, y = NULL) +
  facet_wrap(~species_code, ncol = 2)
map

ggsave(ts_plot, filename = here::here("output", "rf_cog_ts.png"), 
       width = 190, height = 100, unit = "mm", dpi = 300)
ggsave(sparkle, filename = here::here("output", "rf_cog_sparkle.png"), 
       width = 160, height = 100, unit = "mm", dpi = 300)
ggsave(map, filename = here::here("output", "rf_cog_map.png"), 
       width = 150, height = 100, unit = "mm", dpi = 300)
