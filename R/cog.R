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
library(ggsidekick)
theme_set(theme_sleek())

cogs <- read.csv(here::here("output", "rf_cogs.csv"))

cogs_plot <- cogs %>%
  mutate(species_code = factor(species_code)) %>%
  mutate(species_code = case_when(
    species_code == 30020 ~ "Shortspine Thornyhead",
    species_code == 30050 ~ "Rougheye & Blackspotted Rockfish",
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

ggplot(cogs_plot, aes(x = year, y = est)) +
  geom_line(aes(color = species_code), alpha = 0.8) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = species_code), alpha = 0.3) +
  xlab("Year") + ylab("Weighted Mean") +
  theme(legend.title = element_blank()) +
  facet_wrap(~ metric, scales = "free_y")

