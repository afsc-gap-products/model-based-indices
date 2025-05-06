
# Install packages if needed
if (!requireNamespace("gapindex", quietly = TRUE)) {
  devtools::install_github("afsc-gap-products/gapindex")
}

# Connect to Oracle -----------------------------------------------------------
# Load Oracle credentials file, if you've made one, or enter & save credentials
if(file.exists("~/oracle_credentials.R")) { 
  source("~/oracle_credentials.R")
} else {
  oracle_user <- rstudioapi::showPrompt(title = "Username",
                                        message = "Oracle Username",
                                        default = "")
  oracle_pw <- rstudioapi::showPrompt(title = "Password",
                                      message = "Oracle Password",
                                      default = "")
}

# Two different options for connecting to Oracle
channel <- RODBC::odbcDriverConnect(
  connection = paste0("Driver=/opt/oracle/instantclient_12_2/libsqora.so.12.1;DBQ=raja.afsc.noaa.gov:1521/afsc;UID=", 
                      oracle_user, ";PWD=", oracle_pw),
  rows_at_time = 1
)

# Calculate Center of Gravity -------------------------------------------------
calc_cpue <- function(yr_goa) {
  ## Define species and species groupings
  rf_groups <- data.frame(
    GROUP_CODE   = c(30060, 30050, 30050, 30050, 30576, 30420, 30152, 30020),
    SPECIES_CODE = c(30060, 30050:30052, 30576, 30420, 30152, 30020)
  )
  
  ## Pull data
  gp_data <- 
    gapindex::get_data(year_set = c(seq(from = 1990, to = 1999, by = 3),
                                    seq(from = 2003, to = yr_goa, by = 2)),
                       survey_set = "GOA",
                       spp_codes = rf_groups,
                       channel = channel
    )
  
  ## Calculate cpue
  gp_cpue <- gapindex::calc_cpue(gapdata = gp_data) |> as.data.frame()
  
  
  calc_weighted_mean <- function(x, w, lwr_p = 0.025, upr_p = 0.975) {
    ## Count the number of records that have a positive weight (w) 
    ## and a metric (some metrics are missing from some hauls) 
    n <- length(x = w[w > 0 & !is.na(x = w) & !is.na(x)])
    w <- w / sum(w, na.rm = TRUE) # weight add up to 1
    
    ## Calculate weighted mean and standard error
    weighted_mean <- weighted.mean(x = x, 
                                   w = w, 
                                   na.rm = TRUE)
    
    weighted_var <- sum(w * ((x - weighted_mean)^2), na.rm = TRUE) * (n / (n - 1))
    weighted_se <- sqrt(weighted_var / n)
    
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
  
  return(cogs)
}

cogs <- calc_cpue(2023)

# Plot ------------------------------------------------------------------------
plot_cogs <- function(df) {
  if (!requireNamespace("ggsidekick", quietly = TRUE)) {
    devtools::install_github("seananderson/ggsidekick")
  }
  library(ggsidekick)
  theme_set(theme_sleek())
  
  # Special color palette for time series plot!
  if (!requireNamespace("nationalparkcolors", quietly = TRUE)) {
    devtools::install_github("katiejolly/nationalparkcolors")
  }
  
  # Update dataframe to have common names, cleaner labels, correct axis
  cogs_plot <- df %>%
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
      metric == "BOTTOM_TEMPERATURE_C" ~ "Bottom Temp (C)",
      metric == "DEPTH_M" ~ "Depth (m)",
      metric == "LATITUDE_DD_START" ~ "Latitude",
      metric == "LONGITUDE_DD_START" ~ "Longitude"
    )) %>%
    # Remove first two dusky points (1990 & 1993; data should start in 1996)
    filter(!(species_code == "Dusky Rockfish" & year <= 1996)) %>% 
    # Make depth estimates negative for inverted axis when plotting
    mutate(est = if_else(metric == "Depth (m)", -est, est),
           upr = if_else(metric == "Depth (m)", -upr, upr),
           lwr = if_else(metric == "Depth (m)", -lwr, lwr))
  
  # Time series plot ------------------------------------------------------------
  # Transform latitude & longitude to UTM (for estimates, upper & lower bounds)
  utm_transform <- function(column) {
    utm <- data.frame(Y = cogs_plot[cogs_plot$metric == "Latitude", column],
                      X = cogs_plot[cogs_plot$metric == "Longitude", column],
                      species_code = cogs_plot[cogs_plot$metric == "Longitude", "species_code"], 
                      year = cogs_plot[cogs_plot$metric == "Longitude", "year"])
    
    sp::coordinates(utm) <- ~ X + Y
    sp::proj4string(utm) <- sp::CRS("+proj=longlat +datum=WGS84")
    utm <- as.data.frame(sp::spTransform(utm, sp::CRS("+proj=utm +zone=5 +units=km")))
    colnames(utm)[c(3:4)] <- c("Eastings (km)", "Northings (km)")
    
    utm <- reshape2::melt(utm, 
                          id.vars = c("species_code", "year"), 
                          variable.name = "metric",
                          value.name = column)
    return(utm)
  }
  
  utm_out <- cbind.data.frame(utm_transform("est"),
                              se = utm_transform("se")$se,
                              lwr = utm_transform("lwr")$lwr,
                              upr = utm_transform("upr")$upr)
  
  pal <- nationalparkcolors::park_palette("Saguaro")  # special color palette!
  
  ts_plot <- rbind.data.frame(cogs_plot %>% filter(!metric %in% c("Latitude", "Longitude")),
                              utm_out[, c(3, 1, 2, 4:7)]) %>%  # Combine original & UTM dataframes
    ggplot(., aes(x = year, y = est)) +
    geom_line(aes(color = species_code)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = species_code), alpha = 0.4) +
    xlab("Year") + ylab("Weighted Mean") +
    theme(legend.title = element_blank()) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(labels = function(est) abs(est)) +
    facet_wrap(~ metric, scales = "free_y") 
  
  return(ts_plot)
}

