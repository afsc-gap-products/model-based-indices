# compare between years

library(here)
library(tidyverse)
library(ggplot2)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

# index_21 <- read_csv(here("VAST_results", "2021_Index_Table_for_SS3.csv")) %>% 
#   dplyr::filter(Fleet == "Both")
#   # dplyr::select(Year, Estimate_metric_tons)
# index_22_h <- read_csv(here("Index.csv")) %>% 
#   dplyr::filter(Stratum == "Stratum_1") %>% 
#   mutate(index_mt = Estimate * 0.001) %>% 
#   rename(Year = Time)
# 
# compare <- full_join(dplyr::select(index_21, Year, Estimate_metric_tons),
#                      dplyr::select(index_22_h, Year, index_mt)) %>% 
#   mutate(diff_mt = Estimate_metric_tons - index_mt) %>% 
#   rename(index_both_21 = Estimate_metric_tons, index_both_22_h = index_mt)
# 
# write.csv(compare, here("VAST_results", "index_bridge_2022.csv"))

### Plot this year on top of hindcast -----------------------------------------
index_22 <- read.csv("Index_2022.csv")
index_22$version <- "2022 hindcast"
index <- read.csv("Index.csv")
index$version <- "2023"

indices <- rbind.data.frame(index_22, index)
colnames(indices)[6] <- "error"

index_comp <- ggplot(indices, aes(x = Time, y = (Estimate / 1000000000), color = version)) +
  geom_line(linetype = "dashed", alpha = 0.8) +
  geom_pointrange(aes(ymin = (Estimate / 1000000000) - (error / 1000000000),
                      ymax = (Estimate / 1000000000) + (error / 1000000000),
                      shape = version), alpha = 0.8) +
  scale_shape(solid = FALSE) +
  xlab("Year") + ylab("Index (Mt)") +
  facet_wrap(~ Stratum, ncol = 1)

ggsave(index_comp, filename = "index_comparison_2023.png",
       width=130, height=160, units="mm", dpi=300)
