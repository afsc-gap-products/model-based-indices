# compare between years

library(here)
library(tidyverse)

index_21 <- read_csv(here("VAST_results", "2021_Index_Table_for_SS3.csv")) %>% 
  dplyr::filter(Fleet == "Both")
  # dplyr::select(Year, Estimate_metric_tons)
index_22_h <- read_csv(here("Index.csv")) %>% 
  dplyr::filter(Stratum == "Stratum_1") %>% 
  mutate(index_mt = Estimate * 0.001) %>% 
  rename(Year = Time)

compare <- full_join(dplyr::select(index_21, Year, Estimate_metric_tons),
                     dplyr::select(index_22_h, Year, index_mt)) %>% 
  mutate(diff_mt = Estimate_metric_tons - index_mt) %>% 
  rename(index_both_21 = Estimate_metric_tons, index_both_22_h = index_mt)

write.csv(compare, here("VAST_results", "index_bridge_2022.csv"))
