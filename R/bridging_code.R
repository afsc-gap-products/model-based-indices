# Bridging code to check between changes in data set
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2022.08.24
# Date updated: 2022.08.24

options(scipen=999)

# library -----------------------------------------------------------------
library(tidyverse)
library(here)
library(janitor)
library(gt)

# folders -----------------------------------------------------------------

dir.create(here::here("output"), showWarnings = F)
dir.create(here::here("output", "bridging"), showWarnings = F)

# functions ---------------------------------------------------------------

functions <- list.files(here::here("functions"))
purrr::walk(functions, ~ source(here::here("functions", .x)))

# raw data ----------------------------------------------------------------

new_alk <- readRDS(here::here("data", "BS_Pacific_Cod", "unstratified_alk.RDS"))

old_alk <- readRDS(here::here("data", "BS_Pacific_Cod", "unstratified_alk_2021.RDS"))

new_agecomps <- read_csv(here::here("bridging", "Data_Geostat_2022Pacific_Cod_Age_hindcast_check.csv")) #need to get age comps from VAST code

old_agecomps <- read_csv(here::here("bridging", "Data_Geostat_2021Pacific_Cod_Age_hindcast_check.csv")) #need to get age comps from VAST code; need to load combined EBS NBS


# clean data --------------------------------------------------------------

n_alk <- new_alk %>% 
  as_tibble() %>% 
  clean_names() %>% 
  rename(probability_new = probability)
  # mutate(source = 'new_alk')

o_alk_ebs <- old_alk$EBS %>% 
  as_tibble() %>% 
  clean_names() %>% 
  mutate(region = "EBS")

o_alk_nbs <- old_alk$NBS %>% 
  as_tibble() %>% 
  clean_names() %>% 
  mutate(region = "NBS")

o_alk <- bind_rows(o_alk_ebs, o_alk_nbs) %>% 
  arrange(length, age, sex, year) %>% 
  rename(probability_old = probability) 
  # mutate(source = "old_alk")

n_agecomps <- new_agecomps %>% 
  as_tibble() %>% 
  clean_names()

o_agecomps <- old_agecomps %>% 
  as_tibble() %>% 
  clean_names()


# * filter for species ----------------------------------------------------

n_alk <- n_alk %>% 
  dplyr::filter(species_code == 21720)
o_alk <- o_alk %>% 
  dplyr::filter(species_code == 21720)
# n_agecomps <- n_agecomps %>% 
#   dplyr::filter(species_code == 21720)
# o_agecomps <- o_agecomps %>% 
#   dplyr::filter(species_code == 21720)


# duplicate keys check ----------------------------------------------------

dupes_new <- n_alk %>% 
  # dplyr::select(year, sex, length, age, region) %>% 
  duplicated() 
table(dupes_new) #288
n_alk <- n_alk %>% bind_cols(dupes = dupes_new) %>% 
  dplyr::filter(dupes == FALSE) %>% 
  dplyr::select(-dupes)

dupes_old <- o_alk %>% 
  # dplyr::select(year, sex, length, age, region) %>%
  duplicated() 
table(dupes_old) #576
o_alk <- o_alk %>% bind_cols(dupes = dupes_old) %>% 
  dplyr::filter(dupes == FALSE) %>% 
  dplyr::select(-dupes)


# sample size table -------------------------------------------------------

# UPDATE: add table for sanity check to confirm same numbers tows at age


# alk ---------------------------------------------------------------------


# * tables: check differences alk -----------------------------------------------
check_alk <- full_join(n_alk, o_alk) %>% #, by = c( "species_code", "length", "age", "sex", "year", "region")) %>% 
  mutate(abs_diff = round(abs(probability_new - probability_old), 4),
         abs_diff = if_else(is.na(probability_new), probability_old, abs_diff),
         abs_diff = if_else(is.na(probability_old), probability_new, abs_diff))

check_alk2 <- bind_rows(n_alk %>% mutate(source = "new_alk", probability = probability_new), 
                        o_alk %>% mutate(source = "old_alk", probability = probability_old)) %>% 
  dplyr::select(-probability_new, -probability_old)

keys_affected <- check_alk %>% filter(abs_diff > 0)

range(keys_affected$year)
unique(keys_affected$year) %>% sort() #every year
range(keys_affected$abs_diff)         #0.0001 to 0.5

# check for duplicated combos
dupes <- check_alk %>% dplyr::select(year, sex, length, age, region) %>% duplicated() 
table(dupes)
bind_cols(check_alk, dupes = dupes) %>% 
  dplyr::filter(dupes == TRUE) %>% 
  arrange(year, region, sex, length, age) #%>% View()
  
alk_check_table <- check_alk %>% 
  rename(prob_new = probability_new, prob_old = probability_old) %>% 
  dplyr::select(year, length, age, sex, region, prob_old, prob_new) %>% 
  group_by(sex, length, age) %>% 
  summarise(mean_prob_new = mean(prob_new, na.rm = T),
            mean_prob_old = mean(prob_old, na.rm = T)) %>% 
  # dplyr::filter(!is.na(mean_prob_new), !is.na(mean_prob_old),
                # mean_prob_new > 0, mean_prob_old >0) %>% 
  mutate(diff_age = abs(mean_prob_new - mean_prob_old),
         diff_age = if_else(is.na(mean_prob_new), mean_prob_old, diff_age),
         diff_age = if_else(is.na(mean_prob_old), mean_prob_new, diff_age)) %>% 
  rename(new_mean_prob_age = mean_prob_new,
         old_mean_prob_age = mean_prob_old) %>% 
  pivot_wider(names_from = age, values_from = c(new_mean_prob_age, old_mean_prob_age, diff_age), 
              names_sort = TRUE, names_vary = "slowest") #%>% 
  # gt()
write_csv(alk_check_table, file = here("output", "bridging", "alk_comparison_table.csv"))

alk_check_table_yrs <- check_alk %>% 
  dplyr::select(year, length, age, sex, region, probability_new, probability_old, abs_diff) %>% 
  group_by(year, sex, length, age, region) %>% 
  rename(new_prob_age = probability_new,
         old_prob_age = probability_old,
         abs_diff_age = abs_diff) %>% 
  pivot_wider(names_from = age, values_from = c(new_prob_age, old_prob_age, abs_diff_age), 
              names_sort = TRUE, names_vary = "slowest") %>%
  arrange(year, region, sex, length)
  # unlist()
# gt()
write_csv(alk_check_table_yrs, file = here("output", "bridging", "alk_comparison_table_yrs.csv"))

# * plots: check differences alk ------------------------------------------------

hist(check_alk$abs_diff)
hist(keys_affected$abs_diff)

for(i in unique(check_alk$sex))
{
  for(j in unique(check_alk$region))
  {
    # p0 <- check_alk %>% 
    #   dplyr::filter(sex == i, region == j) %>%   
    #   group_by(year, age, length) %>% 
    #   mutate(tf_diff = if_else(abs_diff > 0 && !is.na(abs_diff), TRUE, FALSE)) %>% 
    #   summarize(n_in_group = n(),
    #             n_affected = sum(tf_diff == TRUE),
    #             mean_prop_diff = mean(abs_diff)) %>% 
    # ggplot() +
    #   geom_col(aes(x = length, y = mean_prop_diff)) +
    #   # geom_bar(aes(x = length, y = n_affected), stat = "identity") +
    #   facet_wrap(age~year)
    #   facet_grid(age ~ year)
    #   # facet_grid(rows = vars(age), cols = vars(year))
    # ggsave(p0, filename = paste0("props_affected_sex_", i, "_region_",j, ".png"), 
    #        path = here("output", "bridging"), width = 20, height = 20)
    # 
    p1 <- check_alk %>% 
      dplyr::filter(sex == i, region == j) %>%   
      group_by(age, length) %>% 
      mutate(tf_diff = if_else(abs_diff > 0 && !is.na(abs_diff), TRUE, FALSE)) %>% 
      summarize(n_in_group = n(),
                n_affected = sum(tf_diff == TRUE),
                mean_prop_diff = mean(abs_diff)) %>% 
      ggplot() +
      geom_col(aes(x = length, y = mean_prop_diff)) +
      facet_wrap(~age) +
      labs(x = "mean proportion difference", 
           title = "ALK: Mean absolute difference in proportions at length and age",
           subtitle = "panels = age")
    ggsave(p1, filename = paste0("alk_mean_prop_change_sex", i, "_region_",j, ".png"), 
           path = here("output", "bridging"),
           height = 10, width = 10)
    
    p2 <- check_alk2 %>% 
      dplyr::filter(sex == i, region == j) %>%   
      group_by(age, length, source) %>% 
      summarize(mean_prop = mean(probability, na.rm = TRUE)) %>% 
      ggplot() +
      geom_col(aes(x = length, y = mean_prop, fill = source), 
               alpha = 0.6, position = "identity") +
      scale_fill_manual(values = c("black", "goldenrod1"))+
      # geom_col(aes(x = length, y = mean_old_prop), fill = "darkblue") +      # geom_col(aes(x = length, y = mean_old_prop), color = "orange", alpha = 0.1) +
      facet_wrap(~age)+
      labs(x = "mean proportion", 
           title = "ALK: Mean proportion at length and age",
           subtitle = "panels = age")
    # p2
    ggsave(p2, filename = paste0("alk_mean_prop_bydata_sex", i, "_region_",j, ".png"), 
           path = here("output", "bridging"),
           height = 10, width = 10)
    
  }
}

# UPDATE: tile plot of prob of age at length
# # 1) prob age at length new
# # 2) prob age at length old


# plot_alk <- bind_rows(n_alk, o_alk)

# need raw specimen data for this, not alk
# for(i in 1:unique(plot_alk$sex))
# {
#   p0 <- plot_alk %>% 
#     dplyr::filter(sex == i) %>% 
#     ggplot() + 
#     geom_point(aes(x = length, y = age, colour = source))
#   p1 <- plot_alk %>% 
#     dplyr::filter(sex == i) %>% 
#     ggplot() + 
#     geom_point(aes(x = age, y = length, colour = source))
#   
#   ggsave(p0, filename = paste0("age_at_length_sex_", i, ".png"), path = here("output", "bridging"))
#   ggsave(p1, filename = paste0("length_at_age_sex_", i, ".png"), path = here("output", "bridging"))
# }
# 
# for(i in 1:unique(plot_alk$sex))
# {
#   p0 <- plot_alk %>%
#     dplyr::filter(sex == i) %>% #, source == "new_alk") %>%
#     group_by(age)
#     summarize(mean_laa = mean(as.numeric(length))) %>% 
#     ggplot() +
#     geom_point(aes(x = length, y = probability, colour = factor(age))) +
#       facet_wrap(~year)
#   p1 <- plot_alk %>%
#     dplyr::filter(sex == i) %>%
#     ggplot() +
#     geom_point(aes(x = age, y = probability, colour = factor(length)))
# 
#   ggsave(p0, filename = paste0("age_at_length_sex_", i, ".png"), path = here("output", "bridging"))
#   ggsave(p1, filename = paste0("length_at_age_sex_", i, ".png"), path = here("output", "bridging"))
# }


# age comps ---------------------------------------------------------------

n_agec <- n_agecomps %>% 
  group_by(year, age) %>% 
  summarize(catch_tot = sum(catch_kg, na.rm = T),
            catch_mean = mean(catch_kg, na.rm = T)) %>% 
  mutate(comp_source = "alk_2022")

o_agec <- o_agecomps %>% 
  group_by(year, age) %>% 
  summarize(catch_tot = sum(catch_kg, na.rm = T),
            catch_mean = mean(catch_kg, na.rm = T)) %>% 
  mutate(comp_source = "alk_2021")
  

age_comp_set <- bind_rows(n_agec, o_agec)


# * tables ----------------------------------------------------------------
agecomp_tab <- full_join(
  n_agecomps %>% rename(catch_kg_new = catch_kg),
  o_agecomps %>% rename(catch_kg_old = catch_kg)
) %>% 
  mutate(abs_diff_catch = round(abs(catch_kg_new - catch_kg_old), 4))

write_csv(agecomp_tab, file = here("bridging", "full_diff_agecomps.csv"))

agecomp_tab_age <- agecomp_tab %>% 
  group_by(age) %>% 
  summarise(#mean_of_abs_diffs = mean(abs_diff_catch),
            diff_mean_catch_at_age = round(abs(mean(catch_kg_new) - mean(catch_kg_old)), 4))
write_csv(agecomp_tab, file = here("bridging", "diff_btw_mean_catch_at_age.csv"))

# * figures ---------------------------------------------------------------

p0 <- ggplot(age_comp_set, aes(x = year, y = catch_tot, color = comp_source)) +
  geom_line()+
  facet_wrap(~ age, ncol = 5, scales = "free") +
  theme_bw()

p0

ggsave(p0, filename = paste0("total_catch_at_age.png"), path = here("bridging"),
       height = 5, width = 15)

p0_1 <- ggplot(age_comp_set, aes(x = year, y = catch_mean, color = comp_source)) +
  geom_line()+
  facet_wrap(~ age, ncol = 5, scales = "free") +
  theme_bw()

p0_1

ggsave(p0_1, filename = paste0("mean_catch_at_age.png"), path = here("bridging"),
       height = 5, width = 15)

max_age <- max(age_comp_set$age)

for(i in 0:max_age)
{
  p1 <- age_comp_set %>% 
    dplyr::filter(age == i) %>% 
    ggplot(aes(x = year, y = catch_tot, color = comp_source)) +
    geom_line()+
    # facet_wrap(~ age, ncol = 5, scales = "free") +
    theme_bw() +
    ggtitle(paste("Age", i)) +
    ylab("total catch")
  ggsave(p1, filename = paste0("total_catch_age_", i, ".png"), 
         path = here("bridging"), height = 10, width = 13)
  
  p2 <- age_comp_set %>%
    dplyr::filter(age == i) %>%
    ggplot(aes(x = year, y = catch_mean, color = comp_source)) +
    geom_line()+
    # facet_wrap(~ age, ncol = 5, scales = "free") +
    theme_bw()+
    ggtitle(paste("Age", i)) +
    ylab("mean catch")
  ggsave(p2, filename = paste0("mean_catch_age_", i, ".png"), 
         path = here("bridging"), height = 10, width = 13)

}

# sample size
age_comp_set %>% 
  janitor::tabyl(age, comp_source, year) 

p3 <- age_comp_set %>% 
  group_by(year, comp_source) %>% 
  summarise(n = n()) %>% 
  ggplot() +
  geom_line(aes(x = year, y = n, color = comp_source, linetype = comp_source)) +
  theme_bw() +
  ylab("sample size") #+
# scale_linetype_manual(values = c("dotted", "twodash", "longdash", "solid"))
ggsave(p3, filename = paste0("sample_size_comparison.png"), path = here("bridging"))


# * proportions -----------------------------------------------------------



# * index.csv -------------------------------------------------------------

# sanity check

