
# Notes -------------------------------------------------------------------
# Formerly:  Run_comps_2019-09-27.R from Jim Thorson's email
# For 2021, running this code again, as wrapper functions have proven problematic.

#library(googledrive)
library(tidyverse)
library(rgdal)
library(VAST)
library(tictoc)

# Set species -------------------------------------------------------------

species_code <- 21720
species_name <- "Pacific_Cod_Age_hindcast_check"
species_data <- "BS_Pacific_Cod"
compare <- FALSE # If compare = TRUE, useing 2021 alk

# Set up folder to store species specific results
folder <- here::here("results",species_name)
dir.create(folder, recursive = TRUE, showWarnings = F)


# Settings ----------------------------------------------------------------

Version <- "VAST_v13_1_0"
Region <- c("Eastern_Bering_Sea","Northern_Bering_Sea")
Method <- "Mesh"
knot_method <- "grid"
grid_size_km <- 25
n_x <- 50 #n_x <- 50   # Specify number of stations (a.k.a. "knots")
FieldConfig <- c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig <- c("Beta1"=0, "Beta2"=0, "Epsilon1"=4, "Epsilon2"=4)
OverdispersionConfig <- c("Eta1"=0, "Eta2"=0)
ObsModel <- c(2,4)
Options <- c("Calculate_Range" = FALSE, "Calculate_effective_area" = FALSE, "treat_nonencounter_as_zero" = TRUE )
Aniso <- FALSE
Npool <- 100 #20
fine_scale <- TRUE
BiasCorr <- TRUE
max_cells <- 2000

strata.limits <- data.frame('STRATA'="All_areas")

# Make settings 
settings <- make_settings( 
                      n_x = n_x,
                      Region = Region,
                      purpose = "index2",
                      fine_scale = fine_scale,
                      strata.limits = strata.limits,
                      ObsModel = ObsModel,
                      FieldConfig = FieldConfig,
                      RhoConfig = RhoConfig,
                      OverdispersionConfig = OverdispersionConfig,
                      Options = Options,
                      use_anisotropy = Aniso,
                      Version = Version,
                      max_cells = max_cells,
                      knot_method = knot_method,
                      bias.correct = BiasCorr
)

strata_names = c("Both","EBS","NBS")

#### Explore the data ####

    # Age Composition ---------------------------------------------------------
    plus_group <- 12
    min_year <- 1994
    
    # Load age-length keys produced by sumfish
    #alk_all <- readRDS(here::here("data","unstratified_alk_2021.RDS") )
    alk <- readRDS(here::here("data",species_data,"unstratified_alk.RDS") )
    if(compare == TRUE)
    {
      alk_all <- readRDS(here::here("data",species_data,"unstratified_alk_2021.RDS") ) #old
      alk_ebs <- alk_all$EBS %>%
          filter(SPECIES_CODE == species_code) %>%
          mutate(REGION = "EBS")
      alk_nbs <- alk_all$NBS %>%
          bind_rows( filter(alk_ebs, YEAR == 2018) ) %>%   # Use EBS ALK for 2018 ad hoc sampling in NBS
          filter(SPECIES_CODE == species_code) %>%
          mutate(REGION = "NBS")

      alk <- bind_rows(alk_ebs, alk_nbs)
    }
    
    sizeComp <- readRDS(here::here("data",species_data,"EBS_NBS_SizeComp.RDS") ) %>% #readRDS(here::here("data","EBS_NBS_SizeComp.RDS") )
        dplyr::filter(YEAR >= min_year,
                      SPECIES_CODE == species_code,
                      !is.na(EFFORT)
        )

    haulData <- readRDS(here::here("data",species_data,"EBS_NBS_Index.RDS") ) %>% 
        dplyr::filter(YEAR >= min_year,
                      SPECIES_CODE == species_code,
                      !is.na(EFFORT)
        )
    
    # Get summary calculations - this also fills zeroes within year - add sex if generating sex-specific agecomps
    allCats <- expand.grid(HAULJOIN=unique(haulData$HAULJOIN), AGE = unique(alk$AGE[alk$AGE<=plus_group]), noAge = 0) %>%
    inner_join(haulData, by = c("HAULJOIN")) 
    
    
    # Aggregate by Age key
    Data <- sizeComp %>%
    left_join(alk, by = c("YEAR", "REGION", "LENGTH","SEX","SPECIES_CODE")) %>%
    mutate(ageCPUE = nSizeCPUE * probability,
           AGE = ifelse(AGE > plus_group,plus_group, AGE)) %>% 
    group_by(YEAR,REGION,HAULJOIN,STRATUM,START_LONGITUDE, START_LATITUDE,nCPUE, AGE) %>%
    summarize(ageCPUE = sum(ageCPUE),
              count=n()) %>%
    ungroup() %>%
    select(HAULJOIN,AGE, ageCPUE, count) %>%
    right_join(allCats, by= c("HAULJOIN","AGE")) %>%
    mutate(ageCPUE = ifelse(is.na(ageCPUE), noAge, ageCPUE)
           ) %>%
      filter(YEAR < 2021)
    
    # Test CPUE
    # checkData <- Data %>%
    #   group_by(YEAR,HAULJOIN, nCPUE) %>%
    #   summarize(sum_age_cpue = sum(ageCPUE)) %>%
    #   mutate(diff = nCPUE - sum_age_cpue) 

    # dbSummary <- Data %>%
    #   group_by(YEAR, STRATUM, AGE) %>%
    #   summarize(meanAgeCPUE = mean(ageCPUE)) %>%
    #   inner_join(bind_rows(NBS$stratum,EBS$stratum), by="STRATUM") %>%
    #   mutate(agePopStratum = meanAgeCPUE * STRATUM_AREA) %>%
    #   group_by(YEAR, AGE) %>%
    #   summarize(agePopTotal = sum(agePopStratum)) %>%
    #   ungroup()

    # write.csv(dbSummary, "design-estimate.csv")
    
    # Format the data for VAST
    Data_Geostat <-  transmute(Data,
                             Catch_KG = ifelse(is.na(ageCPUE),0,ageCPUE),
                             Year = YEAR,
                             Vessel = "missing",
                             Age = AGE,
                             AreaSwept_km2 = 0.01, # Converts CPUE to km^2
                             Lat = START_LATITUDE,
                             Lon = START_LONGITUDE,
                             Pass = 0
    ) %>%
    data.frame()
    
    saveRDS(Data_Geostat, file = here::here("data",species_data,
                                            paste0("Data_Geostat_",species_name,".RDS")))
    
    write_csv(Data_Geostat, file = here::here("bridging", paste0("Data_Geostat_2022",species_name,".csv")))
    

# Run Analysis ------------------------------------------------------------
    # Data_Geostat <- readRDS(file = here::here("species_specific_code","BS",species_name,"data",
    #                                           paste0("Data_Geostat_",species_code,".RDS")))

  # Run model
  tic("running model")
  fit = fit_model( "settings"=settings, 
                   "Lat_i"=Data_Geostat[,'Lat'], 
                   "Lon_i"=Data_Geostat[,'Lon'], 
                   "t_i"=Data_Geostat[,'Year'],  # "t_i"=rep(2019,nrow(Data_Geostat)),
                   "c_i"=Data_Geostat[,'Age'], 
                   "b_i"=Data_Geostat[,'Catch_KG'], 
                   "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                   "v_i"=Data_Geostat[,'Vessel'],
                   Npool = Npool, 
                   test_fit=FALSE, 
                   create_strata_per_region=TRUE,
                   "working_dir" = paste0(getwd(),"/species_specific_code/BS/",species_name,"/results/"),
                   "CompileDir" = paste0(getwd(),"/species_specific_code/BS/",species_name,"/results/") )

  toc()  
  # Save results
  #saveRDS(fit, file = "VASTfit.RDS")
  saveRDS(fit, file = here::here("results",species_name,
                                paste0(species_name,"_VASTfit.RDS")))
  
  
  # Plots -------------------------------------------------------------------
  # If you need to load a fit in a new session:
  dyn.load(dynlib("VAST_v13_1_0"))
  
  # Record package versions
  sink("session_info.txt", type = "output")
  sessionInfo()
  sink()
  
  # Plot results
  results <- plot_results( fit, zrange = c(-3,3), n_cells = 600, strata_names = strata_names, check_residuals = FALSE )
  #saveRDS(results, file = "VASTresults.RDS")
  saveRDS(results,file = here::here("data",species_name,
                                    paste0(species_name,"_results.RDS")))
  
  
  # # If residual plots don't... uh... plot...
  # plot_quantile_residuals( fit=fit ) 
  # 
  # map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
  # plot_maps( Obj=fit$tmb_list$Obj, PlotDF=map_list[["PlotDF"]] )
  
 

# Expand to proportional population numbers -------------------------------
  VASTfit <- fit
  # results = plot_results( settings=settings, fit=VASTfit, check_residuals=FALSE )

  
  ## Variable names have changed in FishStatsUtils check for consistency ?calculate_proportion
  Year_Set <- VASTfit$year_labels
  Years2Include <- which(VASTfit$year_labels != 2020)
  proportions <- calculate_proportion( TmbData=VASTfit$data_list, 
                                      Index=results$Index, 
                                      year_labels = Year_Set, 
                                      years_to_plot = Years2Include,
                                      strata_names=strata_names)
  
  prop <- data.frame(t(data.frame(proportions$Prop_ctl))) %>%
    # drop_na() %>%
    rename("age_0"=1,"age_1"=2,"age_2"=3,"age_3"=4,"age_4"=5,"age_5"=6,"age_6"=7,"age_7"=8,"age_8"=9,"age_9"=10,"age_10"=11,"age_11"=12,"age_12+"=13) %>%
    bind_cols(data.frame(Year = rep(Year_Set,3) ), 
              data.frame(Region = c(rep(strata_names[1],length(Year_Set)),
                                    rep(strata_names[2],length(Year_Set)),
                                    rep(strata_names[3],length(Year_Set)) ) 
                         )
              )
  
  write.csv(prop, here::here("results",species_name,"proportions.csv"))
  saveRDS(prop, here::here("results",species_name,paste0(species_name,"_prop.RDS")))
  