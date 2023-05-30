library(dplyr)
library(VAST)
library(tictoc)

# Set species, model -------------------------------------------------------

which_model <- c("hindcast", "production")[1]
compare <- FALSE # If compare = TRUE, using prior year's alk
species <- 21720
species_name <- "pacific_cod"

workDir <- paste0(getwd(),"/species_specific_code/BS/", 
                  species_name, "/", which_model, "/")
if(!dir.exists(workDir))
  dir.create(path = workDir, recursive = TRUE)
if(!dir.exists(paste0(workDir, "results_age/")))
  dir.create(path = paste0(workDir, "results_age/"), recursive = TRUE)

# Record sessionInfo -------------------------------------------------------
sink(file = paste0(workDir, "results_age/session_info.txt"), 
     type = "output")
sessionInfo()
sink()

# Make sure package versions are correct for current year ------------------
current_year <- 2023
prev_year <- current_year-1
VAST_cpp_version <- "VAST_v14_0_1"
pck_version <- c("VAST" = "3.10.0",
                 "FishStatsUtils" = "2.12.0",
                 "Matrix" = "1.5-3",
                 "TMB" = "1.9.2",
                 "DHARMa" = "0.4.6")

for (pck in 1:length(pck_version)) {
  temp_version <- packageVersion(pkg = names(pck_version)[pck])
  
  if(temp_version == pck_version[pck])
    message(paste0("The version of the '", names(pck_version)[pck], 
                   "' package (", temp_version, ") is consistent",
                   " with the ", current_year, " TOR."))
  
  if(!temp_version == pck_version[pck])
    message(paste0("WARNING: ", 
                   "The version of the '", names(pck_version)[pck], 
                   "' package (", temp_version, ") is NOT consistent",
                   " with the ", current_year, " TOR. Please update the '", 
                   names(pck_version)[pck], "' package to ", 
                   pck_version[pck]))
  
  rm(pck, temp_version)
}

# Settings ----------------------------------------------------------------
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
                      Version = VAST_cpp_version,
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
    alk <- readRDS(paste0(workDir,"data/unstratified_alk.RDS"))
    if(compare == TRUE)
    {
      alk_all <- readRDS(paste0(workDir,"data/",prev_year,"_production/unstratified_alk.RDS") )
      alk_ebs <- alk_all %>% filter(REGION == "EBS")
      alk_nbs <- alk_all %>% filter(REGION == "NBS") %>%
          bind_rows( filter(alk_ebs, YEAR == 2018) )   # Use EBS ALK for 2018 ad hoc sampling in NBS

      alk <- bind_rows(alk_ebs, alk_nbs)
    }
    
    sizeComp <- readRDS(paste0(workDir,"data/EBS_NBS_SizeComp.RDS") ) %>% 
      filter(YEAR >= min_year, !is.na(EFFORT))

    haulData <- readRDS(paste0(workDir,"data/EBS_NBS_Index.RDS") ) %>% 
        dplyr::filter(YEAR >= min_year, !is.na(EFFORT))
    
    # Get summary calculations - this also fills zeroes within year - add sex if generating sex-specific agecomps
    allCats <- expand.grid(HAULJOIN=unique(haulData$HAULJOIN), AGE = unique(alk$AGE[alk$AGE<=plus_group]), noAge = 0) %>%
    inner_join(haulData, by = c("HAULJOIN")) 
    
    
    # Aggregate by Age key
    Data <- sizeComp %>%
    left_join(alk, by = c("YEAR", "REGION", "LENGTH","SEX","SPECIES_CODE"), relationship = "many-to-many") %>%
    mutate(ageCPUE = nSizeCPUE * probability,
           AGE = ifelse(AGE > plus_group,plus_group, AGE)) %>% 
    group_by(YEAR,REGION,HAULJOIN,STRATUM,START_LONGITUDE, START_LATITUDE,nCPUE, AGE) %>%
    summarize(ageCPUE = sum(ageCPUE),
              count=n()) %>%
    ungroup() %>%
    select(HAULJOIN,AGE, ageCPUE, count) %>%
    right_join(allCats, by= c("HAULJOIN","AGE")) %>%
    mutate(ageCPUE = ifelse(is.na(ageCPUE), noAge, ageCPUE)) # %>% filter(YEAR < 2021) 
    
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
    
    saveRDS(Data_Geostat, file = paste0(workDir,"data/Data_Geostat_comps_",species_name,".RDS"))
    

# Run Analysis ------------------------------------------------------------
    # Data_Geostat <- readRDS(file = here::here("species_specific_code","BS",species_name,which_model,"data",
    #                                           paste0("Data_Geostat_comps_",species_code,".RDS")))

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
                   test_fit=TRUE, # set to FALSE if want to avoid interruption due to marginal fit issues
                   create_strata_per_region=TRUE,
                   "working_dir" = paste0(workDir,"/results_age/"),
                   "CompileDir" = paste0(workDir,"/results_age/")
  )

  toc()  
  # Save results
  #saveRDS(fit, file = "VASTfit.RDS")
  saveRDS(fit, file = paste0(workDir,"/results_age/",species_name,"_VASTfit.RDS"))
  
  
  # Plots -------------------------------------------------------------------
  # If you need to load a fit in a new session:
  #dyn.load(dynlib("VAST_v14_0_1"))
  
  # Plot results
  results <- plot_results( fit, zrange = c(-3,3), n_cells = 2000, strata_names = strata_names, 
                           check_residuals = TRUE,  working_dir =  paste0(workDir,"/results_age/"))
  #saveRDS(results, file = "VASTresults.RDS")
  saveRDS(results,file = paste0(workDir,"/results_age/",species_name,"_results.RDS"))
  
  
  # # If residual plots don't... uh... plot...
  # plot_quantile_residuals( fit=fit)
  # 
  # map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )
  # plot_maps(fit=fit, n_cells = 2000, Obj=fit$tmb_list$Obj, PlotDF=map_list[["PlotDF"]])
  
 

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
  
  write.csv(prop, paste0(workDir,"/results_age/",species_name,"_proportions.csv"))
  saveRDS(prop, paste0(workDir,"/results_age/",species_name,"_proportions.RDS"))
  