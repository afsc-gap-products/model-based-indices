# Get data from Google drive ----------------------------------------------------------------

drive_auth()
1

# Summary catch data for EBS and NBS with NBS 2018 included
googledrive::drive_download(file=as_id("1TctmzLjuFUopvdD9jqBnhNxw1NuAi87c"),
                            path=here::here("data","EBS_NBS_Index.RDS"),
                            overwrite = TRUE)

# Cold Pool Area covariate
googledrive::drive_download(file=as_id("119pehim03WxFjGH9tzjUF7gZDcHJeUe8"),
                            path=here::here("data","cpa_areas2019.csv"),
                            overwrite = TRUE)

# EBS Strata
googledrive::drive_download(file=as_id("1UZ4OBGuwqSpPhTD3idDUNqO57Fm1GYnL"),
                            path=here::here("data", "shapefiles","EBS_NBS_2019_1983.cpg"),
                            overwrite = TRUE)
googledrive::drive_download(file=as_id("1GhH47aoQ42kx3TYqMo_w0zTV4UeXYWsN"),
                            path=here::here("data", "shapefiles","EBS_NBS_2019_1983.dbf"),
                            overwrite = TRUE)
googledrive::drive_download(file=as_id("1SLOH6Ggp8PufL0XZPXLa8ZzPrwIgz68S"),
                            path=here::here("data", "shapefiles","EBS_NBS_2019_1983.sbn"),
                            overwrite = TRUE)
googledrive::drive_download(file=as_id("1Hk_t3RvwqwHp6ypL4yfUsACXfxq9IynZ"),
                            path=here::here("data", "shapefiles","EBS_NBS_2019_1983.prj"),
                            overwrite = TRUE)
googledrive::drive_download(file=as_id("16GPjJfiWL5ZNCHdimKvoQVp63PGvrVmn"),
                            path=here::here("data", "shapefiles","EBS_NBS_2019_1983.sbx"),
                            overwrite = TRUE)
googledrive::drive_download(file=as_id("1Pfg-HxarbSU2M0u-HzSHAIj7E8v-MeO4"),
                            path=here::here("data", "shapefiles","EBS_NBS_2019_1983.shp"),
                            overwrite = TRUE)
googledrive::drive_download(file=as_id("1Ke_9cy5wwXolzx34TBM1gbRPGVaS2g7b"),
                            path=here::here("data", "shapefiles","EBS_NBS_2019_1983.shp.xml"),
                            overwrite = TRUE)
googledrive::drive_download(file=as_id("1wqdoTKjVSziRdQnUQa7COuYU_2H_TyAt"),
                            path=here::here("data", "shapefiles","EBS_NBS_2019_1983.shx"),overwrite = TRUE)


# Format catch data -------------------------------------------------------

sumAll <- read_rds(here::here("data","EBS_NBS_Index.RDS"))
# test = filter(sumAll, HAULJOIN %in% c(-13853,-13852,-13851,-13850,-13849,-13848,-13847,-13846))

Data <- sumAll %>%
  dplyr::filter(SPECIES_CODE==species)

# Format the data for VAST
Data_Geostat <-  dplyr::transmute(Data,
                                  Catch_KG = wCPUE*100, # sumfish calculates CPUE in kg/ha this converts it to kg/km^2
                                  Year = YEAR,
                                  Vessel = "missing",
                                  AreaSwept_km2 = 1, # area swept is 1 when using CPUE instead of observed weight
                                  Lat = START_LATITUDE,
                                  Lon = START_LONGITUDE,
                                  Pass = 0
) %>%
  data.frame()

saveRDS(Data_Geostat, file = "Data_Geostat.RDS")
