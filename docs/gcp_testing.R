# Code for testing packages / connections needed for Modsquad activities on GCP

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

con <- DBI::dbConnect(
  odbc::odbc(),
  .connection_string = paste0("Driver=/opt/oracle/instantclient_12_2/libsqora.so.12.1;DBQ=raja.afsc.noaa.gov:1521/afsc;UID=", 
                              oracle_user, ";PWD=", oracle_pw)
)

# Connect to Google Drive -----------------------------------------------------
# googledrive::drive_auth(path="/etc/sa_key.json")  # to connect to the default google drive account associated with the instance
library(gargle)
library(googledrive)

# Connect to google drive using your (probs NOAA) email
gdrive_email <- rstudioapi::showPrompt(title = "Email",
                                       message = "Email for Google Drive",
                                       default = "")

drive_auth(token = credentials_user_oauth2(
  scopes = "https://www.googleapis.com/auth/drive", 
  email = gdrive_email))

drive_user()  # check user account

# # download test
# googledrive::drive_download(
#   file = googledrive::as_id(
#   x = "https://drive.google.com/file/d/1Y1fbaesHacPoqgloD350XJq2VF2Jy51-/view?usp=drive_link"))
# 
# # upload test
# write.csv(x = "hello world", file = "./ex.txt")
# googledrive::drive_upload(media = "./ex.txt")

# Test sdmTMB, from their website: --------------------------------------------
# https://pbs-assess.github.io/sdmTMB/index.html#basic-use
library(sdmTMB)

mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)

fit <- sdmTMB(
  density ~ s(depth),
  data = pcod,
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on"
)

fit  # print fit - if this prints, everything should be working!

# Test tinyVAST - from this vignette: -----------------------------------------
# https://vast-lib.github.io/tinyVAST/articles/web_only/age_composition_expansion.html
library(tinyVAST)
library(fmesher)
library(sf)

format_data <- function() {
  # Pull & format data
  data(bering_sea_pollock_ages)
  Data <- subset(bering_sea_pollock_ages, Year >= 2021)
  Data$Age <- factor(paste0("Age_",Data$Age))
  Data$Year_Age <- interaction(Data$Year, Data$Age)
  
  # Project to UTM
  Data <- st_as_sf(Data, 
                   coords = c('Lon','Lat'),
                   crs = st_crs(4326))
  Data <- st_transform(Data, crs = st_crs("+proj=utm +zone=2 +units=km"))
  Data <- cbind(st_drop_geometry(Data), st_coordinates(Data))
  
  return(Data)
}

Data <- format_data()

# Set up tinyVAST settings
sem <- ""

dsem <- "
  Age_1 -> Age_1, 1, lag1
  Age_2 -> Age_2, 1, lag1
  Age_3 -> Age_3, 1, lag1
  Age_4 -> Age_4, 1, lag1
  Age_5 -> Age_5, 1, lag1
  Age_6 -> Age_6, 1, lag1
  Age_7 -> Age_7, 1, lag1
  Age_8 -> Age_8, 1, lag1
  Age_9 -> Age_9, 1, lag1
  Age_10 -> Age_10, 1, lag1
  Age_11 -> Age_11, 1, lag1
  Age_12 -> Age_12, 1, lag1
  Age_13 -> Age_13, 1, lag1
  Age_14 -> Age_14, 1, lag1
  Age_15 -> Age_15, 1, lag1
"

mesh <- fm_mesh_2d(loc = Data[,c("X","Y")],
                   cutoff = 50)
control <- tinyVASTcontrol(getsd = FALSE,
                           profile = c("alpha_j"),  
                           trace = 0)
family <- list(
  Age_1 = tweedie(),
  Age_2 = tweedie(),
  Age_3 = tweedie(),
  Age_4 = tweedie(),
  Age_5 = tweedie(), 
  Age_6 = tweedie(),
  Age_7 = tweedie(),
  Age_8 = tweedie(),
  Age_9 = tweedie(),
  Age_10 = tweedie(),
  Age_11 = tweedie(),
  Age_12 = tweedie(),
  Age_13 = tweedie(),
  Age_14 = tweedie(),
  Age_15 = tweedie()     
)

# Fit tinyVAST model
myfit <- tinyVAST(
  data = Data,
  formula = Abundance_per_hectare ~ 0 + Year_Age,
  sem = sem,
  dsem = dsem,
  family = family,
  space_column = c("X", "Y"), 
  variable_column = "Age",
  time_column = "Year",
  distribution_column = "Age",
  spatial_graph = mesh,
  control = control
)

myfit  # Print fit - if this works everything is groovy!