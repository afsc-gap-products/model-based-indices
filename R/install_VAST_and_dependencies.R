## install VAST and dependencies with versions as in 2023 TOR

.rs.restartR() # restart R before doing these updates - these are delicate installs
library(devtools)

install_version("Matrix", version = "1.5-3", repos = "http://cran.us.r-project.org")
install_version("TMB", version = "1.9.2", repos = "http://cran.us.r-project.org")
install_version("DHARMa", version = "0.4.6", repos = "http://cran.us.r-project.org")
install_github("james-thorson/FishStatsUtils@2.12.0")
install_github("James-Thorson-NOAA/VAST@3.10.0") 
#devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

# or try this if you are having trouble with the VAST installation
#devtools::install_github("james-thorson/VAST@3.10.0", INSTALL_opts="--no-staged-install")

#generate in github and load here token when API limit exceeded
#Sys.setenv(GITHUB_PAT = 'long_key_generated_by_github')

##installs from local files if your API limit is exceeded
##download the .tar.gz to your local downloads from the CRAN or github webpage
## create path to downloadeds
#wd <- getwd()
#wd_split <- strsplit(wd, "\\/")
#usr <- paste(wd_split[[1]][1:3], collapse = "/")

##examples of downloading from your local computer downloads using .tar.gz
#devtools::install_local(paste0(usr,"/Downloads/FishStatsUtils-2.12.0.tar.gz"), dependencies=FALSE) 
#devtools::install_local(paste0(usr,"/Downloads/VAST-3.10.0.tar.gz"), dependencies=FALSE)
#devtools::install_local(paste0(usr,"/Downloads/Matrix_1.4-0.tar.gz"), dependencies=FALSE)
#devtools::install_local(paste0(usr,"/Downloads/DHARMa_0.4.5.tar.gz"), dependencies=FALSE)


# Check R packages ------------------------------------------------------------------------------------

current_year <- 2023
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
