## install VAST and dependencies

#generate in github and load here token when API limit exceeded
#Sys.setenv(GITHUB_PAT = 'long_key_generated_by_github')

##installs from github examples when trying to install a specific version of the package
#devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
#devtools::install_github("james-thorson/FishStatsUtils@2.11.0")
#devtools::install_github("james-thorson/VAST@3.9.0", INSTALL_opts="--no-staged-install")

##installs from local files if your API limit is exceeded
##download the .tar.gz to your local downloads from the CRAN or github webpage
## create path to downloadeds
wd <- getwd()
wd_split <- strsplit(wd, "\\/")
usr <- paste(wd_split[[1]][1:3], collapse = "/")

##examples of downloading from your local computer downloads using .tar.gz
devtools::install_local(paste0(usr,"/Downloads/FishStatsUtils-2.12.0.tar.gz"), dependencies=FALSE) 
devtools::install_local(paste0(usr,"/Downloads/VAST-3.10.0.tar.gz"), dependencies=FALSE)
#devtools::install_local(paste0(usr,"/Downloads/Matrix_1.4-0.tar.gz"), dependencies=FALSE)
#devtools::install_local(paste0(usr,"/Downloads/DHARMa_0.4.5.tar.gz"), dependencies=FALSE)
