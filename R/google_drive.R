#' Script for: 
#' (1) downloading files from google drive; 
#' (2) uploading the contents of a folder to google drive.
#' The download function is largely for pollock, if the data were prepared
#' elsewhere. The upload function can be used for any data or model products 
#' created on GCP (or locally). 
#' 
#' IMPORTANT: for uploading, you must create the upload folder in google drive
#' before copying in the ID string for the folder.

library(gargle)
library(googledrive)
library(purrr)

# Authorize and connect to google drive using NOAA email ----------------------
gdrive_email <- rstudioapi::showPrompt(title = "Email",
                                       message = "Email for Google Drive",
                                       default = "")

drive_auth(token = credentials_user_oauth2(
  scopes = "https://www.googleapis.com/auth/drive", 
  email = gdrive_email))

drive_user()  # check user account

# Download file from google drive ---------------------------------------------
# Use the string to avoid problems with duplicate file names
# Tip: copy in the share link and then remove everything but the long string
target_file <- drive_get(as_id("1NGKLq1__ZpnUD7qZpqqEUS5uxp5KCKle"))  

# Set path for where to download the file. Create folder if it doesn't exist
folder_path <- here(
  "species_specific_code",
  "BS",
  "pollock",
  "hindcast",
  "results"
)
if(!dir.exists(folder_path)) {
  dir.create(folder_path)
}

# Download file
drive_download(
  file = target_file,
  path = here(
    folder_path, 
    target_file$name
  ),
  overwrite = TRUE
)

# Upload contents of a folder to google drive ---------------------------------
# Access drive folder via the string at the end of the URL (click into it in google drive)
drive_folder <- as_id("1NtUfZilNZH0homjvbEDO8YK2YajsMNGg")  

# List local files in the results directory
# DEFINE DIRECTORY HERE
folder_path <- here(
  "species_specific_code",
  "BS",
  "pollock",
  "production",
  "results"
)

results_files <- list.files(folder_path, full.names = TRUE)

# Upload all files to google drive
walk(results_files, ~ drive_upload(
  media = .x,
  path = drive_folder,
  overwrite = TRUE # replace existing files?
))
