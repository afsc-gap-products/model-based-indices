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
target_file <- drive_get(as_id("1NI3OVgFgT932LUnZHVgMusYt8yn-V3d-"))  

# Set path for where to download the file. Create folder if it doesn't exist
folder_path <- here(
  "species_specific_code",
  "BS",
  "pollock",
  "data"
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
drive_folder <- as_id("1U_RpXBnILwWoEWVDmi1ctSScZZ6sXqX5")  

# List local files in the results directory
# results_dir <-   # DEFINE DIRECTORY HERE
results_files <- list.files(results_dir, full.names = TRUE)

# Upload all files to google drive
walk(results_files, ~ drive_upload(
  media = .x,
  path = drive_folder,
  overwrite = TRUE # replace existing files?
))
