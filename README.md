# VAST_GOA
code used for model based index production in GOA

# Mission
Provide VAST estimates of abundance and their standard error from GAP survey data for stock assessment authors in conjunction with traditional design-based estimators. TOR 2021.

# Resources 
Installation Guide for R 4.0.02 here
VASTGAP/ModSquad folder here
VAST GOA template code here
VAST example from Jason here
VAST example from Cecilia here
VAST common troubleshooting here
Link to notes document for 2020 here

# Work Plan 
- Run VAST with previous year survey data included ( March - April, code frozen by May )
- Create species folder in the appropriate region results folder here & upload all VAST output files for hindcast. Per Jason’s request, please also upload your VAST R code and     data as an .Rdat that you used to produce the indices.
- Notify the assessment lead & Cecilia that you’ve completed hindcasts ( March - April )
- Fill in the details of the settings on the project issue and provide any notes ( April )
- Re-run test runs on hindcasts up to 2019 data as requested by SA author & upload to same google drive (April 2021)
- Freeze code & versions after agreement reached ( May 01 )
- Run hindcasts again after survey with updated data (August - Sept., completed by 30 September) and upload results to the appropriate region production folder on google drive, including data and code used
- Plan Team Presentations by SA authors ( September )

# Data Reminders
- fill in zeros for tows where none of that species was observed
- For catch input you can either: (a) input as raw weight (‘Catch_KG’) with ‘AreaSwept_km2’ remaining as area swept values OR (b) Input as wCPUE (‘Catch_KG’ divided by effort)   with ‘AreaSwept_km2’ set to 1 
- Exclude areas > 700 m from the extrapolation grid area (but leave observations at this depth in the data used to fit the model) 
  to do this, follow the template code that Cecilia sent out by setting Region = "User" in make_settings(), load your extrapolation grid  'GOAThorsonGrid_Less700m' file as a     csv named ‘input_grid’ with ‘Lat’, ‘Lon’, and ‘Area_km2’ column headers, and include "input_grid"=input_grid in fit_model()
- Years to Include: 1984* - 2019 & 1984 - 2021 *if there are ID issues for your species, make note of them in your documentation. When we move to AI, we will discuss what to     do about data before 1993
- Combine catches of dusky rockfish and  "dusky and dark rockfishes unid."

## Note: if a request to only include fish west of 140 degrees, just note that there are two components to this request:  (1) excluding data from east of 170 longitude and (2) specifying this as a "strata" boundary in the model settings in the VAST code
###
#### 2020 settings: VAST v3.3.0, FishStatsUtils v2.5.0, cpp VAST_v8_2_0
#### 2021 settings: Rv4.0.2 VAST v3.6.1, FishStatsUtils v2.8.0, cpp VAST_v12_0_0, TMB v1.7.18, Matrix v1.2.18
