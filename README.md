# model-based-indices
code used for annual model based index production by GAP in GOA and the Bering Sea
using the spatiotemporal delta-glmm developed implemented in [VAST](https://github.com/James-Thorson-NOAA/VAST) by [Jim Thorson](https://github.com/James-Thorson-NOAA)

## Authors
Primary contact: [@coleary-noaa](https://github.com/coleary-noaa)

# Mission
Provide VAST estimates of abundance and their standard error from GAP survey data for stock assessment authors in conjunction with traditional design-based estimators.
- [TOR 2022](https://docs.google.com/document/d/1t-pIruLZ-F_iNzCysWLH8cdsM0gZpuUb/edit?usp=sharing&ouid=102897708184166605880&rtpof=true&sd=true)
- [Annual stock assessment, ESP, and ESR request form](https://docs.google.com/spreadsheets/d/18gr3owj5iAq1iCDX4wpQPUC9ldLz-YTsCBIfnkHqibo/edit?usp=sharing)
- [ESP submission tool](https://apex.psmfc.org/akfin/f?p=140:LOGIN_DESKTOP:4779711459935:::::)

# Work Plan 
- Species requests & alternate model-based runs/settings in from SSMA, ESP, ESR by 01 March
- Run VAST with previous year survey data included as requested by SA author( March - April, code frozen by May )
- Merge new frozen code (index & data retrival code) in git repo
- Upload all VAST output files for hindcast and data (as an .Rdat) used to produce the indices to [google drive folder] (https://drive.google.com/drive/folders/1yxn02yF0V1PNVw0_HpqeSK_XxgOy_LAT)
- Notify the assessment lead & Cecilia (GOA) /Jason (Bering) that you’ve completed hindcasts ( March - April )
- Run hindcasts again after survey with updated data (August - Sept., completed by 30 September) and upload results to the appropriate region production folder [on google drive](https://drive.google.com/drive/folders/1yxn02yF0V1PNVw0_HpqeSK_XxgOy_LAT), including data used
- merge final index code & data pull code onto repo

# Data Reminders
- fill in zeros for tows where none of that species was observed
- For catch input you can either: (a) input as raw weight (‘Catch_KG’) with ‘AreaSwept_km2’ remaining as area swept values OR (b) Input as wCPUE (‘Catch_KG’ divided by effort)   with ‘AreaSwept_km2’ set to 1 
- For GOA: exclude areas > 700 m from the extrapolation grid area (but leave observations at this depth in the data used to fit the model) 
  to do this, follow the template code that Cecilia sent out by setting Region = "User" in make_settings(), load your extrapolation grid  'GOAThorsonGrid_Less700m' file as a     csv named ‘input_grid’ with ‘Lat’, ‘Lon’, and ‘Area_km2’ column headers, and include "input_grid"=input_grid in fit_model()
- For GOA: years to Include: 1984* - 2019 & 1984 - 2021 *if there are ID issues for your species, make note of them in your documentation/discuss with SA author.

### **Note**: if a request to only include fish west of 140 degrees, just note that there are two components to this request:  (1) excluding data from east of 170 longitude and (2) specifying this as a "strata" boundary in the model settings in the VAST code

# VAST Resources 
- VASTGAP/ModSquad folder [here](https://drive.google.com/drive/folders/1yxn02yF0V1PNVw0_HpqeSK_XxgOy_LAT)
- VAST example from Jason [here](https://drive.google.com/file/d/1GupAajXozp6afnlO3a_8I6sHC0Ev0R-b/view)
- VAST example from Cecilia [here](https://drive.google.com/file/d/1eNUXhVuezqWYQx0GHoKcKuHyqTssY_BC/view)
- VAST common troubleshooting [here](https://docs.google.com/document/d/1j3Li2aacvy7d4FJxLlGDctHlDJWQgZzI5HzyMujwf8Y/edit?usp=sharing)
- Slack VAST support group (ask Cecilia for an invite)

# Software versions used for model-based index production
- **2022**: Rv4.0.2 or later VAST v3.8.2, FishStatsUtils v2.10.0, cpp VAST_v13_1_0, TMB v1.7.22, Matrix v1.4-0, DHARMa 0.4.5
- **2021**: Rv4.0.2 VAST v3.6.1, FishStatsUtils v2.8.0, cpp VAST_v12_0_0, TMB v1.7.18, Matrix v1.2.18
- **2020**: VAST v3.3.0, FishStatsUtils v2.5.0, cpp VAST_v8_2_0

# Previous TORs
- [TOR 2021](https://docs.google.com/document/d/19gFkuNcJ_ezXzKqqOS1k5YnXyj3Tm_LyTMWGWyhy8ec/edit?usp=sharing)

# GOA specific links & files
- Link to GOA notes document for 2021 [here](https://docs.google.com/document/d/1fWEA8jftM7IRRwnCMtSjKqGRhM2Vzq7DgCCebPDc3ic/edit?usp=sharing)
- Link to GOA notes document for 2020 [here](https://docs.google.com/document/d/1M6SnI6bN16kZCuFu0Crl2BqpGB8D_CZshlqvW9TFCGY/edit?usp=sharing)

# Bering specific links & files

# R support
- Installation Guide for R 4.0.02 [here](https://docs.google.com/document/d/1tjAjvVsYbRBYLWwVdQ-Bs7GYALiUn2xFcUgcP8mQCHw/edit?usp=sharing)
