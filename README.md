# VAST_GOA
code used for model based index production in GOA

# Mission
Provide VAST estimates of abundance and their standard error from GAP survey data for stock assessment authors in conjunction with traditional design-based estimators.
- [TOR 2021.](https://docs.google.com/document/d/19gFkuNcJ_ezXzKqqOS1k5YnXyj3Tm_LyTMWGWyhy8ec/edit?usp=sharing)
- [Annual assessment and ESP/ESR reqest form](https://docs.google.com/spreadsheets/d/18gr3owj5iAq1iCDX4wpQPUC9ldLz-YTsCBIfnkHqibo/edit?usp=sharing)

# Work Plan 
- Run VAST with previous year survey data included ( March - April, code frozen by May )
- Upload all VAST output files for hindcast,R code and data as an .Rdat that you used to produce the indices
- Notify the assessment lead & Cecilia that you’ve completed hindcasts ( March - April )
- Fill in the details of the settings on the project issue and provide any notes ( April )
- Re-run test runs on hindcasts up to previous year's data as requested by SA author & upload to same google drive (April )
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

### **Note**: if a request to only include fish west of 140 degrees, just note that there are two components to this request:  (1) excluding data from east of 170 longitude and (2) specifying this as a "strata" boundary in the model settings in the VAST code

# Resources 
- Installation Guide for R 4.0.02 [here](https://docs.google.com/document/d/1tjAjvVsYbRBYLWwVdQ-Bs7GYALiUn2xFcUgcP8mQCHw/edit?usp=sharing)
- VASTGAP/ModSquad folder [here](https://drive.google.com/drive/folders/1yxn02yF0V1PNVw0_HpqeSK_XxgOy_LAT)
- VAST GOA template code [here](https://drive.google.com/drive/folders/1eV5CfsVH7b2UzVEuTavHfRQGJnHcqgNy)
- VAST example from Jason [here](https://drive.google.com/file/d/1GupAajXozp6afnlO3a_8I6sHC0Ev0R-b/view)
- VAST example from Cecilia [here](https://drive.google.com/file/d/1eNUXhVuezqWYQx0GHoKcKuHyqTssY_BC/view)
- VAST common troubleshooting [here](https://docs.google.com/document/d/1j3Li2aacvy7d4FJxLlGDctHlDJWQgZzI5HzyMujwf8Y/edit?usp=sharing)
- Link to notes document for 2021 [here](https://docs.google.com/document/d/1fWEA8jftM7IRRwnCMtSjKqGRhM2Vzq7DgCCebPDc3ic/edit?usp=sharing)
- Link to notes document for 2020 [here](https://docs.google.com/document/d/1M6SnI6bN16kZCuFu0Crl2BqpGB8D_CZshlqvW9TFCGY/edit?usp=sharing)
- Slack VAST support group (ask me for an invite)

# Software versions
- **2021**: Rv4.0.2 VAST v3.6.1, FishStatsUtils v2.8.0, cpp VAST_v12_0_0, TMB v1.7.18, Matrix v1.2.18
- **2020**: VAST v3.3.0, FishStatsUtils v2.5.0, cpp VAST_v8_2_0
