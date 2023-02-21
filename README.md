# model-based-indices
code used for annual model based index production by GAP in GOA and the Bering Sea. Model-based-indices are produced
using the spatiotemporal delta-glmm implemented in [VAST](https://github.com/James-Thorson-NOAA/VAST), software developed by [Jim Thorson](https://github.com/James-Thorson-NOAA)

## Authors
Primary contact: [@coleary-noaa](https://github.com/coleary-noaa)

# Mission
Provide VAST estimates of abundance and their standard error from GAP survey data for stock assessment authors in conjunction with traditional design-based estimators.
- [TOR 2023](https://drive.google.com/file/d/1Jog0lZx3Y9D0yRqDSl2X-BSqR_zyEJUw/view)
- [Annual stock assessment request form](https://docs.google.com/spreadsheets/d/18gr3owj5iAq1iCDX4wpQPUC9ldLz-YTsCBIfnkHqibo/edit?usp=sharing)
- [Annual ESR request form](https://docs.google.com/spreadsheets/d/1SC-KzRmng0c2e1GpvqYQz4Ijr0_hlbGZbWiislXGGcw/edit?usp=sharing)
- [ESP submission tool](https://apex.psmfc.org/akfin/f?p=140:LOGIN_DESKTOP:4779711459935:::::)
- [2023 GOA notes](https://docs.google.com/document/d/1u97oBL31v5llxO1-5GLu1-x7k0zbU0LgviFEkWYdSCs/edit)

# GAP Work Plan 
- Species requests & alternate model-based runs/settings in from SSMA, ESP, ESR by 15 January
- Run VAST with previous year survey data included as requested by SA author( March - April, code frozen by May )
- Merge new frozen code (index & data retrival code) in git repo
- Upload all VAST output files for hindcast and data (as an .Rdat) used to produce the indices to [google drive folder](https://drive.google.com/drive/folders/1yxn02yF0V1PNVw0_HpqeSK_XxgOy_LAT)
- Notify the assessment lead & Cecilia (GOA) /Lewis (Bering) that you’ve completed hindcasts ( March - April )
- Run hindcasts again after survey with updated data (August - Sept., completed by 30 September) and upload results to the appropriate region production folder [on google drive](https://drive.google.com/drive/folders/1yxn02yF0V1PNVw0_HpqeSK_XxgOy_LAT), including data used
- merge final index code & data pull code onto repo

# General Production Timeline Estimate
**March - April is the time window for any model exploration/iteration, as requested**
  - December - ModSquad leadership meets to prioritize species + research plans
  - January 15 - Deadline for SAFE authors to submit product requests
  - late January - ModSquad planning meeting with program leads
  - February 15 - SAFE author's revised requests given ModSquad capacity and discussion from REFM meeting
  - March 15 - TOR memo finalized and posted on the web
  - April 1 - Crab hindcasts completed
  - May 1 - Groundfish hindcasts completed (stop all research, decide on which to not recommend for use in base model)
  - August 25 - Model-based estimates for EBS crabs completion target date
  - September 25 - Model-based estimates for NBS crabs target date
  - September 30 - Model-based estimates for groundfishes completion target date, pending bottom trawl data QAQC timeline
  - October 15 - Final deadline for completing groundfish model-based estimates 

# VAST Resources 
- VASTGAP/ModSquad folder [here](https://drive.google.com/drive/folders/1yxn02yF0V1PNVw0_HpqeSK_XxgOy_LAT)
- ModSquad training materials [here](https://drive.google.com/drive/folders/1TZRrEwka7OEICC6D81LiuGtYJMSxdaFe?usp=sharing)
- VAST example from Cecilia [here](https://drive.google.com/file/d/1eNUXhVuezqWYQx0GHoKcKuHyqTssY_BC/view)
- VAST common troubleshooting [here](https://docs.google.com/document/d/1j3Li2aacvy7d4FJxLlGDctHlDJWQgZzI5HzyMujwf8Y/edit?usp=sharing)
- Identification of valid hauls for Bering Sea VAST indices, from Jason [here](https://docs.google.com/spreadsheets/d/1-z7AFYoTM0-RApsW9APfXX4CoA5mbk6k/edit#gid=1419989689)
- Slack VAST support group (ask Cecilia for an invite)

# Minimum standard of versions of software and key dependencies (i.e., using these versions or later)
- **2023**: 
  - Rv4.0.2: VAST v3.9.0, FishStatsUtils v2.10.0, cpp VAST_v13_1_0, TMB v1.7.22, Matrix v1.4-0, DHARMa 0.4.5; or,
  - MRAN v4.0.2: VAST v3.9.0, FishStatsUtils v2.10.0, cpp VAST_v13_1_0, TMB v1.7.16, Matrix v1.2-18, DHARMa v0.3.2
- **2022**: 
  - Rv4.0.2: VAST v3.8.2, FishStatsUtils v2.10.0, cpp VAST_v13_1_0, TMB v1.7.22, Matrix v1.4-0, DHARMa 0.4.5; or,
  - MRAN v4.0.2: VAST v3.9.0, FishStatsUtils v2.10.0, cpp VAST_v13_1_0, TMB v1.7.16, Matrix v1.2-18, DHARMa v0.3.2
- **2021**: Rv4.0.2: VAST v3.6.1, FishStatsUtils v2.8.0, cpp VAST_v12_0_0, TMB v1.7.18, Matrix v1.2.18
- **2020**: VAST v3.3.0, FishStatsUtils v2.5.0, cpp VAST_v8_2_0

# Initial model settings to generate model-based indices
| Initial Model Setting  | Suggested alternative setting (if needed) |
| :---         | :--- |
| purpose = "index2” in make_settings()  | NA  |
| knots = 750 in make_settings()  | knots = 500, 1000  |
| Poisson-link delta-gamma observation model[^1]: <br/> ObsModel = c(2,1) in make_settings()  | option 2: Tweedie ObsModel = c(10,2)[^2] <br/> option 3: delta-lognormal ObsModel = c(1,1)  |
| knot_method = ‘grid’ in fit_model()  | knot_method = ‘samples’ if necessary to aid convergence or for comparison to a previous model fit  |
| fine_scale = TRUE in make_settings()  | NA  |
| bias.correct = TRUE in make_settings()  | NA  |
| refine = TRUE in fit_model()  | refine = FALSE  |
| spatiotemporal fields: “IID” default settings for FieldConfig in make_settings()  | model spatiotemporal components (epsilon) as a first-order autoregressive process “AR1” (required for extremely unbalanced data) or “0” (if necessary to aid convergence)  |
| anisotropy is on (use_anisotropy = TRUE) in make_settings()  | anisotropy off (use_anisotropy = FALSE) if necessary to aid convergence  |
| no vessel effects, catchability or density covariates in fit_model()  | may include covariates in cases where their incorporation has been previously demonstrated to improve model fit (e.g., a spatially varying response to cold-pool extent when generating abundance indices combining the EBS and NBS); in these cases, covariates will be centered and scaled prior to fitting  |
[^1]: Noting that for species with 100% encounters in any year we will use c(2,4) instead of c(2,1), or the equivalent setting for the lognormal <br/>
[^2]: Tweedie also involves additional changes to RhoConfig and FieldConfig to ensure that there is only a single linear predictor being estimated, as documented elsewhere

# Previous TORs
- [TOR 2022](https://drive.google.com/file/d/1t13fVai8HOo3xtwxn_vazo9NDT8MrVar/view?usp=sharing)
- [TOR 2021](https://docs.google.com/document/d/19gFkuNcJ_ezXzKqqOS1k5YnXyj3Tm_LyTMWGWyhy8ec/edit?usp=sharing)

# Data Reminders
- fill in zeros for tows where none of that species was observed
- For catch input you can either: (a) input as raw weight (‘Catch_KG’) with ‘AreaSwept_km2’ remaining as area swept values OR (b) Input as wCPUE (‘Catch_KG’ divided by   effort) with ‘AreaSwept_km2’ set to 1 
- For GOA: exclude areas > 700 m from the extrapolation grid area (but leave observations at this depth in the data used to fit the model) 
  to do this, follow the template code that Cecilia sent out by setting Region = "User" in make_settings(), load your extrapolation grid  'GOAThorsonGrid_Less700m'       file as a csv named ‘input_grid’ with ‘Lat’, ‘Lon’, and ‘Area_km2’ column headers, and include "input_grid"=input_grid in fit_model()

### **Note**: if a request to only include fish west of 140 degrees, just note that there are two components to this request:  (1) excluding data from east of 170 longitude and (2) specifying this as a "strata" boundary in the model settings in the VAST code

# GOA specific links & files
- Link to GOA notes document for 2022 [here](https://docs.google.com/document/d/1TOMz9zaRaRd2xLsaE8J9KLn1wPYqtAwgAprJrhsLtSE/edit?usp=sharing)
- Link to GOA notes document for 2021 [here](https://docs.google.com/document/d/1fWEA8jftM7IRRwnCMtSjKqGRhM2Vzq7DgCCebPDc3ic/edit?usp=sharing)
- Link to GOA notes document for 2020 [here](https://docs.google.com/document/d/1M6SnI6bN16kZCuFu0Crl2BqpGB8D_CZshlqvW9TFCGY/edit?usp=sharing)

# Bering specific links & files

# R support
- Installation Guide for R 4.0.02 [here](https://docs.google.com/document/d/1tjAjvVsYbRBYLWwVdQ-Bs7GYALiUn2xFcUgcP8mQCHw/edit?usp=sharing)
