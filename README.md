# model-based-indices
Code used for annual model based index production by the AFSC GAP ModSquad in GOA and the Bering Sea. Model-based-indices are produced
using the spatiotemporal delta-glmm implemented in [VAST](https://github.com/James-Thorson-NOAA/VAST), software developed by [Jim Thorson](https://github.com/James-Thorson-NOAA)

## Authors
Primary contact: [@Lewis-Barnett-NOAA](https://github.com/Lewis-Barnett-NOAA)

# Mission
Provide VAST estimates of abundance and their standard error from GAP survey data for stock assessment authors in conjunction with traditional design-based estimators.
- [TOR or Terms Of Reference](https://drive.google.com/file/d/1TmW9gySrWVGxv5OXpTbkMml6vAbi14Xl/view?usp=drive_link)
- [Annual stock assessment request form](https://docs.google.com/spreadsheets/d/18gr3owj5iAq1iCDX4wpQPUC9ldLz-YTsCBIfnkHqibo/edit?usp=sharing)
- [Annual ESR request form](https://docs.google.com/spreadsheets/d/1SC-KzRmng0c2e1GpvqYQz4Ijr0_hlbGZbWiislXGGcw/edit?usp=sharing)
- [ESP submission tool](https://apex.psmfc.org/akfin/f?p=140:LOGIN_DESKTOP:4779711459935:::::)

# GAP Work Plan 
- Species requests & alternate model-based runs/settings in from SSMA, ESP, ESR by 15 January
- Run VAST with previous year survey data included as requested by SA author (March - April, code frozen by May)
- Merge new frozen code (index & data retrival code) in git repo
- Upload all VAST output files for hindcast and data (as an .Rdat) used to produce the indices to the appropriate region hindcast folder in the ModSquad [google drive folder](https://drive.google.com/drive/folders/1yxn02yF0V1PNVw0_HpqeSK_XxgOy_LAT)
- Notify the assessment lead & ModSquad lead (Lewis) that you’ve completed hindcasts (March - April)
- Run hindcast model structure again after survey with updated data (August - Sept., completed by 30 September) and upload results and data (as an .Rdat) used to produce the indices to the appropriate region production folder [on google drive](https://drive.google.com/drive/folders/1yxn02yF0V1PNVw0_HpqeSK_XxgOy_LAT)
- Merge final index code & data pull code onto repo
- Notify the assessment lead & ModSquad lead (Lewis) that you’ve completed the production estimates

# General Production Timeline Estimate
**March - April is the time window for any model exploration/iteration, as requested**
  - December - ModSquad leadership meets/communicates to prioritize species + research plans
  - January 15 - Deadline for SAFE authors to submit product requests
  - late January - ModSquad planning meeting with program leads as needed
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
- VAST common troubleshooting [here](https://docs.google.com/document/d/1j3Li2aacvy7d4FJxLlGDctHlDJWQgZzI5HzyMujwf8Y/edit?usp=sharing)
- Identification of valid hauls for Bering Sea VAST indices, from Jason [here](https://docs.google.com/spreadsheets/d/1-z7AFYoTM0-RApsW9APfXX4CoA5mbk6k/edit#gid=1419989689)

# Minimum standard of versions of software and key dependencies (i.e., using these versions or later)
- **2024**: 
  - Rv4.0.2: VAST v3.9.0, FishStatsUtils v2.10.0, cpp VAST_v13_1_0, TMB v1.7.22, Matrix v1.4-0, DHARMa 0.4.5; or,
  - MRAN v4.0.2: VAST v3.9.0, FishStatsUtils v2.10.0, cpp VAST_v13_1_0, TMB v1.7.16, Matrix v1.2-18, DHARMa v0.3.2
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
| knots = 750 in make_settings()  | knots = 500, 1000; 50 for age comps |
| Poisson-link delta-gamma observation model[^1]: <br/> ObsModel = c(2,1) in make_settings()  | option 2: Tweedie ObsModel = c(10,2)[^2] <br/> option 3: delta-lognormal ObsModel = c(1,1)  |
| knot_method = ‘grid’ in fit_model()  | knot_method = ‘samples’ if necessary to aid convergence or for comparison to a previous model fit  |
| fine_scale = TRUE in make_settings()  | fine_scale = FALSE  |
| max_cells = 2000 in make_settings()  | max_cells = n_x * 10  |
| Npool = 100 in make_model()  | Npool = 20 to 99  |
| bias.correct = TRUE in make_settings()  | NA  |
| refine = TRUE in fit_model()  | refine = FALSE  |
| spatiotemporal fields: “IID” default settings for FieldConfig in make_settings()  | model spatiotemporal components (epsilon) as a first-order autoregressive process “AR1” (required for extremely unbalanced data) or “0” (if necessary to aid convergence)  |
| anisotropy is on (use_anisotropy = TRUE) in make_settings()  | anisotropy off (use_anisotropy = FALSE) if necessary to aid convergence  |
| no vessel effects, catchability or density covariates in fit_model()  | may include covariates in cases where their incorporation has been previously demonstrated to improve model fit (e.g., a spatially varying response to cold-pool extent when generating abundance indices combining the EBS and NBS); in these cases, covariates will be centered and scaled prior to fitting  |
[^1]: Noting that for species with 100% encounters in any year we will use c(2,4) instead of c(2,1), or the equivalent setting for the lognormal <br/>
[^2]: Tweedie also involves additional changes to RhoConfig and FieldConfig to ensure that there is only a single linear predictor being estimated, as documented elsewhere

### **Note**: if a request to only include fish west of X degrees longitude, just note that there are two components to this request:  (1) excluding data from east of X degrees longitude and (2) specifying this as a "strata" boundary in the model settings in the VAST code

# R support
- Installation Guide for R 4.0.2, required to run MRAN [here](https://docs.google.com/document/d/1tjAjvVsYbRBYLWwVdQ-Bs7GYALiUn2xFcUgcP8mQCHw/edit?usp=sharing)
