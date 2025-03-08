<img src="https://github.com/user-attachments/assets/e0fbc2b1-aa23-42fe-8a71-2b0862f9ede2" width="200" height="200" />

# model-based-indices
Code used for model-based indices of abundance and distribution from bottom trawl survey data, produced by the AFSC GAP ModSquad in the Gulf of Alaska and Bering Sea. Model-based indices of abundance and distribution are produced using spatiotemporal delta-GLMMs implemented in [sdmTMB](https://pbs-assess.github.io/sdmTMB/). Model-based age composition estimates follow a similar approach, but employ a multivariate extension and are implemented in [tinyVAST](https://vast-lib.github.io/tinyVAST/).

## Authors
Primary contact: [@Lewis-Barnett-NOAA](https://github.com/Lewis-Barnett-NOAA)

# Mission
Provide model-based estimates of abundance and their standard error from GAP survey data for stock assessment authors in conjunction with traditional design-based estimators.
- [TOR or Terms Of Reference](https://drive.google.com/file/d/1aVv-lbh-l_Q_zCnjbyr_38jH5EuVos1K/view?usp=drive_link)
- [Annual product request form](https://docs.google.com/spreadsheets/d/18gr3owj5iAq1iCDX4wpQPUC9ldLz-YTsCBIfnkHqibo/edit?usp=sharing)
- [ESP submission tool](https://apex.psmfc.org/akfin/f?p=140:LOGIN_DESKTOP:4779711459935:::::)

# GAP Work Plan 
- Species requests & alternate model-based runs/settings in from leads from MESA, SSMA, ESP, ESR by 15 January
- Run "hindcasts" with previous years of survey data to finalize models prior to inclusion of new data (February 15 - April, code frozen by May)
- Merge new frozen code (index & data retrieval code) in git repo
- Upload all output files for hindcast and data used to produce the indices to the appropriate region hindcast folder in the corresponding year of the ModSquad [google drive folder](https://drive.google.com/drive/folders/1yxn02yF0V1PNVw0_HpqeSK_XxgOy_LAT)
- Notify the assessment lead & ModSquad lead (Lewis) that you’ve completed hindcasts (March - April)
- Run hindcast model structure again after survey with updated and validated data (August - Sept, completed by 30 September) and upload results and data used to produce the indices to the appropriate region production folder within the corresponding year folder [on google drive](https://drive.google.com/drive/folders/1yxn02yF0V1PNVw0_HpqeSK_XxgOy_LAT)
- Notify the assessment lead & ModSquad lead (Lewis) that you’ve completed the production estimates

# General Production Timeline Estimate
**Feb 15 - April is the time window for any model exploration/iteration, as requested**
  - December - ModSquad leadership meets/communicates to prioritize species + research plans
  - January 15 - Deadline for SAFE authors to submit product requests
  - late January - ModSquad planning meeting with program leads as needed
  - February 15 - SAFE author's revised requests given ModSquad capacity and discussion from REFM meeting
  - March 1 - TOR memo finalized and posted on the web
  - April 1 - Crab hindcasts completed
  - May 1 - Groundfish hindcasts completed (stop all research, decide on which to recommend or not recommend for use in base stock assessment model)
  - August 25 - Model-based estimates for EBS crabs completion target date
  - September 25 - Model-based estimates for NBS crabs target date
  - September 30 - Model-based estimates for groundfishes completion target date, pending bottom trawl data QAQC timeline
  - October 15 - Final deadline for completing groundfish model-based estimates 

# Resources 
- [ModSquad folder](https://drive.google.com/drive/folders/1yxn02yF0V1PNVw0_HpqeSK_XxgOy_LAT)
- [ModSquad training materials](https://drive.google.com/drive/folders/1TZRrEwka7OEICC6D81LiuGtYJMSxdaFe?usp=sharing)
- [Common troubleshooting guidance](https://docs.google.com/document/d/1j3Li2aacvy7d4FJxLlGDctHlDJWQgZzI5HzyMujwf8Y/edit?usp=sharing)
- Also see vignettes and workshop materials for sdmTMB and tinyVAST

# Minimum standard of versions of software and key dependencies (i.e., using these versions or later)
- **2025**: 
  - R v4.0.2: sdmTMB v0.6.0, tinyVAST v0.7.0-alpha, TMB v1.7.22, Matrix v1.4-0, DHARMa v0.4.5; or,
  - MRAN v4.0.2: sdmTMB v0.5.0, tinyVAST v0.1.0-alpha, TMB v1.7.16, Matrix v1.2-18, DHARMa v0.3.2
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

## NOAA README

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.

## NOAA License

Software code created by U.S. Government employees is not subject to
copyright in the United States (17 U.S.C. §105). The United
States/Department of Commerce reserve all rights to seek and obtain
copyright protection in countries other than the United States for
Software authored in its entirety by the Department of Commerce. To this
end, the Department of Commerce hereby grants to Recipient a
royalty-free, nonexclusive license to use, copy, and create derivative
works of the Software outside of the United States.

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" alt="NOAA Fisheries" height="75"/>

[U.S. Department of Commerce](https://www.commerce.gov/) \| [National
Oceanographic and Atmospheric Administration](https://www.noaa.gov) \|
[NOAA Fisheries](https://www.fisheries.noaa.gov/)
