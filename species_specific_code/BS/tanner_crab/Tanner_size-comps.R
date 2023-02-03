

library(VAST)
library(googledrive)

root_dir = R'(C:\Users\James.Thorson\Desktop)'
#root_dir = R'(C:\Users\James.Thorson\Desktop\Work files\Collaborations\2023 -- Stockhausen tanner crab)'
data_dir = root_dir

dataset = c("Males")
# CPUE in numbers is in no/(sq. nm.)
# CPUE in weight is in mt/(sq. nm.)

Date = Sys.Date()
run_dir = paste0(root_dir, "/", Date, "_", dataset, "/" )
  dir.create(run_dir)

if( dataset=="Males" ){
  Data_drive <- drive_get("TannerCrab_CPUEbySizeBin_Males.csv", corpus = "allDrives" )
  drive_download( Data_drive$id, path=paste0(run_dir,'/TannerCrab_CPUEbySizeBin_Males.csv') )
  Data = read.csv( paste0(run_dir,'/TannerCrab_CPUEbySizeBin_Males.csv') )
}

#
category_names = sort(unique(Data$SIZE))

# Check balance by year-size ... looks good
table( Data[,'YEAR'],Data[,'SIZE'] )

# Check encounter probability by year-size
tapply( Data[,'numCPUE'], INDEX=list(Data[,'YEAR'],Data[,'SIZE']), FUN=function(x){sum(x>0)} )
tapply( Data[,'numCPUE'], INDEX=list(Data[,'SIZE']), FUN=function(x){sum(x>0)} )

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 50,
  Region = "eastern_bering_sea",
  purpose = "index2",
  ObsModel = c(2,4), # Deal with 0% encounters
  max_cells = 500 )

# Run model
fit = fit_model( settings = settings,
  Lat_i = Data[,'LATITUDE'],
  Lon_i = Data[,'LONGITUDE'],
  t_i = Data[,'YEAR'],
  b_i = as_units(Data[,'numCPUE'],"count"),
  a_i = set_units(as_units(rep(1,nrow(Data)),"nmile2"),"km2"),
  c_i = as.numeric(factor(Data[,'SIZE'],levels=category_names))-1,
  Npool = 2000,   # Will mirror hyperparameters for top-2 bins for Males
  newtonsteps = 0, # Speed up but looser convergence
  category_names = category_names,
  test_fit = FALSE, # Avoid checks for many betas
  working_dir = run_dir )

# Plot results
out = plot( fit,
  working_dir = run_dir )

# save proportions
Prop_ct = out$Proportions$Prop_ctl[,,1]
dimnames(Prop_ct) = list( fit$category_names, fit$year_labels )
  write.csv( Prop_ct, file=paste0(run_dir,"Prop_ct.csv") )
  
###############
# Unstratified area-swept
###############

Dens_ct = tapply( Data[,'numCPUE']/rep(1,nrow(Data)), 
                  INDEX = list( factor(Data[,'SIZE'],levels=category_names), Data[,'YEAR'] ),
                  FUN = mean )
Dens_ct = as_units(Dens_ct, "count / nmile2" )  
Dens_ct = set_units(Dens_ct, "count / km2" )

# 
Area = make_extrapolation_info( Region = "eastern_bering_sea" )
N_ct = set_units(Dens_ct * sum(Area$Area_km2), "count")

# 
dimnames(N_ct) = list( "Size"=category_names, "Year"=sort(unique(Data[,'YEAR'])) )
