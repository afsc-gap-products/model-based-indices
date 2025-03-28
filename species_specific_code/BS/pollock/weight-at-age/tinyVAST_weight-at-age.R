
# https://github.com/jindivero/bs-pollock-weight/blob/main/scripts/2_model.R#L38-L71
# devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST)', force=TRUE )
# remotes::install_github( "vast-lib/tinyVAST@dev", force=TRUE )

library(tinyVAST)
library(fmesher)

root_dir = R'(C:\Users\James.Thorson\Desktop\Work files\Collaborations\2024 -- pollock weight-at-age)'

#Date = Sys.Date()
Date = "2024-10-16"
date_dir = file.path(root_dir, paste0(Date,"") )
  dir.create(date_dir)

#Data = readRDS( "data_combined.rds" )
Data = readRDS( file.path(root_dir, "Data_2024", "data_combined.rds") )

# Potentially subset to save time when exploring
#Data = subset( Data, age_bin==8 )
#Data = subset( Data, year>=2010 )

Data$data_type = ifelse( !is.na(Data$age_cpue_sum), "D", "W")
# table(Data$age_bin, Data$data_type)

# Add columns to define variables
Data$var = paste0( Data$data_type, Data$age_bin )
Data$Response = ifelse( Data$data_type=="D", Data$age_cpue_sum, Data$weight_combined )
Data$var_year = interaction( Data$var, Data$year, drop=TRUE )

# Define SEM
dsem = "
  D1 -> D1, 1, lagD
  D2 -> D2, 1, lagD
  D3 -> D3, 1, lagD
  D4 -> D4, 1, lagD
  D5 -> D5, 1, lagD
  D6 -> D6, 1, lagD
  D7 -> D7, 1, lagD
  D8 -> D8, 1, lagD
  D9 -> D9, 1, lagD
  D10 -> D10, 1, lagD
  D11 -> D11, 1, lagD
  D12 -> D12, 1, lagD
  D13 -> D13, 1, lagD
  D14 -> D14, 1, lagD
  D15 -> D15, 1, lagD

#  D1 -> D2, 1, cohortD
#  D2 -> D3, 1, cohortD
#  D3 -> D4, 1, cohortD
#  D4 -> D5, 1, cohortD
#  D5 -> D6, 1, cohortD
#  D6 -> D7, 1, cohortD
#  D7 -> D8, 1, cohortD
#  D8 -> D9, 1, cohortD
#  D9 -> D10, 1, cohortD
#  D10 -> D11, 1, cohortD
#  D11 -> D12, 1, cohortD
#  D12 -> D13, 1, cohortD
#  D13 -> D14, 1, cohortD
#  D14 -> D15, 1, cohortD

  W1 -> W1, 1, lagW
  W2 -> W2, 1, lagW
  W3 -> W3, 1, lagW
  W4 -> W4, 1, lagW
  W5 -> W5, 1, lagW
  W6 -> W6, 1, lagW
  W7 -> W7, 1, lagW
  W8 -> W8, 1, lagW
  W9 -> W9, 1, lagW
  W10 -> W10, 1, lagW
  W11 -> W11, 1, lagW
  W12 -> W12, 1, lagW
  W13 -> W13, 1, lagW
  W14 -> W14, 1, lagW
  W15 -> W15, 1, lagW

#  W1 -> W2, 1, cohortW
#  W2 -> W3, 1, cohortW
#  W3 -> W4, 1, cohortW
#  W4 -> W5, 1, cohortW
#  W5 -> W6, 1, cohortW
#  W6 -> W7, 1, cohortW
#  W7 -> W8, 1, cohortW
#  W8 -> W9, 1, cohortW
#  W9 -> W10, 1, cohortW
#  W10 -> W11, 1, cohortW
#  W11 -> W12, 1, cohortW
#  W12 -> W13, 1, cohortW
#  W13 -> W14, 1, cohortW
#  W14 -> W15, 1, cohortW
"

spatial_graph = fm_mesh_2d( loc = Data[,c('start_longitude','start_latitude')],
                            cutoff = 1.5 )
control = tinyVASTcontrol( getsd = FALSE,
                           nlminb_loops = 1,
                           newton_loops = 0,
                           profile = c("alpha_j"),
                           trace = 1 )
#family_link = cbind( rep(c(1,2),each=length(unique(Data$age_bin))), 1 )
#  rownames(family_link) = paste0( rep(c("D","W"),each=length(unique(Data$age_bin))),
#                                  rep(sort(unique(Data$age_bin)),2) )
family = NULL
for( i in seq_along(unique(Data$age_bin)) ){
  family[[i]] = tweedie()
}
for( i in (length(unique(Data$age_bin))+seq_along(unique(Data$age_bin))) ){
  family[[i]] = lognormal()
}
names(family) = paste0( rep(c("D","W"),each=length(unique(Data$age_bin))),
                                  rep(sort(unique(Data$age_bin)),2) )

if( !("myfit.RDS" %in% list.files(date_dir)) ){
  myfit = tinyVAST(
    data = Data,
    #formula = Response ~ 0 + factor(var) + s(year, by=factor(var), bs="ts", k=9),   # bs="ts" is more robust
    formula = Response ~ 0 + var_year,
    dsem = dsem,
    family = family,
    space_columns = c("start_longitude", "start_latitude"),
    variable_column = "var", 
    variables = names(family),
    time_column = "year", 
    distribution_column = "var",
    spatial_graph = spatial_graph,
    control = control
  )
  
  saveRDS( myfit, 
        file = file.path(date_dir,"myfit.RDS") )
  capture.output( myfit, 
                  file = file.path(date_dir,"myfit.txt") )
}else{
  myfit = readRDS( file = file.path(date_dir,"myfit.RDS") )
  myfit = reload_model(myfit)
}

# DHARMa residuals
# simulate new data conditional on fixed and random effects
y_ir = replicate( n = 100,
           expr = myfit$obj$simulate()$y_i )

#
res = DHARMa::createDHARMa( simulatedResponse = y_ir,
                            observedResponse = Data$Response,
                            fittedPredictedResponse = fitted(myfit) )
png( file.path(date_dir,"residuals=data_type.png"), width=8, height=4, res=200, units="in")
  par( mar=c(3,3,2,2), mgp=c(2,0.5,0), tck=-0.02 )
  plot( res, form = Data$data_type )
dev.off()
png( file.path(date_dir,"residuals=data_type+year.png"), width=8, height=4, res=200, units="in")
  par( mar=c(3,3,2,2), mgp=c(2,0.5,0), tck=-0.02 )
  plot( res, form = interaction(Data$data_type,Data$year) )
dev.off()

################
# Calculate indices
###############

# Get polygon
library(sf)
ebs = st_read( R'(C:\Users\James.Thorson\Desktop\Git\FishStatsUtils\inst\region_shapefiles\EBSshelf)' )
nbs = st_read( R'(C:\Users\James.Thorson\Desktop\Git\FishStatsUtils\inst\region_shapefiles\NBS)' )
bs = st_geometry(st_union( ebs, nbs ))
bs = st_transform( bs, crs=st_crs("+proj=longlat +datum=WGS84") )

# Make extrapolation grid
grid = st_make_grid( bs, n=c(10,10) )
grid = st_intersection( grid, bs )
grid = st_make_valid( grid )
loc_gz = st_coordinates(st_centroid( grid ))
colnames(loc_gz) = c("start_longitude", "start_latitude")

# Get area for extrapolation grid
library(units)
areas = set_units(st_area(grid), "hectares") #  / 100^2 # Hectares

# Get abundance
W_jz = expand.grid( Age=sort(unique(Data$age_bin)), year=myfit$internal$times )
W_jz = cbind( W_jz, "Weight"=NA, "SE"=NA )
for( j in seq_len(nrow(W_jz)) ){
  gc()
  if( is.na(W_jz[j,'Weight']) ){
    message( "Integrating ", W_jz[j,'year'], " ", W_jz[j,'Age'], ": ", Sys.time() )
    # Manual epsilon
    newdata = rbind(
      data.frame( loc_gz, year=W_jz[j,'year'], var=paste0("D",W_jz[j,'Age']), var_year=paste0("D",W_jz[j,'Age'],".",W_jz[j,'year']) ),
      data.frame( loc_gz, year=W_jz[j,'year'], var=paste0("W",W_jz[j,'Age']), var_year=paste0("W",W_jz[j,'Age'],".",W_jz[j,'year']) )
    )
    # predict( myfit, newdata = newdata )
    if( paste0("W",W_jz[j,'Age'],".",W_jz[j,'year']) %in% Data$var_year ){
      area = c(as.numeric(areas), 0*as.numeric(areas))
      type = rep(c(0,3), each=length(areas))
      weighting_index = c( rep(0,length(areas)), seq_along(areas)-1 )
      #source( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\R\internal.R)' )
      index1 = integrate_output( myfit,
                      area = area,
                      newdata = newdata,
                      type = type,
                      weighting_index = weighting_index,
                      apply.epsilon = TRUE,
                      bias.correct = FALSE,
                      intern = TRUE )
      W_jz[j,'Weight'] = ifelse( is.na(index1[3]), index1[1], index1[3] ) # Depends on bias-correction settings above
    }
  }
}
W_ct = array( W_jz$Weight, dim=c("Age"=length(unique(Data$age_bin)),"year"=length(myfit$internal$times)),
              dimnames=list(sort(unique(Data$age_bin)),myfit$internal$times) )
#image( z=t(W_ct), x=myfit$internal$times, y=seq_along(myfit$internal$variables) )
write.csv( W_ct, file=file.path(date_dir,"W_ct.csv") )
write.csv( W_jz, file=file.path(date_dir,"W_jz.csv") )

#
Wmean_ct = tapply( Data$weight_combined, INDEX=list("Age"=Data$age_bin,"Year"=Data$year), FUN=mean, na.rm=TRUE )
Wmean_jz = data.frame( expand.grid(dimnames(Wmean_ct)), "Weight"=as.vector(Wmean_ct) )
write.csv( Wmean_jz, file="Wmean_jz.csv" )

# Plot densities
library(viridisLite)
source( R'(C:\Users\James.Thorson\Desktop\Git\Spatio-temporal-models-for-ecologists\Shared_functions\add_legend.R)' )
for(data_type in c("D","W") ){
  var_set = paste0(data_type,1:15)
  year_set = sort(unique(Data$year))
  #
  X_jz =  expand.grid( "cell"=1:nrow(loc_gz), var=var_set, year=myfit$internal$times )
  X_jz$var_year = paste0( X_jz$var, ".", X_jz$year )
  X_jz = subset( X_jz, var_year %in% unique(subset(Data,data_type==data_type)$var_year) )
  X_jz = cbind( X_jz, loc_gz[X_jz$cell,] )
  X_jz$density = predict( myfit, newdata = X_jz, remove_origdata = FALSE )
  gc()
  # Plot
  mfrow = ceiling(sqrt(length(year_set)))
    mfrow = c( mfrow, ceiling(length(year_set)/mfrow) )
  for( var in var_set ){
    png( file=file.path(date_dir,paste0(var,".png")), width=mfrow[2]*3, height=mfrow[1]*3, res=200, units="in" )
      par( mfrow=mfrow, mar=c(0,0,2,0) )
      #d_vec = X_jz[which(X_jz$var==var),'density']
      #dims = c( nrow(loc_gz), length(year_set) )
      #if(length(d_vec) != prod(dims)) stop("Check dimensions")
      #d_gt = array( d_vec, dim=dims )
      tmp_jz = X_jz[which(X_jz$var==var),]
      d_gt = tapply( tmp_jz$density, INDEX=list(tmp_jz$cell,factor(tmp_jz$year,year_set)), FUN=mean )
      d_gt = ifelse( d_gt < max(d_gt,na.rm=TRUE)/1000, NA, d_gt )
      breaks = seq( min(log(d_gt),na.rm=TRUE), max(log(d_gt),na.rm=TRUE), length=11 )
      for( tI in seq_along(year_set) ){
        plotgrid = st_sf(grid, log(d_gt[,tI]), crs=st_crs(grid) )
        plot(plotgrid, max.plot=20, border=NA, key.pos=NULL, reset=FALSE,
             pal=viridis, main=year_set[tI], breaks=breaks )
        plot( bs, add=TRUE )
      }
      add_legend( round(range(breaks),2), legend_y=c(0.6,1), legend_x=c(1,1.05), col=viridis(10) )
    dev.off()
  }
}
