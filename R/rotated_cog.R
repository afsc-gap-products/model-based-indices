#' Calculated COG on a rotated axis for GOA ESP requests
#' 
#' @param prediction object from predict.sdmTMB()
#' @param rotated_radians numeric, rotation in radians, default at -pi/8
#' 

rotated_cog <- function(prediction, rotation_radians = - pi / 8){
  do.call(
    what = rbind, 
    args = lapply(
      X = split(x = prediction$data, 
                f = prediction$data$year_f),
      FUN = function(df) {
        
        rotated_grid <- data.frame(
          X = scale(x = cos(x = rotation_radians) * df$X -
                      sin(x = rotation_radians) * df$Y,
                    scale = FALSE),
          Y = scale(x = sin(x = rotation_radians) * df$X +
                      cos(x = rotation_radians) * df$Y,
                    scale = FALSE),
          X = df$X,
          Y = df$Y,
          area_km2 = df$area_km2,
          est1 = df$est1,
          est2 = df$est2
        )
        
        ## Average eastings and northings weighted by predicted biomass
        apply(X = rotated_grid[, c("X", "Y")],
              MARGIN = 2,
              FUN = weighted.mean,
              w = with(rotated_grid, 
                       exp(est1 + est2) * area_km2 ))
      }
    )
  )
}
