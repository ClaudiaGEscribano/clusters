## pruebas con las climatologias de aerosoles para hacer cluster con las wavelets. Pruebo con los datos de la malla LCC y sólo para representar los proyectos.

library(raster)
library(rasterVis)
library(waveslim)
library(zoo)

###############################################
MODWT <- function(x, filter = 'haar', n.levels = 6){
    y <- modwt(x, wf = filter, n.levels = n.levels)
    ## Circular shift: needed to align them with the time series
    y <- phase.shift(y, filter)
    ## modwt provides a list: convert it into a data.frame
    w <- do.call(cbind, y)
    w
}

## Standard deviation with sign (average sign in the interval)
sdSign <- function(x, na.rm = TRUE){
    s <- sign(mean(x, na.rm = na.rm))
    s * sd(x, na.rm = na.rm)
}

wavVar <- function(x, ...)
{
    wx <- MODWT(x, ...)
    apply(wx[, -ncol(wx)], 2, sdSign)
}


#################################################
## extraigo cada uno de las especies de aerosoles por separado

## bc

## donde están los datos
## "/aerosoles_DATA/climatologia_AOD/AOD/bc"
old <- setwd("/home/claudia/aerosoles_DATA/climatologia_AOD/AOD/bc")

listFich <- dir(pattern='\\.nc')

bc <- stack(listFich, varname='aero')

bc_features<- calc(bc, wavVar)

save(bc_features, file='bc_features.Rdata')

setwd(old)
