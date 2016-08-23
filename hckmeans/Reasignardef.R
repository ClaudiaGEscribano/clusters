## reasigno las celdas aisladas después de calcular los índices. de 16 a 22 clusters

library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(rgdal)
library(mapdata)
library(zoo)
library(ncdf4)

#############################################################################
## 1. LOAD THE DATA
###########################################################################

load("kmeansexpafterhc.Rdata")

s <- stack("SIS_complete_83_13.nc")

#############################################################################

P <- raster(s)
P <- setValues(P, kmeansexp[[21]]$cluster)

## Calculo los vecinos

ad <- adjacent(P, 1:ncell(P), directions=8, pairs=TRUE, id=TRUE, sorted=TRUE)

## cluster al que pertenecen los vecinos

vecinosValor <- lapply(seq(1:ncell(P)), FUN= function(x) P[ad[ad[,1] == x, 3]])

## Para comparar con el valor del cluster de la celda que estoy analizando utilizo una comparación lógica. Cuento cuantos de los vecinos cumplen la condición y en función de eso (if) reasigno. El valor al que reasigno viene dado por el valor más común que tengan los vecinos (utilizo la función table)


Reasignar <- lapply(seq(1:ncell(P)),
	FUN=function(i)  if (length(which(vecinosValor[[i]] != P[i])) >= 6) {
                 P[i] <- sort(vecinosValor[[i]], decreasing= TRUE)[1]
             } else { P[i] })


## Save the cluster partition with the hc+kmeans method

R <- raster(P)
R <- setValues(R, unlist(Reasignar))

writeRaster(R, filename='hckmeanspartition21.grd', overwrite=TRUE)

