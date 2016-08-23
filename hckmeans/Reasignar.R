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

load("kmeansexpafterhcmask.Rdata")

s <- stack("/home/claudia/Documentos/satelite/SIS_complete_83_13.nc")

ext <- as.vector(extent(s))
boundaries <- map('worldHires', region=c('Spain','Portugal', 'Andorra'),fill=TRUE, exact=TRUE, xlim=ext[1:2], ylim= ext[3:4], plot=FALSE)
boundaries$names
IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
boundaries_sp<- map2SpatialPolygons(boundaries, IDs=IDs, proj4string=CRS(projection(s)))

S <- subset(s, 1)
S <- mask(S, boundaries_sp)

#############################################################################


data <- as.data.frame(S)
namesNNA <- which(!is.na(data[,1]))
clusters_df<- as.data.frame(kmeansexp[[18]]$cluster)

data[namesNNA,] <- clusters_df[,1]

## data contiene ahora los clusters

P <- raster(S)
P <- setValues(P, data[,1])

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

writeRaster(R, filename='hckmeanspartition20mask.grd', overwrite=TRUE)

