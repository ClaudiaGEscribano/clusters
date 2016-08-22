## do the kmeans  with the centers of the hierarchical


library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(rgdal)
library(mapdata)
library(zoo)
library(ncdf4)

load("centroidshc.Rdata")
load("datafromPCA.Rdata")

############################################################################################
## 3. CENTERS MATRIX TO KMEANS
##########################################################################################

dataMatrix <- as.matrix(datahc)

centros <- lapply(seq(from=1, to=69),
                  FUN=function(i) do.call(rbind, centroids[[i]]))


kmeansexp <- lapply(centros,
                    FUN=function(i) kmeans(dataMatrix, centers=i, iter.max=3000))

#############################################################################################

save(kmeansexp, file='kmeansexpafterhc.Rdata')

