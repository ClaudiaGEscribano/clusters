## sCRIOT PARA DETERMINAR LA PARTICIÓN ÓPTIMA CON HC+KMEANS

library(raster)
library(rasterVis)
library(ncdf4)
library(maps)
library(maptools)
library(rgdal)
library(ncdf4)

#######################################################################
## 1. LOAD THE DATA
#######################################################################

load("kmeansexpafterhcmask.Rdata")

#######################################################################
## function to analyse the correct number of cluters.
######################################################################

library(clusterCrit)


load('datahcmaskPCA.Rdata')
datakm_matrix <- as.matrix(datahcmask)


criterioDBhc <- lapply(seq(from=1, to=69), 
	FUN= function(x) intCriteria(datakm_matrix, kmeansexp[[x]]$cluster, 'Davies_Bouldin')
                     )

save(criterioDBhc, file='criterioDBhcmask.RData')
