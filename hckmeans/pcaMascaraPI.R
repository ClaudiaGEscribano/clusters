## Este script hace un pca para los datos de radiaición de la península haciendo una máscara con los bordes

library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(rgdal)
library(mapdata)
library(zoo)
library(ncdf4)

#################################################################
## 1. LOAD THE DATA
###############################################################

SISd <- stack("/home/claudia/Documentos/satelite/SIS_complete_83_13.nc")

################################################################
## 2. MÁSCARA DE LA PI
################################################################

ext <- as.vector(extent(SISd))
boundaries <- map('worldHires', region=c('Spain','Portugal', 'Andorra'),fill=TRUE, exact=TRUE, xlim=ext[1:2], ylim= ext[3:4], plot=FALSE)
boundaries$names
IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
boundaries_sp<- map2SpatialPolygons(boundaries, IDs=IDs, proj4string=CRS(projection(SISd)))

SISd <- mask(SISd, boundaries_sp)

################################################################
## 2. PCA
###############################################################

pca2cluster <- function(x) {
    data <- as.data.frame(x)
    
    ## elimino los NA del dataframe

	namesNNA <- which(!is.na(data[,1]))
	dataNNA <- data[namesNNA,]
	##dataNNA <- as.data.frame(dataNNA)
	##rownames(dataNNA) <- namesNNA
    
   ## Hago la descomposición sólo para los datos con valores disintos a na

    datapca <- prcomp(dataNNA, center=TRUE, scale=TRUE) ## PCA

    ## I select the number of PC that I consider

    cumvar <- cumsum(datapca$sdev^2)/sum(datapca$sdev^2)
    b <- which(cumvar < 0.95)
    c <- length(b)
    datapca2 <- data.frame(datapca$x[,1:c])
    return(datapca2)
}

datahcmask <- pca2cluster(SISd)

save(datahcmask, file='datahcmaskPCA.Rdata')
