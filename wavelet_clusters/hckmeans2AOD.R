## Este script hace el algoritmo de clustering pca + hierarchical + kmeans sobre la climatología de aerosoles

library(raster)
library(rasterVis)

## trying the kmeans algorithm first:
 
set.seed(10)

mync <- stack("data/AOD_total_monthly20032009.grd")

## pca y hc clustering para inicializar el kmeans

################################################################
## 2. MÁSCARA DE TIERRA
################################################################

mascara <- raster("data/masque_terre_mer.nc", varname='zon_new')

## mask the .nc aod with the mask of the land

mycrs <- CRS("+proj=lcc +lat_1=43 +lat_0=43 +lon_0=15 +k=0.684241 +units=m +datum=WGS84 +no_defs")

maslat <- raster("data/masque_terre_mer.nc", varname='lat')
maslon <- raster("data/masque_terre_mer.nc", varname='lon')

pmaslat <- rasterToPoints(maslat)
pmaslon <- rasterToPoints(maslon)
maslonlat <- cbind(pmaslon[,3], pmaslat[,3])

## Specify the lonlat as spatial points with projection as long/lat
maslonlat <- SpatialPoints(maslonlat, proj4string = CRS("+proj=longlat +datum=WGS84"))

maslonlat
extent(maslonlat)

pmaslonlat <- spTransform(maslonlat, CRSobj = mycrs)
## Take a look                                                                                                                                                       
pmaslonlat
extent(pmaslonlat)

projection(mascara) <- mycrs
extent(mascara) <- extent(pmaslonlat)

######################################################################
## 3. Para pasar la mascara a los datos AOD
#######################################################################

newproj <- projectExtent(mascara, mycrs)
aodproj <- projectRaster(mync, newproj)

data <- mask(aodproj, mask=mascara, maskvalue=0)

############################################################
## 4. MAPA CON LOS BORDES
##########################################################

library(maps)
library(maptools)
library(mapdata)

data(worldMapEnv)

crslonlat <- CRS("+proj=longlat +datum=WGS84")
crslcc <- CRS("+proj=lcc +lat_1=43 +lat_0=43 +lon_0=15 +k=0.684241 +units=m +datum=WGS84 +no_defs")

ext <- as.vector(extent(projectExtent(aodproj, crs=crslonlat)))
boundaries <- map('world', fill=TRUE, exact=FALSE, xlim=ext[1:2], ylim= ext[3:4], plot=FALSE)
                                        #boundaries$names
IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
boundaries_sp<- map2SpatialPolygons(boundaries, IDs=IDs, proj4string=crslonlat)

border <- as( boundaries_sp, 'SpatialLines')

boundaries_lcc <- spTransform(boundaries_sp, crslcc)
border_lcc <- spTransform(border, crslcc)

## representación del raster AOD con la mascara tierra-mar y los bordes de continente:

levelplot(data, layers=1)+layer(sp.lines(border_lcc))

#################################################################
## 5. pca
##################################################################

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

datahcmask <- pca2cluster(data)
 
#save(datahcmask, file='datahcmaskPCA.Rdata')

###############################################################
## 6. hclust
#############################################################

dataMatrix <- as.matrix(datahcmask)
dataDist <- dist(dataMatrix) ## computes the distance matrix in order to use it in the hc algorithm

datahclus <- hclust(dataDist, method="ward.D2")

## I get the cut the clustering tree for k=2:70 clusters. I will decide the optimun number in next steps.

kut <- cutree(datahclus, k=2:70)

#################################################################
## 7. CENTROIDS
###############################################################

## I have to select the centroids of the clusters in order to use them in the kmeans clustering.
 
cluster.centroid <- function(datahc, data, n){ ## n es el numero de clusters max. menos uno.
     lapply(seq(from=1, to=n),
            FUN=function(x) lapply(seq(from=1, to=x+1),
                FUN=function(i){
                    ind  <- which(datahc[,x]==i) ## cambio aqui  i por x
                    r <- data[ind,]
                    c <- colMeans(r)
                    return(c)
                })
            )
 }


centroids <- cluster.centroid(kut, dataMatrix, 69) ## contiene una lista de 70 elementos. Cada uno de los elementos es otra lista con los valores de de los centroides, que son vectores de dimensión d igual a las columnas de data.
## saving centroids, just if you need to check

save(centroids, file='centroidshcmask.Rdata')

############################################################################################
## 8. CENTERS MATRIX TO KMEANS
##########################################################################################

centros <- lapply(seq(from=1, to=69),
                  FUN=function(i) do.call(rbind, centroids[[i]]))


kmeansexp <- lapply(centros,
                    FUN=function(i) kmeans(dataMatrix, centers=i, iter.max=3000))

#############################################################################################

## 8.b KMEANS WITHOUT CENTERS:

kmeansexp <- function(x, n, k){
    km <- lapply(seq(1:k),
                 FUN=function(i) kmeans(x, i, nstart=n, iter.max=1000))
    return(km)
}

a <- kmeansexp(datahcmask, 100, 50)
a <- kmeansexp(data2, 100, 50)

#################################################
## SELECT SOME PARTITIONS TO PUT THEM ON A RASTER

## Supongo la particion optima la de 20 rasters
 
##cluster20 <- raster(data)
##cluster20 <- setValues(cluster20, kmeansexp[[17]]$cluster)

##pdf("clusters17_2.pdf")
##levelplot(cluster20)
##dev.off()

##################################################
## 9. CRITERIO PARA NUMERO DE CLUSTERS
##################################################
 
library(clusterCrit)

datakmmatrix <- as.matrix(datahcmask)
  
criterioCHhc <- lapply(seq(from=1, to=69),
                       FUN= function(x) intCriteria(datakmmatrix, kmeansexp[[x]]$cluster, 'Calinski_Harabasz')
                       )

################################################
## 10. LMETHOD
################################################

rmse <- list()

ajuste  <- function(x){
    for(j in 1:67){
        r1L <- lm(x[1:j+1]~c(1:j+1))
        r1R <- lm(x[j+2:70]~c(j+2:70))
        rmse[[j]] <- c(sqrt(sum((r1L$residuals)^2)), sqrt(sum((r1R$residuals)^2)))}
    return(rmse)
}

## I apply this 'ajuste' function to every element of 'Indice' (they are lists of kmeans results)

rmse_exp <- ajuste(as.vector(unlist(criterioCHhc)))

## Ponderate the rmse from lm results.
ajuste_ponderado  <- function(x){
    rmseT <- c()
    for(i in 1:67)
        rmseT[i] <-((i+1)-1)/((70)-1)*(x[[i]][1])+(70-(i+1))/(70-1)*(x[[i]][2])
    return(rmseT)
}

rmse_expP <- ajuste_ponderado(rmse_exp)

minimo <- which(rmse_expP == min(rmse_expP))

###############################################################################
## 11. REASIGNAR
###########################################################################

## 11.a necesito que aparezcan todas las celdas, incluidas las que tenían NA.

dataN <- as.data.frame(data)
namesNNA <- which(!is.na(dataN[,1]))
clusters_df<- as.data.frame(kmeansexp[[15]]$cluster) ## tomo uno menos del óptimo porque kmeansexp[[1]] coniene 2 clusters

dataN[namesNNA,] <- clusters_df[,1]

cluster16 <- raster(data)
cluster16 <- setValues(cluster16, dataN[,1])

## con kmeans:

## 25 DB
## 17 CH

## minimo con mascara
## CH 16

## 11.b reasigno

## Calculo los vecinos

ad <- adjacent(cluster16, 1:ncell(cluster16), directions=8, pairs=TRUE, id=TRUE, sorted=TRUE)

## cluster al que pertenecen los vecinos

vecinosValor <- lapply(seq(1:ncell(cluster16)), FUN= function(x) cluster16[ad[ad[,1] == x, 3]])

## Para comparar con el valor del cluster de la celda que estoy analizando utilizo una comparación lógica. Cuento cuantos de los vecinos cumplen la condición y en función de eso (if) reasigno. El valor al que reasigno viene dado por el valor más común que tengan los vecinos (utilizo la función table)


Reasignar <- lapply(seq(1:ncell(cluster16)),
	FUN=function(i)  if (length(which(vecinosValor[[i]] != cluster16[i])) >= 5) {
                 cluster16[i] <- sort(vecinosValor[[i]], decreasing= TRUE)[1]
             } else { cluster16[i] })


## Save the cluster partition with the hc+kmeans method

R <- raster(cluster16)
R <- setValues(R, unlist(Reasignar))

## PRIMERA APROX: kmeans, ch=16 + reasignar

pdf("aod_clusters_16_ch.pdf")
levelplot(R)+layer(sp.lines(border_lcc, lwd=0.2))
dev.off()

####################
## 12 VISUAL
####################


Polygons <- list()
for (i in 1:16) Polygons[[i]] <- rasterToPolygons(R, fun=function(x) {x==i})
 
levelplot(mask(R, Polygons[[15]]))+layer(sp.lines(border_lcc, lwd=0.3)) ## visualiza los clusters uno a uno.

ksB <- lapply(seq(1:20), FUN=function(x) mask(R, Polygons[[x]]))
            
save(ksB, file='mascaraClusters20enTierra.Rdata')


