## pruebas con las climatologias de aerosoles para hacer cluster con las wavelets. Pruebo con los datos de la malla LCC y sólo para representar los proyectos.

library(raster)
library(rasterVis)
library(waveslim)
library(zoo)

###############################################
## Sin proyectar a lat lon el nc primero.

wave2cluster <- function(x) {
    data <- as.data.frame(x)

    dataZoo <- zoo(t(data))

    waveDecomposition <- apply(dataZoo, 2, FUN=function(x) modwt(x, n.levels=6))
    waveVariance <- lapply(waveDecomposition, FUN=function(x) wave.variance(x))

    features <- lapply(1:ncol(dataZoo), FUN=function(x)
                                 waveVariance[[x]]$wavevar[-7]/var(dataZoo[,x]))
    return(features)
}

################################################
## Para descomponer las series temporales con wavelet es necesario centrarlas en el 0, por lo que calcularemos las anomalías del periodo y haremos la descomposición sobre eso

wave2cluster <- function(x) {
    data <- as.data.frame(x)
    dataZoo <- zoo(t(data))

    ## series de anomalías
    M <- colMeans(dataZoo)
    #dataZooM <- apply(dataZoo, 2, FUN=function(x) x- M[[x]])
    a <-    lapply(seq(1:ncol(dataZoo)), FUN= function(x) dataZoo[,x]-M[[x]])
    al <- do.call(cbind, a)
    ## Aplico la descomposición wavelet

    waveDecomposition <- apply(al, 2, FUN=function(x) modwt(x, n.levels=6)) ## el numero de niveles depende de la longitud
    waveVariance <- lapply(waveDecomposition, FUN=function(x) wave.variance(x))
 
    features <- lapply(1:ncol(dataZoo), FUN=function(x)
                                 waveVariance[[x]]$wavevar[-7]/var(al[,x]))
    return(features)
}

#################################################
## extraigo cada uno de las especies de aerosoles por separado

## bc

## donde están los datos
## "/aerosoles_DATA/climatologia_AOD/AOD/bc"

listFich <- dir(path="/data/bc", pattern='\\.nc')

bc <- stack(listFich, varname='aero')

listafeatures_bc<- wave2cluster(bc)

bc_features <- do.call(rbind, listafeatures_bc)
colnames(bc_features) <- c("f1","f2","f3", "f4","f5", "f6")
save(bc_features, file='bc_features.Rdata')

## sd

## donde están los datos
## "/aerosoles_DATA/climatologia_AOD/AOD/sd"

listFich <- dir(path="/data/sd", pattern='\\.nc')

sd <- stack(listFich, varname='aero')

listafeatures_sd<- wave2cluster(sd)

sd_features <- do.call(rbind, listafeatures_sd)
colnames(sd_features) <- c("f1","f2","f3", "f4","f5", "f6")
save(sd_features, file='sd_features.Rdata')

## or

## done están los datos
## "/aerosoles_DATA/climatologia_AOD/AOD/or"

listFich <- dir(path="/data/or", pattern='\\.nc')

or <- stack(listFich, varname='aero')

listafeatures_or<- wave2cluster(or)

or_features <- do.call(rbind, listafeatures_or)
colnames(or_features) <- c("f1","f2","f3", "f4","f5", "f6")
save(or_features, file='or_features.Rdata')

## ss

## donde están los datos
## "/aerosoles_DATA/climatologia_AOD/AOD/ss"

listFich <- dir(path="/data/ss", pattern="^macc_regcm_2.*\\.nc$")

ss <- stack(listFich, varname='aero')

listafeatures_ss<- wave2cluster(ss)

ss_features <- do.call(rbind, listafeatures_ss)
colnames(ss_features) <- c("f1","f2","f3", "f4","f5", "f6")
save(ss_features, file='ss_features.Rdata')

## su

## donde estan los datos:
## "aerosoles_DATA/climatologia_AOD/AOD/su"

listFich <- dir(path="/data/su", pattern="\\.nc")

su <- stack(listFich, varname='aero')

listafeatures_su<- wave2cluster(su)

su_features <- do.call(rbind, listafeatures_su)
colnames(su_features) <- c("f1","f2","f3", "f4","f5", "f6")
save(su_features, file='su_features.Rdata')
############################################################################

## pruebo a aplicar wave2features a AOD.

aod <- stack("data/AOD_total_monthly20032009.grd")

listafeaturesaod <- wave2cluster(aod)
aod_features <- do.call(rbind, listafeaturesaod)
colnames(aod_features) <- c("f1","f2","f3", "f4","f5", "f6")
save(aod_features, file='aod_features.Rdata')

## aod_features es el objeto que utilizo para hacer los clusters.
