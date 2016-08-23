## This script will do the Lmethod to find the optimun partition after the index and clustering

## tengo que correr este script para cada uno de los resultados que he sacado de índice

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

load('criterioDBhc.RData') 

###################################################################
## 2. Función ajuste
#################################################################

rmse <- list()

ajuste  <- function(x){
       for(j in 1:67){
        r1L <- lm(x[1:j+1]~c(1:j+1))
        r1R <- lm(x[j+2:70]~c(j+2:70))
	rmse[[j]] <- c(sqrt(sum((r1L$residuals)^2)), sqrt(sum((r1R$residuals)^2)))}
	return(rmse)
}

## I apply this 'ajuste' function to every element of 'Indice' (they are lists of kmeans results)

rmse_exp <- ajuste(as.vector(unlist(criterioDBhc)))

## Ponderate the rmse from lm results.
ajuste_ponderado  <- function(x){
	rmseT <- c()
	for(i in 1:67)
	rmseT[i] <-((i+1)-1)/((70)-1)*(x[[i]][1])+(70-(i+1))/(70-1)*(x[[i]][2])
	return(rmseT)
}

rmse_expP <- ajuste_ponderado(rmse_exp)

minimo <- which(rmse_expP == min(rmse_expP))

save(minimo, file='minimoDBhc.Rdata')
