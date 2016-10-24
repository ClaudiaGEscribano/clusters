library(raster)
library(zoo)

## datos de radiación diarios

SIS <- stack("/home/claudia/clusters/SIS_cmsaf30_complete_005.grd")

load("/home/claudia/clusters/mascaraClusters19enTierra_resolucion.Rdata")
## Convierto la lista de Raster a un RasterLayer
ksB$fun <- min
ksB <- do.call(mosaic, ksB)
## sequencia de dias

tt <- seq(as.Date("1983-01-01"), as.Date("2013-12-31"), 'day')
names(SIS) <- tt
## selecciono el mes que mas em interesa: Noviembre
month <- function(x) as.numeric(format(tt, "%m"))
idxNov <- which(month(tt) == 11)
ttNov <- tt[idxNov]
## tomo lo dias de noviembre de todo el stack de los 30 años
noviembreDiasSerie <- subset(SIS, idxNov)

## medias por clusters
meanCluster <- zonal(noviembreDiasSerie, ksB, 'mean')
## meanCluster es una matriz. La primera columna (zones) corresponde a los clusters.
## Genero una serie temporal, una fila por instante y una columna por zona
zNov <- zoo(t(meanCluster[, -1]), order = ttNov)
## Los nombres de esta serie los obtengo de la primera columna de la matriz meanCluster
names(zNov) <- paste0('C', meanCluster[, 1])

## Media por zonas
avgNov <- colMeans(zNov)
## Anomalía definida como la diferencia diaria en cada zona respecto de la media mensual
anomNov <- sweep(zNov, 2, avgNov)
## También en valores relativos a la media zonal
rAnomNov <- sweep(anomNov, 2, avgNov, '/')

saveRDS(anomNov, file = 'anomalias_noviembre.Rds')
saveRDS(rAnomNov, file = 'anomalias_rel_noviembre.Rds')

## Otro enfoque: descarto el año, y calculo la media de cada día del mes, y rehago cálculos. Por tanto, me fijo en el 1 de noviembre, el 2 de noviembre, etc.
day <- function(x) as.numeric(format(x, '%d'))
zDailyNov <- aggregate(zNov, by = day, FUN = 'mean')
avgDailyNov <- colMeans(zDailyNov)
anomDailyNov <- sweep(zDailyNov, 2, avgDailyNov)


## Compare anomaly with radiation
G <- stack(as.data.frame(zNov))
names(G) <- c('Rad', 'cluster')
A <- stack(as.data.frame(rAnomNov))
names(A) <- c('Anom', 'cluster')

M <- cbind(G, Anom = A[, "Anom"])

xyplot(Anom ~ Rad, groups = cluster, data = M)
