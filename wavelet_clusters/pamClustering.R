## Aplicar algoritmo PAM a las características obtenidas con wavelet analysis

library(raster)
library(rasterVis)

## library of clustering

library(cluster)

## ¿aplico a una de las especies de aerosoles (or)?

pamx <- pam(or_features, k=20, diss=FALSE)

## kmeans sore las características

kmeansexp <- function(x, n, k){
    km <- lapply(seq(1:k),
                 FUN=function(i) kmeans(x, i, nstart=n, iter.max=1000))
    return(km)
}


## 2 especies juntas. algoritmo CLARA

features <- cbind(or_features, ss_features)

clarax <- clara(or_features, 20)

claraClusters <- clarax$clustering

