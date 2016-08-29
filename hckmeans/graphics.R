library(raster)
library(rasterVis)
library(maptools)

## Read results
load('mascaraClusters20enTierra.Rdata')
## It's a list, convert it to a RasterStack
ksB <- stack(ksB)
## Collapse the RasterStack into a RasterLayer
r <- calc(ksB, sum, na.rm = TRUE)
r[r == 0] <- NA
## Define it as a categorical variable
r <- ratify(r)
rat <- levels(r)[[1]]
rat$cluster <- 1:20
levels(r) <- rat

## Download boundaries and explanatory variables from DIVA-GIS
old <- setwd(tempdir())
## Administrative boundaries
download.file('http://biogeo.ucdavis.edu/data/diva/adm/ESP_adm.zip',
              destfile = 'ESP_adm.zip')
unzip('ESP_adm.zip')
admin <- readShapeLines('ESP_adm0')

## Digital Elevation Model
download.file('http://biogeo.ucdavis.edu/data/diva/alt/ESP_alt.zip',
              destfile = 'ESP_alt.zip')
unzip('ESP_alt.zip')
dem <- raster('ESP_alt')
dem <- resample(dem, r)

## Land Cover
download.file('http://biogeo.ucdavis.edu/data/diva/cov/ESP_cov.zip',
              destfile = 'ESP_cov.zip')
unzip('ESP_cov.zip')
land <- raster('ESP1_cov')

## Use nearest neighbour because land cover is a categorical variable
land <- resample(land, r, method = 'ngb')
setwd(old)

## Display clusters (as a categorical variable)
levelplot(r) + layer(sp.lines(admin))
## Define an alternative colour palette
myPal <- brewer.pal(9, 'Set1')
myTheme <- rasterTheme(region = myPal)
levelplot(r, par.settings = myTheme) + layer(sp.lines(admin))

## Define a RasterStack with the explanatory variables 
s <- stack(r, dem, land) 
names(s) <- c('clusters', 'DEM', 'land')

## DEM and clusters
levelplot(dem)

xyplot(DEM ~ clusters, data = s)

bwplot(DEM ~ clusters, data = s)

histogram(~DEM|clusters, data = s)

## Land Cover and clusters
levelplot(land)
histogram(~land|clusters, data = s)

