## This file is trying to read the .nc files with the LCC projection not defined. 
library(raster)
library(rasterVis)
library(ncdf4)


###########################################################################
## 1. project to lat lon grid all the files
##########################################################################

## First, estimate the extent.
## We start with the lat and long, then project it to the Lambert Conformal projection

inputfile <- "../macc_regcm_2003_01_bc.nc"

## Grab the lat and lon from the data
lat <- raster(inputfile, varname="lat")
lon <- raster(inputfile, varname="lon")

## Convert to points and match the lat and lons
plat <- rasterToPoints(lat)
plon <- rasterToPoints(lon)
lonlat <- cbind(plon[,3], plat[,3])

## Specify the lonlat as spatial points with projection as long/lat
lonlat <- SpatialPoints(lonlat, proj4string = CRS("+proj=longlat +datum=WGS84"))

## Need the rgdal package to project it to the original coordinate system
library("rgdal")

## My best guess at the proj4 string from the information given
mycrs <- CRS("+proj=lcc +lat_1=43.f +lat_0=43.f +lon_0=15.f +k=0.684241 +units=m +datum=WGS84 +no_defs")
plonlat <- spTransform(lonlat, CRSobj = mycrs)
## Take a look
plonlat
extent(plonlat)

## Now we can properly set the coordinate information for the rasterStack with all the time-steps.

##########################################################################
## 2. different aerosols
##########################################################################

## 2.1 BC ##

updir <- "/disco2/aerosoles_DATA/climatologia_AOD/AOD/bc"
direc <- "/disco2/aerosoles_DATA/climatologia_AOD/AOD/bc/calc"

## List with the files I want to merge in a stack, not neccessary previous work with cdo because I have a list of files one per time step.

listFich <- dir(path="/disco2/aerosoles_DATA/climatologia_AOD/AOD/bc", pattern='\\.nc')
setwd(updir) 

## pr has the rasterStack with the monthly data of bc aerosols.
pr <- stack(listFich, varname="aero")
setwd(direc)

## Fix the projection and extent
projection(pr) <- mycrs
extent(pr) <- extent(plonlat)
## Take a look
pr
plot(pr)

## Project to long lat grid
bc <- projectRaster(pr, crs=CRS("+proj=longlat +datum=WGS84"))
writeRaster(bc, filename="bc_monthlyclim.grd")
## Take a look
r
plot(r)
## Add contours
contour(r, add=TRUE)

## Add country lines
library("maps")
map(add=TRUE, col="blue")


## 2.2 SC ## 

updir <- "/disco2/aerosoles_DATA/climatologia_AOD/AOD/or"
direc <- "/disco2/aerosoles_DATA/climatologia_AOD/AOD/bc/calc"

## List with the files I want to merge in a stack, not neccessary previous work with cdo because I have a list of files one per time step.

listFich <- dir(path="/disco2/aerosoles_DATA/climatologia_AOD/AOD/or", pattern='\\.nc')
setwd(updir) 

## pr has the rasterStack with the monthly data of bc aerosols.
pr <- stack(listFich, varname="aero")
setwd(direc)

## Fix the projection and extent
projection(pr) <- mycrs
extent(pr) <- extent(plonlat)
## Take a look
pr
plot(pr)

## Project to long lat grid
or <- projectRaster(pr, crs=CRS("+proj=longlat +datum=WGS84"))
writeRaster(or, filename="or_monthlyclim.grd")
## Take a look

## 2.3 SD ##

updir <- "/disco2/aerosoles_DATA/climatologia_AOD/AOD/sd"
direc <- "/disco2/aerosoles_DATA/climatologia_AOD/AOD/bc/calc"

## List with the files I want to merge in a stack, not neccessary previous work with cdo because I have a list of files one per time step.

listFich <- dir(path="/disco2/aerosoles_DATA/climatologia_AOD/AOD/sd", pattern='\\.nc')
setwd(updir) 

## pr has the rasterStack with the monthly data of bc aerosols.
pr <- stack(listFich, varname="aero")
setwd(direc)

## Fix the projection and extent
projection(pr) <- mycrs
extent(pr) <- extent(plonlat)


## Project to long lat grid
sd <- projectRaster(pr, crs=CRS("+proj=longlat +datum=WGS84"))
writeRaster(sd, filename="sd_monthlyclim.grd")

## 2.4 SS ##

updir <- "/disco2/aerosoles_DATA/climatologia_AOD/AOD/ss"
direc <- "/disco2/aerosoles_DATA/climatologia_AOD/AOD/bc/calc"

## List with the files I want to merge in a stack, not neccessary previous work with cdo because I have a list of files one per time step.

listFich <- dir(path="/disco2/aerosoles_DATA/climatologia_AOD/AOD/ss", pattern='\\.nc')
setwd(updir) 

pr <- stack(listFich, varname="aero")
setwd(direc)

## Fix the projection and extent
projection(pr) <- mycrs
extent(pr) <- extent(plonlat)

## Project to long lat grid
ss <- projectRaster(pr, crs=CRS("+proj=longlat +datum=WGS84"))
writeRaster(ss, filename="ss_monthlyclim.grd")

## 2.5 SU ##

updir <- "/disco2/aerosoles_DATA/climatologia_AOD/AOD/su"
direc <- "/disco2/aerosoles_DATA/climatologia_AOD/AOD/bc/calc"

listFich <- dir(path="/disco2/aerosoles_DATA/climatologia_AOD/AOD/su", pattern='\\.nc')
setwd(updir) 

pr <- stack(listFich, varname="aero")
setwd(direc)

## Fix the projection and extent
projection(pr) <- mycrs
extent(pr) <- extent(plonlat)

## Project to long lat grid
su <- projectRaster(pr, crs=CRS("+proj=longlat +datum=WGS84"))
writeRaster(su, filename="su_monthlyclim.grd")
