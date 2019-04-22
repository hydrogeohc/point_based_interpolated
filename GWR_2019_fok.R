#### Geographically Weighted Regression (GWR)
### Author: Ying-Jung Deweese
### Date:2019-03-24

## Load library
library(spgwr)
library(sf)
library(sp)
library(spData)
library(Rcpp)
library(spdep)
library(ggplot2)
library(rgdal)
library(lattice)
library(maptools)
library(raster)
library(robustbase)
library(GWmodel)
library(RColorBrewer)


### Read the rainfall information
P=read.csv("allok_WY2005_fok.csv",header=T)

### Read Digital Elevation Model(DEM)
DEM <- raster("all_ok.tif")

##check the file dimension
S=dim(P)[2]-7
year=2005

## generate the output matrix
out=matrix(NA,S,1,byrow=T)

## generate spatialpoint dataframe
locations=cbind(P$lon,P$lat)
P.spdf=SpatialPointsDataFrame(P[,3:4],P)

##Compute basic GWR model and export the output

for (i in 1:length(S)){
    S_index = i + 7
    best.bw =gwr.sel(P[,S_index]~elv,data=P,coords=locations)
    gwr.res=gwr.basic(P[,S_index]~elv,data=P.spdf,bw=best.bw,kernel='gaussian')
    ### estimated coeffient
    elv.coefs <- gwr.res$SDF$elv
    int.coefs=gwr.res$SDF$Intercept
    ME=median(elv.coefs,)
    MEIN=median(int.coefs,)
    ### Apply 10m DEM to the regression model results
    rain_day=DEM*ME+MEIN
    ras2 <- as(rain_day, "SpatialGridDataFrame")
    ras2_name = sprintf("All_%d_S%d_GWR.tiff", year, i)
    writeGDAL(ras2,ras2_name, drivername="GTiff", options=NULL)
}
