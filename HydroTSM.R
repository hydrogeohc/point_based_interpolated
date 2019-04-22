###Kriging with exteranl drift method (KED)
## Author: Ying-Jung Deweese
## Date: 2017-08-24
###load R package
library(hydroTSM)
library(zoo)
library(sp)
library(xts)
library(rgdal)
library(maptools)
library(rgeos)
library(raster)
## read file 
# read storm file information
A=read.csv("all_2005fok_0920.csv",header=T)

## read the station geographic information
year=2005
A1=read.csv("all_2005fok_0920_fok.csv",header=T)
coordinates(A1)=~x+y
## Loading the point gis datdata
B=readOGR(".","coastal_SB")

## Projection for the Subcatchments file
p4s=CRS("+proj=lcc +lat_1=34.03333333333333 +lat_2=35.46666666666667 +lat_0=33.5 +lon_0=-118 +x_0=2000000 +y_0=500000.0000000001 +ellps=GRS80 +units=us-ft +no_defs")

#Loading the DEM
DEM=raster("all_ok_Projec.tif")
ras2 <- as(DEM, "SpatialGridDataFrame")
# #Giving a meaningful name to the predictor
ras2$ELEVATION<-ras2$all_ok_Projec
# # Saving memory
ras2$all_ok_Projec<-NULL

## 4) Kriging with External Drift interpolation and plot
S=dim(A)[1]
for (i in 1:S){
x.ts=as.numeric(A[i,2:ncol(A)])
## Setting the name of the gauging stations
names(x.ts)<-colnames(A[i,2:ncol(A)])
# Computing KED
x.kedok <- hydrokrige(x.ts= x.ts, x.gis=A1, 
                    X="long", Y="lat", sname="ID",
                    p4s=p4s,
                    elevation="elv",
                    type= "cells", 
                    formula=value~ELEVATION,
                    #subcatchments= B,
                    predictors=ras2,
                    cell.size= 3000, 
                    ColorRamp= "Precipitation", 
                    main= "KED Precipitation ",
                    arrow.plot= TRUE, 
                    arrow.offset= c(6200000,5810000), arrow.scale= 20000,
                    scalebar.plot= TRUE, 
                    sb.offset= c(1960000,2030000), sb.scale= 100000)
ras2_name = sprintf("All_%d-S%d_KED.tif", year, i)
writeGDAL(x.kedok, ras2_name, drivername="GTiff", options=NULL)
}

