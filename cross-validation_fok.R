### Cross-validation approach
## author: Ying-Jung Deweese
## Date: 2017-09-10

## load package
library(raster)
library(sp)
library(dismo)
library(rgdal)
library(raster)
library(maptools)
library(Metrics)


## extract rain station point output from interpolation results

PD=readOGR(".","WY2011_vt4")
crs(PD)

# GWR

year=2011
for (i in 1:33)
{
name = sprintf("T1_%d_S%d_GWR.tif",year, i)
#name = sprintf("All_2005_S1_GWR.tiff")
GWR=raster(name)
GWRok=as(GWR,"SpatialGridDataFrame")
Bok[i]=over(PD,GWRok)

}

write.csv(Bok,file="GWR_WY2011_t1.csv")

# spBayes
for (i in 1:33)
{
  year=2011
  name = sprintf("All_%d-S%d_spBayes_t4.tif", year, i)

  spBayes=raster(name)
  
  newproj <- "+proj=lcc +lat_1=34.03333333333333 +lat_2=35.46666666666667 +lat_0=33.5
  +lon_0=-118 +x_0=2000000 +y_0=500000.0000000001 +ellps=GRS80
  +towgs84=0,0,0,0,0,0,0 +units=us-ft +no_defs"
  projection(spBayes) <- newproj

  spBayesok=as(spBayes,"SpatialGridDataFrame")
  Bok[i]=over(PD,spBayesok)

}


write.csv(Bok,file="spBayes_2011_t4.csv")


# KED
Bok=list()
year=2011
for (i in 1:30)
{
  name = sprintf("WY%d_2011-S%d_KED_t4.tiff", year, i)

  KED=raster(name)
  #projection(spBayes) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
#   newproj <- "+proj=lcc +lat_1=34.03333333333333 +lat_2=35.46666666666667 +lat_0=33.5
#   +lon_0=-118 +x_0=2000000 +y_0=500000.0000000001 +ellps=GRS80 +units=us-ft
#   +no_defs"
#   pr1 <- projectRaster(KED, crs=newproj)
  KEDok=as(KED,"SpatialGridDataFrame")
  Bok[i]=over(PD,KEDok)

}

write.csv(Bok,file="KED_WY2011_t4.csv")


####Calculate statistical metrics
## Read data

A=read.csv("WY2011_vt3.csv",header=T)
A=data.matrix(A)
A1=read.csv("spBayes_2011_t3.csv",header=T)
A1=data.matrix(A1)

## MCE calculation
s=dim(A)[2]-6
k=dim(A1)[2]
kok=dim(A1)[1]
out=matrix(NA,kok,1,byrow=T)
for (i in 1:s){
  for (j in 1:k){ 
    out=1-(sum(abs(A1[,j]-A[,i+6]))/sum(abs(A[,i]-A[,6])))  
  }  
}


## Bias calculation, ideal value is zero
s=dim(A)[2]-6
k=dim(A1)[2]
kok=dim(A1)[1]
out2=matrix(NA,kok,1,byrow=T)
for (i in 1:s){
  for (j in 1:k){ 
    out2=sum(abs(A1[,j]-A[,i+6]))/k
  }  
}

### The other statistical metrics
mae(A,A1)
accuracy(A, A1)
rmse(A,A1)
rse(A,A1)
