###Spatial Bayesian regression
## Author: Ying-Jung Deweese
## Date: 2017-08-24
## Load the libraries to be used
library(coda)
library(abind)
library(magic)
library(maps)
library(Matrix)
library(Formula)
library(spBayes)
library(MBA)
library(geoR)
library(spam)
library(grid)
library(fields)
library(sp)
library(maptools)
library(rgdal)
library(classInt)
library(lattice)
library(raster)

### Read file 
year = 2005
filename = sprintf("all_2005fok_0920.csv", year);

A=read.csv(filename,header=T)

## Data preliminaries
## Read data length
S_num = dim(A)[2]-7


## Compute the spatial Bayesian regression
for (i in 1:S_num)
{

print(sprintf("start of loop %d", i))

## Subset data variable for elevation, rain variables
bio <- A[,i+7]
elv <- A[,2]


## Extract the coordinates
coords <- as.matrix(A[,c("long","lat")])


p <- 2 ## This is the number of columns in the design matrix
## Set the prior mean and precision for the regression
beta.prior.mean <- as.matrix(rep(0, times=p))
beta.prior.precision <- matrix(0, nrow=p, ncol=p)

## For use with bayesGeostatExact, do the following
phi <- 0.5 ## Set the spatial range (from the variogram)
alpha <- 1600/0.5 ## Set the nugget/partial-sill ratio
sigma.sq.prior.shape <- 2.0 ## Set IG shape for sigma.sq (partial sill)
sigma.sq.prior.rate <- 0.5 ## Set IG scale for sigma.sq (partial sill)

## Run bayesGeostatExact to deliver exact posterior samples
sp.exact <- bayesGeostatExact(
bio~elv,
coords=coords, n.samples=1000,
beta.prior.mean=beta.prior.mean,
beta.prior.precision=beta.prior.precision,
cov.model="spherical",
phi=phi, alpha=alpha,
sigma.sq.prior.shape=sigma.sq.prior.shape,
sigma.sq.prior.rate=sigma.sq.prior.rate,
sp.effects=FALSE)

##Produce the posterior summaries
round(summary(sp.exact$p.samples)$quantiles,3)

cook=coords
## Run spLM to deliver MCMC samples from marginal posterior distributions
n.samples <- 1000
bef.sp <- spLM(bio~elv,
coords=cook, starting=list("phi"=3/200,"sigma.sq"=0.08,

                           "tau.sq"=0.02), tuning=list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0.05),
               

priors=list("phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(2, 0.08),"tau.sq.IG"=c(2, 0.02)), cov.model="spherical",n.samples=n.samples)

round(summary(mcmc(bef.sp$p.theta.samples))$quantiles,3)

## Recover spatial residuals using spRecover
burn.in <- floor(0.75*n.samples)
bef.sp <- spRecover(bef.sp, start=burn.in, thin=2)

## The posterior samples of the regression coefficients and the spatial effects can then be obtained as
beta.samples = bef.sp$p.beta.recover.samples
w.samples = bef.sp$p.w.recover.samples

## Obtain trace plots for regression coefficients
par(mfrow=c(3,2))
plot(beta.samples, auto.layout=TRUE, density=FALSE)

## Obtain posterior means and sd's of spatial residuals for each location
w.hat.mu <- apply(w.samples,1,mean)
w.hat.sd <- apply(w.samples,1,sd)

## Obtain OLS residuals
lm.test = lm(bio~elv)
test.resid = resid(lm.test)

## Plot the spatial residual mean surface and a map of sd's
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords, test.resid), no.X=100, no.Y=100, extend=FALSE)$xyz.est
z.lim <- range(surf[[3]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", zlim=z.lim, main="LM residuals")
surf <- mba.surf(cbind(coords, w.hat.mu), no.X=100, no.Y=100, extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", zlim=z.lim, main="Mean spatial effects")

## Predictions
BEF.shp <- readShapePoly("coastalbasinsdata.shp")
shp2poly <- BEF.shp@polygons[[1]]@Polygons[[1]]@coords
BEF.poly <- as.matrix(shp2poly)
BEF.grids <- readGDAL("all_ok_Projec.tif")


## Construct the prediction design matrix for the entire grid extent.
pred.covars <- cbind(BEF.grids[["band1"]], BEF.grids[["band2"]], BEF.grids[["band3"]], BEF.grids[["band4"]], BEF.grids[["band5"]])
pred.covars <- cbind(rep(1, nrow(pred.covars)), pred.covars)
p=SpatialPointsDataFrame(coordinates(BEF.grids),data=BEF.grids@data) # convert spatialgrid dataframe to spatial point dataframe

## Extract the coordinates of the BEF bounding polygon vertices and use the pointsInPoly (spBayes) function to obtain the desired subset of the prediction design matrix and associated prediction coordinates (i.e., pixel centroids).
pred.coords <- SpatialPoints(BEF.grids)@coords
pointsInPolyOut <- pointsInPoly(BEF.poly, pred.coords)
pred.coords <- pred.coords[pointsInPolyOut,]
pred.covars <- p[pointsInPolyOut,]
#convert spdf to matrix
DF <- as.matrix(pred.covars@data)
DF <- cbind(rep(1, nrow(DF)), DF)# generate 2 column
bef.bio.pred <- spPredict(bef.sp, start=burn.in, thin=2, pred.coords=pred.coords, pred.covars=DF)

## Mapping the predicted values
bef.bio.pred.mu = apply(bef.bio.pred$p.y.predictive.samples,1,mean)
bef.bio.pred.sd = apply(bef.bio.pred$p.y.predictive.samples,1,sd)
surf <- mba.surf(cbind(coords,bio), no.X=500, no.Y=500, extend=TRUE, sp=TRUE)$xyz.est


## write out raster

ras3_name = sprintf("All_%d-S%d_spBayes.tif", year, i)
writeGDAL(surf, ras3_name, drivername="GTiff", options=NULL)

}
