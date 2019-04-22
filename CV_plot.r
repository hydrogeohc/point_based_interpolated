### Cross-validation plot
## author: Ying-Jung Deweese
## Date: 2017-09-10


### load library
library(purrr)
library(ggplot2)
library(reshape)

## Read data
A=read.csv("CV_plot.csv",header=T)

## plot boxplot

boxplot(MCE~Methods,data=A,xlab=expression(bold("Method")),ylab=expression(bold("MCE")),font.axis=2,las=1,lwd=2,cex.axis = 1.2)
boxplot(Bias~Methods,data=A,xlab=expression(bold("Method")),ylab=expression(bold("Bias")),font.axis=2,las=1,lwd=2,cex.axis = 1.2)
boxplot(MAE~Methods,data=A,xlab=expression(bold("Method")),ylab=expression(bold("MAE")),font.axis=2,las=1,lwd=2,cex.axis = 1.2)
boxplot(RSE~Methods,data=A,xlab=expression(bold("Method")),ylab=expression(bold("RSE")),font.axis=2,las=1,lwd=2,cex.axis = 1.2)
boxplot(RMSE~Methods,data=A,xlab=expression(bold("Method")),ylab=expression(bold("RMSE")),font.axis=2,las=1,lwd=2,cex.axis = 1.2)

# size: 600*500


