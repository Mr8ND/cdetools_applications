library(rhdf5)
library(tidyverse)
library(readr)


zBeforeScale=read.table("data/z.txt",header=F,sep=",")[,1] # redshift
dataSpec=read.table("data/specs.txt",sep=",")

set.seed(401)
n=length(zBeforeScale)
zBeforeScaleNoisy=zBeforeScale+rnorm(n,0,0.02)

zMin=min(zBeforeScaleNoisy)
zMax=max(zBeforeScaleNoisy)
z=(zBeforeScaleNoisy-zMin)/(zMax-zMin)

randomPerm=sample(1:n)
z_train=as.matrix(z[randomPerm[1:2000]])
z_valid=as.matrix(z[randomPerm[2001:2400]])
z_test=as.matrix(z[randomPerm[2401:n]])

x_train=as.matrix(dataSpec[randomPerm[1:2000],])
x_valid=as.matrix(dataSpec[randomPerm[2001:2400],])
x_test=as.matrix(dataSpec[randomPerm[2401:n],])

save(x_train, x_valid, x_test,
     z_train, z_valid, z_test, file = 'data/jcgs_data_rafael.RData')
