## ---------------------------
## Script Name: MDS_Main
## Author: Tomas Miskov
## Date Created: 2022-02-01
## Purpose: Compare own MDS package with `smacof`
## ----------------------------------------------

#--------
# SET UP |
#--------
rm(list=ls())                                         # clean the environment
if (!require("pacman")) install.packages("pacman")    # install pacman
pacman::p_load(ggplot2, tidyverse, smacof)            # pre-load packages
source("MDS_Lib.R")                                   # load local libraries

options(scipen = 6, digits = 4)                       # clean numerical notation

#------
# DATA |
#------
load("basket.Rdata")
basket <- as.matrix(basket)
rownames(basket) <- colnames(basket)

#transform data into dissimilarities using the log transformation
basketDis <- log(outer(diag(basket), diag(basket)) / (basket * t(basket)))

#------------
# BASKET MDS |
#------------
set.seed(42)
iN <- nrow(basketDis)
iP <- 2
mInit <- mInit <- matrix(runif(iN*iP, -1, 1), nrow = iN, ncol = iP) #initial X
mds1 <- mds(basketDis, init = mInit)                                #smacof
mds2 <- fnMDS(basketDis, iP = 2, mInit = mInit, bSilent = TRUE)    #ours

cat("The optimized stress from smacof is:", mds1$stress)
cat("The optimized stress from our implementation is:", mds2$stress)
cat("The difference is:", mds2$stress - mds1$stress)

par(mfrow = c(1, 2)) 
plot(mds1)
plot(mds2$conf, pch = 20, xlab = "Dim 1", ylab = "Dim 2", 
     main = "Configuration Plot", ylim = c(-3,3))
text(mds2$conf, labels=colnames(basket), cex = 0.7, pos = 3)








#--------------------------------------
# TESTS FOR INTERMEDIATE FUNCTIONALITY |
#--------------------------------------
iN <- nrow(basket) #number of items/observations taken from dissimilarity matrix
iP <- 2            #number of output dimensions
mX <- matrix(rnorm(iN*iP), nrow = iN, ncol = iP)          #initial output matrix
mX <- scale(mX, scale = FALSE)                            #column center X
mJ <- diag(iN) - iN^-1 * (rep(1, iN) %*% t(rep(1, iN)))   #centering matrix
mV <- iN * mJ                                             #Laplacian V matrix
dRawStress <- fnRawStress(basketDis, mX)                  #Raw stress score
dKruskalStress <- fnKruskalStress(basketDis, mX)          #Normalized stres score
dEtaSquared <- iN*sum(diag(t(mX) %*% mX))                 #Normalization denominator
isZero <- sqrt(dRawStress/dEtaSquared) - dKruskalStress   #Correctness check

vW <- rep(1, sum(seq(1,nrow(mX)-1)))
mW <- fnLowTriToMatrix(vW, iN)
mF <- mW * basketDis * fnEucDistMatrix(mX, FALSE)^-1
mF[is.na(mF)] <- 0
mBx <- diag(rowSums(mF)) - mF
