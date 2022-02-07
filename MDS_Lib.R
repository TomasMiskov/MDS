## ---------------------------
## Script Name: MDS_Lib
## Author: Tomas Miskov
## Date Created: 2022-02-01
## Purpose: Implement Multi-Dimensional Scaling algorithm by majorization
## ----------------------------------------------------------------------

#' Euclidean Distance Matrix Function
#' 
#' Compute distance matrix N x N of input matrix N x p
#' 
#' @param mX N x p matrix X
#' @param bLowerTriang If TRUE returns only the lower triangular as a vector
fnEucDistMatrix <- function(mX, bLowerTriang = TRUE){
  
  mDistX <- apply(mX, 1, function(x) apply(mX, 1, function(y) sqrt(sum((x-y)^2))))
  
  if(bLowerTriang){
    return(mDistX[lower.tri(mDistX)])
  } else{
    return(mDistX)
  }
}

#--------------------------------------------------------------------------
#' Lower Triangular to symmetric matrix
#' 
#' Compute symmetric matrix with 0s on the diagonal given its lower triangular
#' as a vector
#' 
#' @param vX Lower triangular as a vector
#' @param iN Dimension of the output matrix
fnLowTriToMatrix <- function(vX, iN){
  mX <- matrix(0, ncol = iN, nrow = iN)   #initialize the output matrix
  mX[lower.tri(mX)] <- vX                 #assign lower triangular
  mX <- t(mX)                             #transpose
  mX[lower.tri(mX)] <- vX                 #assign lower triangular (again)
  
  return(mX)
}

#--------------------------------------------------------------------------
#' Raw Stress Function
#' 
#' Compute raw stress of the dissimilarity matrix and the corresponding 
#' distance matrix
#' 
#' @param mDS The input dissimilarity matrix
#' @param mX N x p matrix of MDS distances between observations in X
#' @param vW Vector of weights, defaults to vector of 1s
fnRawStress <- function(mDS, mX, vW = rep(1, sum(seq(1,nrow(mX)-1)))){
  
  #extract the lower triangular as a vector
  mDS.lowTriang <- mDS[lower.tri(mDS)]
  
  #compute euclidean distances of mX
  mX.eucDist.lowTriang <- fnEucDistMatrix(mX, bLowerTriang = TRUE)
  
  #compute the raw stress score
  dRawStress <- t(vW) %*% (mDS.lowTriang - mX.eucDist.lowTriang)^2
  
  return(dRawStress)
}

#--------------------------------------------------------------------------
#' Kruskal's Stress Function
#' 
#' Compute normalized Kruskal's stress
#' 
#' @param mDS The input dissimilarity matrix
#' @param mX N x p matrix of MDS distances between observations in X
#' @param vW Vector of weights, defaults to vector of 1s
fnKruskalStress <- function(mDS, mX, vW = rep(1, sum(seq(1,nrow(mX)-1)))){
  
  #compute raw stress
  dRawStress <- fnRawStress(mDS, mX, vW)
  
  #compute distances in mX
  mX.eucDist.lowTriang <- fnEucDistMatrix(mX, bLowerTriang = TRUE)
  
  dKruskalStress <- sqrt(dRawStress / t(vW) %*% (mX.eucDist.lowTriang^2))
  
  return(dKruskalStress)
}

#--------------------------------------------------------------------------
#' Normalized Stress Function
#' 
#' Compute normalized Normalized stress
#' 
#' @param mDS The input dissimilarity matrix
#' @param mX N x p matrix of MDS distances between observations in X
#' @param vW Vector of weights, defaults to vector of 1s
fnNormalStress <- function(mDS, mX, vW = rep(1, sum(seq(1,nrow(mX)-1)))){
  
  #compute raw stress
  dRawStress <- fnRawStress(mDS, mX, vW)
  dNormalStress <- sqrt(dRawStress / t(vW) %*% (mDS[lower.tri(mDS)]^2))
  
  return(dNormalStress)
}

#--------------------------------------------------------------------------
#' MDS Majorization Function
#' 
#' Majorize the stress loss function to obtain MDS given a dissimilarity matrix
#' 
#' @param mDS The input dissimilarity matrix
#' @param vW Vector of weights, defaults to vector of 1s
#' @param iP Number of dimensions for the MDS output
#' @param dEpsilon Precision parameter
#' @param mInit Optional matrix of the initial configuration
#' @param bSilent If FALSE, the iterations are printed, default is TRUE
fnMDS <- function(mDS, vW = rep(1, sum(seq(1,nrow(mX)-1))), iP, 
                        dEpsilon = 1e-6, mInit = NULL, bSilent = TRUE){
  iN <- nrow(mDS)
  if(is.null(mInit)){                   #initialize X if none is given
    mInit <- matrix(runif(iN*iP, -1, 1), nrow = iN, ncol = iP) 
  }
  mX <- scale(mInit, scale = FALSE)     #column center initial X
  dStress <- fnNormalStress(mDS, mX, vW)   #compute initial stress value
  dStressDelta <- dStress               #temporary delta stress value
  
  mJ <- diag(iN) - iN^-1 * (rep(1, iN) %*% t(rep(1, iN)))  #centering matrix
  mV <- iN * mJ                                            #Laplacian matrix
  mGenInvV <- matrix(iN^-1, nrow = iN, ncol = iN)          #computing the generalized
  mVinv <- solve(mV + mGenInvV) - mGenInvV                 #inverse of the Laplacian
  
  k <- 1
  while((k == 1) || (dStressDelta > dEpsilon)){
    k <- k + 1
    mY <- mX
    
    #COMPUTING THE B MATRIX
    mW <- fnLowTriToMatrix(vW, iN)                 #transforming lower.tri to matrix
    mF <- mW * mDS * fnEucDistMatrix(mY, FALSE)^-1 #computing the intermediate F matrix
    mF[is.na(mF)] <- 0                             #setting NAs to 0
    mBy <- diag(rowSums(mF)) - mF                  #computing B matrix
    
    #COMPUTING UPDATE FOR X MATRIX
    mX <- mVinv %*% mBy %*% mY                     #updating X
    
    #COMPUTING NEW STRESS
    dStressY <- fnNormalStress(mDS, mY, vW)       #normalized stress iter (k-1)
    dStressX <- fnNormalStress(mDS, mX, vW)       #normalized stress iter (k)
    dStressDelta <- dStressY - dStressX            #normalized stress delta
    
    if(!bSilent){
      cat("Normalized stress value in iteration", k, "is:", dStressX)
      cat("\nStress improved from previous interation by:", dStressDelta, "\n\n")
    }
  }
  return(list("conf" = mX, "confdist" = fnEucDistMatrix(mX), 
              "stress" = fnNormalStress(mDS, mX, vW)[[1]],
              "iter" = k, "init" = mInit))
}

