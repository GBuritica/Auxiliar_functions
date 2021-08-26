#######################################################################
#######################################################################
#######################################################################
#######################################################################
#### Gloria Buriticá
#### 
#### Simulation observations from time series models
####
#######################################################################
#######################################################################
#######################################################################
#######################################################################
## 1 - ARMAX-Model
## 2 - Squared ARCH-Model
## 3 - ARCH-Model
## 
require(extRemes)     # For generating Frechet random variables
require(evd)
#######################################################################


#######################################################################
## Functions for creating random paths
#######################################################################
ARMAX1 <- function(lam, pathlength ,al=1){
  frech          <- rfrechet(pathlength,shape=al)
  path           <- vector(mode = "double" , pathlength)
  path[1]        <- frech[1]
  for( i in 2:pathlength) path[i] = max((lam*path[(i-1)]), (((1-lam^al)^(1/al))*frech[i]))
  return(path)
}  
#######################################################################
#######################################################################
AR1 <- function(phi, pathlength ,al=1){
  frech          <- rfrechet(pathlength,0,1,al)
  path           <- vector(mode = "double" , pathlength)
  path[1]        <- frech[1]
  for( i in 2:pathlength) path[i] = (phi*path[(i-1)]) + frech[i]
  return(path)
}  
AR2 <- function(phi2,phi1 ,pathlength  ){
  frech          <- rfrechet(pathlength,0,1,1)
  path           <- vector(mode = "double" , pathlength)
  path[1:2]        <- frech[1:2];
  for( i in 3:pathlength) path[i] = (phi1*path[(i-1)]) + (phi2*path[(i-2)])  + frech[i]
  return(path)
} 
AR5 <- function(phi2,phi1 , phi3,phi4 ,phi5, pathlength  ){
  frech          <- rfrechet(pathlength,0,1,1)
  path           <- vector(mode = "double" , pathlength)
  path[1:5]        <- frech[1:5];
  for( i in 6:pathlength) {
    path[i] <-  (phi1*path[(i-1)]) + (phi2*path[(i-2)]) +(phi3*path[(i-3)]) + (phi4*path[(i-4)]) + (phi5*path[(i-5)])+ frech[i]
  }
  return(path)
} 
#######################################################################
squaredARCH <- function(lambda, pathlength){
  path           <- vector(mode = "double" , 200000)
  norm           <- rnorm(200000,0,1)
  path[1]        <- 2*10^{-5}
  for( i in 2:200000) path[i] <-  (norm[i])^2*( (2*10^(-5)) + lambda*path[(i-1)])
  return(path[(200000-pathlength+1):200000])
}
#######################################################################
ARCHmodel <- function(lambda, pathlength){
  path           <- vector(mode = "double" , 100000)
  norm           <- rnorm(100000,0,1)
  path[1]        <- 0#2*10^{-5}   
  for( i in 2:100000) path[i] <-  (norm[i])*sqrt( (2*10^(-5)) + lambda*(path[(i-1)])^2 )
  return(path[(100000-pathlength+1):100000]) 
}
#######################################################################
GARCH <- function(pathlength){
  #######################################################################
  normal <- rnorm((pathlength+1), mean = 0, sd = 1)
  t <- rt((pathlength+1), 2)
  omega <- 10^{-6}
  alpha <- 1/4
  beta  <- 7/10
  sigma1 <- 0
  #######################################################################
  ### paths -> volatility path GARCH(1,1)
  pathsigmacuadrado1 <- 1:pathlength 
  pathsigmacuadrado1[1] <- (sigma1)^2
  ### paths -> X_t path X_t = sigma_tZ_t
  pathXGARCH <- 1:pathlength
  pathXGARCH[1] <- sigma1*normal[1]
  for(i in 2:pathlength){
   pathsigmacuadrado1[i] <- omega + alpha*(pathXGARCH[(i-1)])^2 + beta*pathsigmacuadrado1[(i-1)]
   pathXGARCH[i] <- sqrt(pathsigmacuadrado1[i])*normal[i]
  }
  return(pathXGARCH)
} 
#######################################################################
GARCH2 <- function(pathlength){
  #######################################################################
  normal <- rnorm((pathlength+1), mean = 0, sd = 1)
  t <- rt((pathlength+1), 2)
  omega <- 10^{-6}
  alpha <- 0.7
  beta  <- 0.3
  sigma1 <- 0
  #######################################################################
  ### paths -> volatility path GARCH(1,1)
  pathsigmacuadrado1 <- 1:pathlength 
  pathsigmacuadrado1[1] <- (sigma1)^2
  ### paths -> X_t path X_t = sigma_tZ_t
  pathXGARCH <- 1:pathlength
  pathXGARCH[1] <- sigma1*normal[1]
  for(i in 2:pathlength){
    pathsigmacuadrado1[i] <- omega + alpha*(pathXGARCH[(i-1)])^2 + beta*pathsigmacuadrado1[(i-1)]
    pathXGARCH[i] <- sqrt(pathsigmacuadrado1[i])*normal[i]
  }
  return(pathXGARCH)
} 
#######################################################################
logAR  <- function(pathlength){
  normal <- rnorm((pathlength+1), mean = 0, sd = 1)
  exponential <- rexp((pathlength+1), rate = 1)
  phi <- 0.5
  sigma1 <- 0
  #######################################################################
  #######################################################################
  pathXlogAR <- 1:pathlength      ## X_t <- \sigma_tZ_t
  pathXlogAR[1] <- sigma1*normal[1]
  
  AR <- 1:pathlength              ## AR <- ln(\sigma^2)
  AR[1] <- 0
  #############################
  for(i in 2:pathlength){
    AR[i] <- phi*AR[(i-1)] + exponential[i]
    pathXlogAR[i] <- exp(AR[i]/2)*normal[i]
  }
  pathsigmacuadrado2 <- exp(AR)
  return(pathXlogAR)
}
#######################################################################
pathiid <- function(lambda, pathlength){
  path           <- vector(mode = "double" , pathlength)
  for(i in 1:pathlength) path[i]<- revd(1,1,1,1)
  return(path)
}
#######################################################################

