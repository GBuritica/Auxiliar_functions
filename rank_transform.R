#######################################################################
#######################################################################
#######################################################################
#######################################################################
#### Gloria Buriticá
####
#### Rank transformation
#######################################################################
r             <- function(x0){
  v0 <- abs(vector(length = length(x0), mode = "integer"))
  y0 <- order(x0, decreasing=TRUE)
  for(i in 1:length(x0))  v0[y0[i]] = i 
  return(v0)
}
ranktransform <- function(x0){
    samplerank0           <- (length(x0)+1)/r(x0)
    samplerank0[x0 == 0]  <-  0                ## zeros in the sample set to zero.
    if(sum(is.na(x0))> 0){
      samplerank0[ is.na(x0)] <- 0             ## I'm setting the NA's to 0.
      print("NA's were found and set to 0")
    }
    return(samplerank0)
}


