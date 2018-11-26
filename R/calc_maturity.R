

####################################################################################
# Returns an nSize nSpecie smatrix indicating the proportion of each size
# class for each species that is mature and contribute to SSB.
calc_maturity <- function(nSize,nSpecies,midBound,scLinfMat,scLinf){

  # creates matrix form. Matrix operations to avoid looping
  kappaMat <- outer(rep(1,nSize),parameterValues$kappa)
  LmatMat <- outer(rep(1,nSize),parameterValues$Lmat)
  midMat <- outer(midBound,rep(1,nSpecies))

  maturity <- 1/(1+exp(-kappaMat*(midMat-LmatMat)))

  # all in laast class are mature
  for (isp in 1:nSpecies){
    maturity[scLinf[isp],isp] <- 1
  }

  maturity <- maturity*scLinfMat # multiplies by binary matrix


  return(maturity)
}

