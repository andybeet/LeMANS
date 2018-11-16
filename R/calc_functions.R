
# utility functions unlikely to change


# Calculates the proportion of individuals that leave
# each size class. Output is scaled relative to the
# species size class combinations which is emptied fastest
# so that this combination empties completely in one timestep.
#
# ie scaled to fasted growing individual
calc_phi <- function(nSize,nSpecies,uBound,lBound){

  #calculate the time (yrs) for an average fish to grow from the lower
  #to the upper limit of a size class (Hilborn & Walters p428, eqn 13.7.2)

  LMat <- outer(rep(1,nSize),parameterValues$Linf,"*")
  uppMat <- outer(uBound,rep(1,nSpecies),"*")
  lowMat <- outer(lBound,rep(1,nSpecies),"*")
  kMat <- outer(rep(1,nSize),parameterValues$k,"*")

  options(warn = -1) # turn warnings off. NaNs produced
  growthRate <- (1/kMat)*log((LMat-lowMat)/(LMat-uppMat))
  growthRate[is.nan(growthRate)] <- 0
  growthRate[growthRate<0] <- 0
  options(warn = 0) # turn warnings back on
  # shortest time
  phiMin <- min(growthRate[growthRate>0])
  # scale by this.
  probGrowOut <- phiMin/growthRate
  probGrowOut[is.infinite(probGrowOut)] <- 0

  return(list(phiMin=phiMin,probGrowOut=probGrowOut))
}

