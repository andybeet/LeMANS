
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

# Calculates the ration for a species in a size class
# The ration is the amount that must be consumed by a predator in size class to account for growth
# see Predation Mortality (M2) p1349 of Hall et al
# uses von bertalanfy growth equation and growth efficiency
calc_ration <- function(nSize,nSpecies,uBound,lBound,midBound,phiMin){
  # dumb loop. Eventually change this

  # find Size class at which each species reaches Linf
  scLinf <- sapply(parameterValues$Linf,function(x) {which((x>lBound) & (x<=uBound))})

  wgt <- matrix(data=0,nrow=nSize,ncol=nSpecies)
  gEff <- matrix(data=0,nrow=nSize,ncol=nSpecies)
  ration <- matrix(data=0,nrow=nSize,ncol=nSpecies)
  # loop over species and theri maximum size class
  for (isp in 1:nSpecies) {
    for (jsc in 1:(scLinf[isp]-1)){

      L1 <- midBound[jsc] # initial Length
      W1 <- parameterValues$wa[isp] * (L1^parameterValues$wb[isp]) # initial weight
      wgt[jsc,isp] <- W1 # return this for other functions

      # change in length in unit interval of time
      deltaL <- (parameterValues$Linf[isp] - midBound[jsc]) * (1 - exp(-parameterValues$k[isp]*phiMin))
      L2 <- midBound[jsc] + deltaL # length after time unit
      W2 <- parameterValues$wa[isp] * (L2^parameterValues$wb[isp]) # weight after time unit
      changeInWeight <- W2-W1
      WInf <- parameterValues$wa[isp] * (parameterValues$Linf[isp]^parameterValues$wb[isp])

      # growth efficiency
      gEff[jsc,isp] <- (1-(W1/WInf)^.11)*0.5 # see paper

      ration[jsc,isp] <- changeInWeight/gEff[jsc,isp]

    }
  }


  return(list(ration=ration,scLinf=scLinf,wgt=wgt,gEff=gEff))
}
