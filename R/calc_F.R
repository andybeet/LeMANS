####################################################################################
# Calculates the fishing mortality for each species in each size
# class using a logistic selectivity curve.
# Returns the nSize x nSpecies matrix of F values.
calc_F <- function(nSize,nSpecies,mBound,lBound,Ffull,Falpha,FL50,scLinf,scLinfMat,phiMin) {

  L <- outer(mBound,rep(1,nSpecies))
  isFished <- outer(rep(1,nSize),parameterValues$IsFished)

  eF <- (Ffull*isFished)/(1+exp(-Falpha*(L-FL50)))

  # last size class need adjusting
  for (isp in 1:nSpecies) {
    midPoint <- lBound[scLinf[isp]]  +  (parameterValues$Linf[isp]-lBound[scLinf[isp]])/2
    eF[scLinf[isp],isp] <- (Ffull*parameterValues$IsFished[isp])/(1+exp(-Falpha*(midPoint-FL50)))
  }
  # zeros in all size classes  species Linf
  eF <- eF*scLinfMat*phiMin


  return(eF)
}
