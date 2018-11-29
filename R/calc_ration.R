#'Calculate species ration, growth efficiency, and size class weight
#'
#' Calculates the ration, weight, and growth efficiency for a species in a size class. Ration is defined as growth increment / growth efficiency
#'
#'
#'@param nSize Number of size class intervals species can grow through
#'@param nSpecies Number of species in the model
#'@param uBound upper bound of each size class interval
#'@param lBound lower bound of each size class interval
#'@param mBound mid point of each size class interval
#'@param phiMin Model timestep (years). See \code{\link{calc_phi}}
#'
#'@return A list is returned
#'
#'    \code{ration} nSize x nSpecies matrix. (growth in time interval)/growth efficiency
#'
#'    \code{wgt} Weight of average sized fish (each size class). Uses length weight relationship. See \code{data_parameterValues}
#'
#'    \code{gEff} nSize x nSpecies matrix. Growth Efficiency of species j in size class i
#'
#'    \code{scLinf} The size class at which each species reaches L_inf (maximum length)
#'
#'
#'@section References:
#'Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63:1344-1359.
#'
#'Rochet et al. (2011). Does selective fishing conserve community biodiversity? Prediction from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68:469-486
#'
#'@export




####################################################################################
# Calculates the ration for a species in a size class
# The ration is the amount that must be consumed by a predator in size class to account for growth
# see Predation Mortality (M2) p1349 of Hall et al
# uses von bertalanfy growth equation and growth efficiency
calc_ration <- function(nSize,nSpecies,uBound,lBound,midBound,phiMin,parameterValues){
  # dumb loop. Eventually change this

  # find Size class at which each species reaches Linf
  scLinf <- sapply(parameterValues$Linf,function(x) {which((x>lBound) & (x<=uBound))})

  # find the change in length in a time interval - based on average length of fish in a bin
  wgt <- matrix(data=0,nrow=nSize,ncol=nSpecies)
  gEff <- matrix(data=0,nrow=nSize,ncol=nSpecies)
  ration <- matrix(data=0,nrow=nSize,ncol=nSpecies)
  # loop over species and their maximum size class
  for (isp in 1:nSpecies) {
    for (jsc in 1:scLinf[isp]){

      if (jsc == scLinf[isp]) {
        # the last size class contains Linf. But Linf may be reached at any point in the size class depending on the species.
        # we need to find the change in weight from thelower bound to Linf, recognizing it will exist part way in the interval
        # point at which midway from lower bound to Linf
        L1 <- lBound[scLinf[isp]]  +  (parameterValues$Linf[isp]-lBound[scLinf[isp]])/2

      } else {
        L1 <- midBound[jsc] # midpoint of interval initial Length.
      }


      W1 <- parameterValues$wa[isp] * (L1^parameterValues$wb[isp]) # initial weight
      wgt[jsc,isp] <- W1 # return this for other functions

      # change in length in unit interval of time
      deltaL <- (parameterValues$Linf[isp] - L1) * (1 - exp(-parameterValues$k[isp]*phiMin))
      L2 <- L1 + deltaL # length after time unit


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
