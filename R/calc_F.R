#'Calculate fishing Mortality
#'
#' Calculates the fishing mortality for each species in each size
#' class using a logistic selectivity curve.
#'
#'
#'@param nSize Number of size class intervals species can grow through
#'@param nSpecies Number of species in the model
#'@param mBound Vector. Mid point of each size class interval
#'@param lBound Vector. Lower bound of each size class interval
#'@param Ffull Scalar. Fishing mortaliy rate for a fully recruited fish
#'@param Falpha Scalar. Number of years to simulate.
#'@param FL50 Scalar. The length at which 50 \% selection by the fishery occurs.
#'@param scLinf Vector. The size class at which each species reaches L_inf (maximum length)
#'@param scLinfMat Binary matrix indicating which size classes each species occupies
#'@param phiMin Scalar. Model time step (years). See  \code{\link{calc_phi}}
#'@param parameterValues Matrix (nSize x nSpecies). See \code{\link{rochet_GB_parameterValues}}
#'
#'@return A matrix is returned
#'
#'    \item{F}{A matrix (size nSize x nSpecies) of F (fishing mortality) values.
#'    Note: F = 0 for size classes in which species do not grow (based on von Bertalanfy growth curve)}
#'
#'
#'@section References:
#'Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63:1344-1359.
#'
#'Rochet et al. (2011). Does selective fishing conserve community biodiversity? Prediction from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68:469-486
#'
#'@export
#'

calc_F <- function(nSize,nSpecies,mBound,lBound,Ffull,Falpha,FL50,scLinf,scLinfMat,phiMin,parameterValues) {

  L <- outer(mBound,rep(1,nSpecies))
  isFished <- outer(rep(1,nSize),parameterValues$IsFished)

  eF <- (Ffull*isFished)/(1+exp(-Falpha*(L-FL50)))

  # last size class need adjusting
  for (isp in 1:nSpecies) {
    midPoint <- lBound[scLinf[isp]]  +  (parameterValues$Linf[isp]-lBound[scLinf[isp]])/2
    eF[scLinf[isp],isp] <- (Ffull*parameterValues$IsFished[isp])/(1+exp(-Falpha*(midPoint-FL50)))
  }
  # zeros in all size classes species Linf
  eF <- eF*scLinfMat*phiMin


  return(eF)
}
