#'Calculates the proportion of individuals that leave a size class
#'
#'The proportion of individuals that leave each size class. These proportion are scaled to the fastest growing species/size class combination.
#'This ensures that no species can grow though more than one class interval in any time step.
#'This is based on the time (yrs) for an average fish to grow from the lower limit of a size class to the
#'upper limit of a size class (Hilborn & Walters p428, eqn 13.7.2).
#'
#'@param nSize Number of size class intervals species can grow through
#'@param nSpecies Number of species in the model
#'@param uBound Upper bound of each size class interval
#'@param lBound Lower bound of each size class interval
#'@param parameterValues A Matrix of species specific parameters. See \code{\link{rochet_GB_parameterValues}}


#'
#'@return A list is returned
#'
#'    \item{probGrowOut}{A matrix of size nSize x nSpecies of proportions. A cell represent the proportion of individuals growing into the following cell.
#'    Note: probGrowOut = 0 for size classes in which a species does not reach (due to its L_inf). probGrowOut = 0 for largest size class for each species}
#'
#'    \item{phiMin}{Model timestep (years). Shhortest time taken by any species to grow into next size class. probGrowOut is normalized by this quantity)}
#'
#'@section References:
#'Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63:1344-1359.
#'
#'Rochet et al. (2011). Does selective fishing conserve community biodiversity? Prediction from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68:469-486
#'
#'@export


# Calculates the proportion of individuals that leave
# each size class. Output is scaled relative to the
# species size class combinations which is emptied fastest
# so that this combination empties completely in one timestep.
#
# ie scaled to fasted growing individual
calc_phi <- function(nSize,nSpecies,uBound,lBound,parameterValues){

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
  # scale to the fastest fish and size class combination to calculate proportion leaving a size class in a time step.
  probGrowOut <- phiMin/growthRate
  probGrowOut[is.infinite(probGrowOut)] <- 0

  return(list(phiMin=phiMin,probGrowOut=probGrowOut))
}
