#'Calculates M1 mortality
#'
#'calculates other sources of natural mortality other than predation.
#'It is assumed that M1 mortality follows a beta distribution and that in general the small and large size classses suffer the greatest mortality.
#'
#'@param nSize Number of size class intervals species can grow through
#'@param nSpecies Number of species in the model
#'@param lBound Lower bound of each size class interval
#'@param mBound Mid point of each size class interval
#'@param alphaM1 Shape parameter (alpha) of the beta distribution
#'@param betaM1 Shape parameter (beta) of the beta distribution
#'@param cM1 Scales the pdf ??????
#'@param scLinf The size class at which each species reaches L_inf (maximum length)
#'@param phiMin Model timestep (years). See \code{\link{calc_phi}}
#'@param parameterValues Matrix of species specific parameters. See \code{\link{rochet_GB_parameterValues}}

#'
#'@return A matrix is returned
#'
#'    \item{M1}{A matrix (size nSize x nSpecies) of M1 (natural mortality) values.
#'    Note: M1 = 0 for size classes in which a species does not reach/grow in to.}
#'
#'
#'@section References:
#'Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63:1344-1359.
#'
#'Rochet et al. (2011). Does selective fishing conserve community biodiversity? Prediction from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68:469-486
#'
#'@export
#'


####################################################################################
# calculates other sources of natural mortality other than predation
# assumption: follows a beta distribution and that the smaller the size class the greater the mortality
# this is based on the ratio of the sizeclass midpoint to the largest sizeclass midpoint (0,1) variable
calc_M1 <- function(nSize,nSpecies,lBound,mBound,alphaM1,betaM1,cM1,scLinfMat,scLinf,phiMin,parameterValues) {

  ## need to work out why scaled by cM1
  x <- outer((mBound / max(parameterValues$Linf)),rep(1,nSpecies))
  M1 <- stats::dbeta(x,alphaM1,betaM1)*cM1

  # largest size class has wrong midpoint since Linf for a species may be < the midpoint of the bin
  # the midpoint used in midway between lowerbounfd and Linf
  for (isp in 1:nSpecies) {
    midPoint <- lBound[scLinf[isp]]  +  (parameterValues$Linf[isp]-lBound[scLinf[isp]])/2
    x <- midPoint/max(parameterValues$Linf)
    M1[scLinf[isp],isp] <- stats::dbeta(x,alphaM1,betaM1)*cM1
  }
  # scaled to time step
  M1 <- M1*scLinfMat*phiMin

  return(M1)
}
