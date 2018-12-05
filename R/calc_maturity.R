#'Calculates Maturity
#'
#'The proportion of each species (in each size class) that are considered mature.
#'Matrue indisviduals contribute to SSB.
#'
#'@param nSize Number of size class intervals species can grow through
#'@param nSpecies Number of species in the model
#'@param mBound Mid point of each size class interval
#'@param scLinfMat A matrix (nSize x nSpecies) of binary values indicating which size classes each species occupies.
#'@param scLinf A vector (length nSpecies) of size classes indicating the maximum size class a species can occupy. Based on L_inf (maximum length)
#'@param parameterValues A Matrix of species specific parameters. See \code{\link{rochet_GB_parameterValues}}
#'
#'@return A matrix is returned
#'
#'    \item{Maturity}{A matrix (nSize x nSpecies) of Maturity proportions.
#'    Note: Maturity = 0 for size classes in which a species does not reach.}
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
# Returns an nSize nSpecie smatrix indicating the proportion of each size
# class for each species that is mature and contribute to SSB.
calc_maturity <- function(nSize,nSpecies,mBound,scLinfMat,scLinf,parameterValues){

  # creates matrix form. Matrix operations to avoid looping
  kappaMat <- outer(rep(1,nSize),parameterValues$kappa)
  LmatMat <- outer(rep(1,nSize),parameterValues$Lmat)
  midMat <- outer(mBound,rep(1,nSpecies))

  maturity <- 1/(1+exp(-kappaMat*(midMat-LmatMat)))

  # all in laast class are mature
  for (isp in 1:nSpecies){
    maturity[scLinf[isp],isp] <- 1
  }

  maturity <- maturity*scLinfMat # multiplies by binary matrix


  return(maturity)
}

