#'Calculates recruits and spawning stock biomass (SSB)
#'
#'The number of recruits are estimated as a function of SSB using the Ricker stock-recruitment relationship.
#'The SSB is estimated as a function of the number of mature individuals.
#'
#'@param N A matrix (nSize x nSpecies) of abundance (number of individuals)
#'@param maturity A matrix (nSize x nSpecies) of the proportion of individuals that are considered mature (\code{\link{calc_maturity}})
#'@param wgt A matrix (nSize x nSpecies) of an individuals weight at the mid point of each size class (Units: grams). See \code{\link{calc_ration}}
#'@param recAlpha A vector (length nSpecies). Scale parameter of the Ricker S-R relationship.
#'@param recBeta A vector (length nSpecies). Shape parameter of the Ricker S-R relationship.

#'
#'@return A list is returned
#'
#'    \item{recruits}{A vector (length nSpecies). Number of recruits (millions).}
#'
#'    \item{SSB}{A vector (length nSpecies). SSB (1000's tonnes).}
#'
#'@section References:
#'Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63:1344-1359.
#'
#'Rochet et al. (2011). Does selective fishing conserve community biodiversity? Prediction from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68:469-486
#'
#'@export

####################################################################################
# Returns SSB and Recruits.
# Calculates the recruits from the SSB
calc_recruits <- function(N,maturity,wgt,recAlpha,recBeta){
  # this will eventually be generalized

  # SSB= proportion of mature individuals*number of individuals * their weight
  SSB <- colSums(maturity*N*wgt)  #(unit = grams)
  SSB <- SSB/1e9 #(1000's tonnes')
  # Ricker Stock recruitment curve
  recruits <- recAlpha*SSB*exp(-recBeta*SSB)
  recruits <- recruits*1e6 # number of individuals.

  return(list(recruits=recruits,SSB=SSB))
}
####################################################################################
####################################################################################
