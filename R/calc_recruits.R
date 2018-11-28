#'Calculates recruits and spawning stock biomass (SSB)
#'
#'The number of recruits are estimated as a function of SSB using the Ricker stock-recruitment relationship.
#'The SSB is estimated as a function of the number of mature individuals.
#'
#'@param N nSize x nSpecies matrix of abundance (number of individuals)
#'@param maturity nSize x nSpecies matrix of the proportion of each size class for each species that is mature (\code{\link{calc_maturity}})
#'@param wgt Weight of species at the mid point of each size class (Units: grams). See \code{\link{calc_ration}}
#'@param recAlpha Scale parameter of the Ricker S-R relationship. A vector of length nSpecies.
#'@param recBeta Shape parameter of the Ricker S-R relationship. A vector of length nSpecies.

#'
#'@return A list is returned
#'
#'    \code{recruits}    Number of recruits (millions). A vector of length nSpecies.
#'
#'    \code{SSB}  SSB (1000's tonnes). A vector of length nSpecies.
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
