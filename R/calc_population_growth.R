#'Calculates the change in the poulation in each time step
#'
#'Updates the species' population size. This change is a function of the number of individuals that grow out of their present size class and the number that remain.
#'
#'@param nSize Number of size class intervals species can grow through
#'@param nSpecies Number of species in the model
#'@param N A matrix (nSize x nSpecies) of species abundance (number of individuals)
#'@param probGrowOut A matrix (nSize x nSpecies) of proportions.
#'    Note: probGrowOut = 0 for size classes in which a species does not reach and probGrowOut = 0 for largest size class. It is not possible to grow out of the largest size class since it is determined by Linf. See \code{\link{calc_phi}}
#'
#'
#'@return A list is returned
#'
#'    \item{N}{A matrix (nSize x nSpecies) of abundance (number of individuals)}
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
# Calculates the numbers for each species in each
#  size class  for the next time step.

calc_population_growth <- function(nSize,nSpecies,N,probGrowOut){
  # N and probGrowout are matrices of size nsize*nSpecies
  # the population size in a class is the sum of the % that stay + % that grow out of previous class
  stay <- (1-probGrowOut)*N
  leave <- probGrowOut*N

  updatedN <- stay + rbind(rep(0,nSpecies),utils::head(leave,-1))
  return(updatedN)
}
