
####################################################################################
# Calculates the numbers for each species in each
#  size class  for the next time step.

calc_population_growth <- function(nSize,nSpecies,N,probGrowOut){
  # N and probGrowout are matrices of size nsize*nSpecies
  # the population size in a class is the sum of the % that stay + % that grow out of previous class
  stay <- (1-probGrowOut)*N
  leave <- probGrowOut*N

  updatedN <- stay + rbind(rep(0,nSpecies),head(leave,-1))
  return(updatedN)
}
