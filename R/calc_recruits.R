

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
