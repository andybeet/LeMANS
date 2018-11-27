#'calculate predation Mortality
#'
#' Calculates the predation mortality for each species in each size
#' class
#'
#'
#'@param nSize Number of size class intervals species can grow through
#'@param nSpecies Number of species in the model
#'@param N an nSize x nSpecies matrix of abundance (number of individuals)
#'@param ration a list of ration, weight by class, The size class at which each species reaches L_inf (maximum length). see \code{\link{calc_sizepref_suitability}}.
#'@param suitability a 3D array of predator prey suitabilities. See \code{\link{calc_ration}}
#'@param phiMin Model timestep (years).
#'@param otherfood The amount of other food available unaccounted for in the model (grams)
#'
#'@return A matrix is returned
#'
#'    \code{M2}    - nSize x nSpecies matrix of M2 (predation mortality) values where nSpecies = number of species and nSize = number of size classes.
#'    Note that M2_i,j = 0 for size classes i in which species j is not preyed upon.
#'
#'@seealso \code{\link{calc_sizepref_suitability}}, \code{\link{calc_ration}}
#'
#'@section References:
#'Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63:1344-1359.
#'
#'Rochet et al. (2011). Does selective fishing conserve community biodiversity? Prediction from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68:469-486
#'
#'@export
#'


# 3D array
# Calculates the predation mortality for each species in each
# size class and returns an nSize x nSpecies matrix of M2 values.
calc_M2_r <- function(nSize,nSpecies,N,ration,suitability,phiMin,otherFood){

  M2 <- matrix(data=0,nrow = nSize,ncol = nSpecies)
  # mortality of prey m in size class n
  for (msp in 1:nSpecies) {  # prey species
    for (nsz in 1:ration$scLinf[msp]) { # prey size class
      for (isp in 1:nSpecies) { # predator
        for (jsz in 1:ration$scLinf[isp]) { # predator size
          index <- ((isp-1)*nSize) + jsz
          # mortality is summed over isp and jsize
          numerator <- ration$ration[jsz,isp]*N[jsz,isp]*suitability[index,nsz,msp]
          denominator <- sum(suitability[index,,] * ration$wgt * N) + otherFood
         # cat("numerator = ",numerator,"\n")
        #  cat("denominator = " ,denominator,"\n")
        #  invisible(readline("press Enter to continue ..."))
          M2[nsz,msp] <- M2[nsz,msp] + numerator/denominator
        }
      }
    }
  }
  M2 <- M2*phiMin

  return(M2)
}

# # Calculates the predation mortality for each species in each
# # size class and returns an nSize x nSpecies matrix of M2 values.
# # Modified by J. Collie on 3-Sep-09.  M2_denom is now used to
# # return the amount of prey consumed by each predator per year
# calc_M2Old <- function(nSize,nSpecies,N,ration,suitability,phiMin,otherFood){
#   # 4D array. denominator takes a lot of time to execute. ~80% of execution time
#   M2 <- matrix(data=0,nrow = nSize,ncol = nSpecies)
#   # mortality of prey m in size class n
#   for (msp in 1:nSpecies) {  # prey species
#     for (nsz in 1:ration$scLinf[msp]) { # prey size class
#       for (isp in 1:nSpecies) { # predator
#         for (jsz in 1:ration$scLinf[isp]) { # predator size
#           # mortality is summed over isp and jsize
#           numerator <- ration$ration[jsz,isp]*N[jsz,isp]*suitability[jsz,isp,nsz,msp]
#           denominator <- sum(suitability[jsz,isp,,] * ration$wgt * N) + otherFood
#           M2[nsz,msp] <- M2[nsz,msp] + numerator/denominator
#         }
#       }
#     }
#   }
#   M2 <- M2*phiMin
#
#   return(M2)
# }
