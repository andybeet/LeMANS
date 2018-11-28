#'calculate size preference and suitability
#'
#'Calculates the lognormal probability functions for prey preferences, based on the predator/prey size (wgt) ratio.
#'Returns a 3D matrix with the (prey/predator) body size ratio for predator of size i, species j, and prey size k species l.
# then calculates standardized suitability based on foodweb matrix
#'
#'
#'@param nSize Number of size class intervals species can grow through
#'@param nSpecies Number of species in the model
#'@param mBound Mid point of each size class interval
#'@param spMu Mean of log normal distribution
#'@param spSigma Standard deviation of the log normal distribution
#'@param wgt Weight of species at the mid point of each size class (Units: grams). See \code{\link{calc_ration}}
#'@param scLinf The size class at which each species reaches L_inf (maximum length)
#'@param FW Food web represented as a binary matrix. Predator (columns) and prey (rows). See \code{\link{data_foodweb}}
#'
#'@return A list is returned
#'
#'    \code{sizePref} 3D array of predator size preference of prey size.
#'
#'    \code{suitability} 3D array of predator size preference of prey size. (sizePref*Food Web)
#'
#'
#'@section References:
#'Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63:1344-1359.
#'
#'Rochet et al. (2011). Does selective fishing conserve community biodiversity? Prediction from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68:469-486
#'
#'@export

####################################################################################
# Calculates the lognormal probability functions for
# prey preferences, based on the predator/prey size(wgt) ratio.
# Returns a 3D matrix with the (prey/predator) body size ratio
# for predator of size i, species j, and prey size k species l.
# then calculates standardized suitability based on foodweb matirx
calc_sizepref_suitability <- function(nSize,nSpecies,mBound,spMu,spSigma,wgt,scLinf,FW) {
  # should vectorize operations. # do later
  # transpose foodweb. predator on rows, prey columns
  FW <- t(FW)
  # 3 dimensional array to utilize Rcpprmadillo's speed
  sizePref <- array(data=0,dim=c(nSize*nSpecies,nSize,nSpecies))
  suitability <- array(data=0,dim=c(nSize*nSpecies,nSize,nSpecies))

  for (isp in 1:nSpecies) { # predator
    for (isize in 1:scLinf[isp]) { # predator
      index <- ((isp-1)*nSize)  + isize

      for (jsp in 1:nSpecies) { # prey.  pair of species to calculate ratios
        for (jsize in 1:min(scLinf[isp],scLinf[jsp])) { # can only eat at least as big as itself
          ratio <- wgt[jsize,jsp]/wgt[isize,isp]
          sizePref[index,jsize,jsp] <- stats::dlnorm(ratio,spMu,spSigma)
          suitability[index,jsize,jsp] <- sizePref[index,jsize,jsp]*FW[isp,jsp]
        }
      }

    }
  }

  # standardize the suitabilirties so sum to 1 (Hall et al reference to Magnuson (multispecies vpa) 1995)
  for (isp in 1:nSpecies) {
    for (isize in 1:scLinf[isp]) {
      index <- ((isp-1)*nSize) + isize
      standardize <- sum(suitability[index,,])
      if (standardize > 0) { # species is a predator of something
        suitability[index,,] <- suitability[index,,]/standardize
      }
    }
  }
  return(list(sizePref=sizePref,suitability=suitability))
}

####################################################################################
# # old version with 4d arrays
# # Calculates the lognormal probability functions for
# # prey preferences, based on the predator/prey size(wgt) ratio.
# # Returns a 3D matrix with the (prey/predator) body size ratio
# # for predator of size i, species j, and prey size k species l.
# #  Modified by J. Collie on 17-juin-09 to omit the standardization
# # then calculates standardized suitability based on foodweb matirx
# calc_sizepref_suitabilityOld <- function(nSize,nSpecies,mBound,spMu,spSigma,wgt,scLinf,FW) {
#
#   # 4 dimensional array
#   sizePref <- array(data=0,dim=c(nSize,nSpecies,nSize,nSpecies))
#   suitability <- array(data=0,dim=c(nSize,nSpecies,nSize,nSpecies))
#
#   for (isp in 1:nSpecies) { # predator
#     for (jsp in 1:nSpecies) { # prey.  pair of species to calculate ratios
#
#       for (isize in 1:scLinf[isp]) { # predator
#         for (jsize in 1:min(scLinf[isp],scLinf[jsp])) { # can only eat at least as big as itself
#           ratio <- wgt[jsize,jsp]/wgt[isize,isp]
#           sizePref[isize,isp,jsize,jsp] <- dlnorm(ratio,spMu,spSigma)
#           suitability[isize,isp,jsize,jsp] <- sizePref[isize,isp,jsize,jsp]*FW[isp,jsp]
#         }
#       }
#     }
#   }
#
#   # standardize the suitabilirties so sum to 1 (Hall et al reference to Magnuson (multispecies vpa) 1995)
#   for (isp in 1:nSpecies) {
#     for (isize in 1:scLinf[isp]) {
#       standardize <- sum(suitability[isize,isp,,])
#       if (standardize > 0) { # species is a predator of something
#         suitability[isize,isp,,] <- suitability[isize,isp,,]/standardize
#       }
#     }
#   }
#
#   return(list(sizePref=sizePref,suitability=suitability))
# }
