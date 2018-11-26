
####################################################################################
# Calculates the lognormal probability functions for
# prey preferences, based on the predator/prey size(wgt) ratio.
# Returns a 3D matrix with the (prey/predator) body size ratio
# for predator of size i, species j, and prey size k species l.
# then calculates standardized suitability based on foodweb matirx
calc_sizepref_suitability <- function(nSize,nSpecies,mBound,spMu,spSigma,wgt,scLinf,FW) {
  # should vectorize operations. # do later

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
# old version with 4d arrays
# Calculates the lognormal probability functions for
# prey preferences, based on the predator/prey size(wgt) ratio.
# Returns a 3D matrix with the (prey/predator) body size ratio
# for predator of size i, species j, and prey size k species l.
#  Modified by J. Collie on 17-juin-09 to omit the standardization
# then calculates standardized suitability based on foodweb matirx
calc_sizepref_suitabilityOld <- function(nSize,nSpecies,mBound,spMu,spSigma,wgt,scLinf,FW) {

  # 4 dimensional array
  sizePref <- array(data=0,dim=c(nSize,nSpecies,nSize,nSpecies))
  suitability <- array(data=0,dim=c(nSize,nSpecies,nSize,nSpecies))

  for (isp in 1:nSpecies) { # predator
    for (jsp in 1:nSpecies) { # prey.  pair of species to calculate ratios

      for (isize in 1:scLinf[isp]) { # predator
        for (jsize in 1:min(scLinf[isp],scLinf[jsp])) { # can only eat at least as big as itself
          ratio <- wgt[jsize,jsp]/wgt[isize,isp]
          sizePref[isize,isp,jsize,jsp] <- dlnorm(ratio,spMu,spSigma)
          suitability[isize,isp,jsize,jsp] <- sizePref[isize,isp,jsize,jsp]*FW[isp,jsp]
        }
      }
    }
  }

  # standardize the suitabilirties so sum to 1 (Hall et al reference to Magnuson (multispecies vpa) 1995)
  for (isp in 1:nSpecies) {
    for (isize in 1:scLinf[isp]) {
      standardize <- sum(suitability[isize,isp,,])
      if (standardize > 0) { # species is a predator of something
        suitability[isize,isp,,] <- suitability[isize,isp,,]/standardize
      }
    }
  }

  return(list(sizePref=sizePref,suitability=suitability))
}
