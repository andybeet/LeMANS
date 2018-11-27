####################################################################################
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
