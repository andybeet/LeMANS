# utility functions unlikely to change
####################################################################################
####################################################################################

# Calculates the proportion of individuals that leave
# each size class. Output is scaled relative to the
# species size class combinations which is emptied fastest
# so that this combination empties completely in one timestep.
#
# ie scaled to fasted growing individual
calc_phi <- function(nSize,nSpecies,uBound,lBound){

  #calculate the time (yrs) for an average fish to grow from the lower
  #to the upper limit of a size class (Hilborn & Walters p428, eqn 13.7.2)

  LMat <- outer(rep(1,nSize),parameterValues$Linf,"*")
  uppMat <- outer(uBound,rep(1,nSpecies),"*")
  lowMat <- outer(lBound,rep(1,nSpecies),"*")
  kMat <- outer(rep(1,nSize),parameterValues$k,"*")

  options(warn = -1) # turn warnings off. NaNs produced
  growthRate <- (1/kMat)*log((LMat-lowMat)/(LMat-uppMat))
  growthRate[is.nan(growthRate)] <- 0
  growthRate[growthRate<0] <- 0
  options(warn = 0) # turn warnings back on
  # shortest time
  phiMin <- min(growthRate[growthRate>0])
  # scale by this.
  probGrowOut <- phiMin/growthRate
  probGrowOut[is.infinite(probGrowOut)] <- 0

  return(list(phiMin=phiMin,probGrowOut=probGrowOut))
}

####################################################################################
# Calculates the ration for a species in a size class
# The ration is the amount that must be consumed by a predator in size class to account for growth
# see Predation Mortality (M2) p1349 of Hall et al
# uses von bertalanfy growth equation and growth efficiency
calc_ration <- function(nSize,nSpecies,uBound,lBound,midBound,phiMin){
  # dumb loop. Eventually change this

  # find Size class at which each species reaches Linf
  scLinf <- sapply(parameterValues$Linf,function(x) {which((x>lBound) & (x<=uBound))})

  # find the change in length in a time interval - based on average length of fish in a bin
  wgt <- matrix(data=0,nrow=nSize,ncol=nSpecies)
  gEff <- matrix(data=0,nrow=nSize,ncol=nSpecies)
  ration <- matrix(data=0,nrow=nSize,ncol=nSpecies)
  # loop over species and their maximum size class
  for (isp in 1:nSpecies) {
    for (jsc in 1:scLinf[isp]){

      if (jsc == scLinf[isp]) {
        # the last size class contains Linf. But Linf may be reached at any point in the size class depending on the species.
        # we need to find the change in weight from thelower bound to Linf, recognizing it will exist part way in the interval
        # point at which midway from lower bound to Linf
        L1 <- lBound[scLinf[isp]]  +  (parameterValues$Linf[isp]-lBound[scLinf[isp]])/2

      } else {
        L1 <- midBound[jsc] # midpoint of interval initial Length.
      }


      W1 <- parameterValues$wa[isp] * (L1^parameterValues$wb[isp]) # initial weight
      wgt[jsc,isp] <- W1 # return this for other functions

      # change in length in unit interval of time
      deltaL <- (parameterValues$Linf[isp] - L1) * (1 - exp(-parameterValues$k[isp]*phiMin))
      L2 <- L1 + deltaL # length after time unit


      W2 <- parameterValues$wa[isp] * (L2^parameterValues$wb[isp]) # weight after time unit
      changeInWeight <- W2-W1

      WInf <- parameterValues$wa[isp] * (parameterValues$Linf[isp]^parameterValues$wb[isp])

      # growth efficiency
      gEff[jsc,isp] <- (1-(W1/WInf)^.11)*0.5 # see paper

      ration[jsc,isp] <- changeInWeight/gEff[jsc,isp]
    }
  }
 return(list(ration=ration,scLinf=scLinf,wgt=wgt,gEff=gEff))
}


####################################################################################
# Returns an nSize nSpecie smatrix indicating the proportion of each size
# class for each species that is mature and contribute to SSB.
calc_maturity <- function(nSize,nSpecies,midBound,scLinfMat,scLinf){

  # creates matrix form. Matrix operations to avoid looping
  kappaMat <- outer(rep(1,nSize),parameterValues$kappa)
  LmatMat <- outer(rep(1,nSize),parameterValues$Lmat)
  midMat <- outer(midBound,rep(1,nSpecies))

  maturity <- 1/(1+exp(-kappaMat*(midMat-LmatMat)))

  # all in laast class are mature
  for (isp in 1:nSpecies){
    maturity[scLinf[isp],isp] <- 1
  }

  maturity <- maturity*scLinfMat # multiplies by binary matrix


  return(maturity)
}


####################################################################################
# calculates other sources of natural mortality other than predation
# assumption: follows a beta distribution and that the small er the size class the greater the mortality
# this is based on the ratio of the sizeclass midpoint to the largest sizeclass midpoint (0,1) variable
calc_M1 <- function(nSize,nSpecies,lBound,mBound,alphaM1,betaM1,cM1,scLinfMat,scLinf,phiMin) {

  ## need to work out why scaled by cM1
  x <- outer((mBound / max(parameterValues$Linf)),rep(1,nSpecies))
  M1 <- stats::dbeta(x,alphaM1,betaM1)*cM1

  # largest size class has wrong midpoint since Linf for a species may be < the midpoint of the bin
  # the midpoint used in midway between lowerbounfd and Linf
  for (isp in 1:nSpecies) {
    midPoint <- lBound[scLinf[isp]]  +  (parameterValues$Linf[isp]-lBound[scLinf[isp]])/2
    x <- midPoint/max(parameterValues$Linf)
    M1[scLinf[isp],isp] <- stats::dbeta(x,alphaM1,betaM1)*cM1
  }
  # scaled to time step
  M1 <- M1*scLinfMat*phiMin

  return(M1)
}


####################################################################################
# Calculates the lognormal probability functions for
# prey preferences, based on the predator/prey size(wgt) ratio.
# Returns a 4D matrix with the (prey/predator) body size ratio
# for predator of size i, species j, and prey size k species l.
#  Modified by J. Collie on 17-juin-09 to omit the standardization
# then calculates standardized suitability based on foodweb matirx
calc_sizePref_suitability <- function(nSize,nSpecies,mBound,spMu,spSigma,wgt,scLinf,FW) {
  # should vectorize operations. # do later

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


####################################################################################
# Calculates the fishing mortality for each species in each size
# class using a logistic selectivity curve.
# Returns the nSize x nSpecies matrix of F values.
calc_F <- function(nSize,nSpecies,mBound,lBound,Ffull,Falpha,FL50,scLinf,scLinfMat,phiMin) {

  L <- outer(mBound,rep(1,nSpecies))
  isFished <- outer(rep(1,nSize),parameterValues$IsFished)

  eF <- (Ffull*isFished)/(1+exp(-Falpha*(L-FL50)))

  # last size class need adjusting
  for (isp in 1:nSpecies) {
    midPoint <- lBound[scLinf[isp]]  +  (parameterValues$Linf[isp]-lBound[scLinf[isp]])/2
    eF[scLinf[isp],isp] <- (Ffull*parameterValues$IsFished[isp])/(1+exp(-Falpha*(midPoint-FL50)))
  }
  # zeros in all size classes  species Linf
  eF <- eF*scLinfMat*phiMin


  return(eF)
}


####################################################################################
# Calculates the predation mortality for each species in each
# size class and returns an nSize x nSpecies matrix of M2 values.
# Modified by J. Collie on 3-Sep-09.  M2_denom is now used to
# return the amount of prey consumed by each predator per year
calc_M2 <- function(nSize,nSpecies,N,ration,suitability,phiMin,otherFood){

  M2 <- matrix(data=0,nrow = nSize,ncol = nSpecies)
  # mortality of prey m in size class n
  for (msp in 1:nSpecies) {  # prey species
    for (nsz in 1:ration$scLinf[msp]) { # prey size class
      for (isp in 1:nSpecies) { # predator
        for (jsz in 1:ration$scLinf[isp]) { # predator size
          # mortality is summed over isp and jsize
          numerator <- ration$ration[jsz,isp]*N[jsz,isp]*suitability[jsz,isp,nsz,msp]
          denominator <- sum(suitability[jsz,isp,,] * ration$wgt * N) + otherFood
          M2[nsz,msp] <- M2[nsz,msp] + numerator/denominator
        }
      }
    }
  }
  M2 <- M2*phiMin

  # place holder for M2_denom in case need to implement this
  M2_denom=NULL

  return(list(M2=M2,M2_denom=M2_denom))
}


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
