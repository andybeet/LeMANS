#'LeMANS program
#'
#'Runs the LeMANS model given the users specification for data need
#'
#'
#'
#'
#'
#'@export



# Note: all data fines are read in to memory when package is loaded
# initialValues, parameterValues, species, foodweb
lemans <- function() {
  # initial set up
  nSizeClass <- dim(initialValues)[2]
  nSpecies <- dim(initialValues)[1]
  otherFood <- 55000000 # (grams)
  convertCatch <- 3046527 # number of tows: converts catch per tow to total numbers
  predationFlag <- 1 # turns off/on predation
  # Fishing Parameters
  Falpha <- 0.25 # Steepness of selctivity curve
  FL50 <- 25 # Length at 50% selection
  #Model Parameters
  alphaInt <- 11 # scaling for Ricker alpha wrt L_inf
  alphaExp <- -2.298 # scaling for Ricker alpha wrt L_inf
  betaInt <- 0.1513 # scaling for Ricker beta wrt Smax
  betaExp <- 0.9484 # scaling for Ricker beta wrt Smax
  SmaxScale <- 1 # Catchability scaling for S.max estimate from survery
  # M1
  alphaM1 <- 0.8 # parametes of beta function used to model M1
  betaM1 <- 0.4 # parametes of beta function used to model M1
  cM1 <- 0.35 #scaling of final M1
  # Calculate upper and lower size class bins
  maxFishSize <- max(parameterValues$Linf) * 1.001
  lowScBin <- seq(from=0,to=maxFishSize-maxFishSize/nSizeClass,length.out = nSizeClass)
  uppScBin <- lowScBin + maxFishSize/nSizeClass
  midScBin <- lowScBin + (uppScBin-lowScBin)/2
  # transpose foodweb. predator on rows, prey columns
  FW <- t(foodweb)*predationFlag
  N <- t(initialValues)
  N <- N*convertCatch
  # size prefernce function Parameters
  spMu <- 0.5
  spSigma <- 2

  #logical matrix reperesnteing size class bins applicable for each species
  scLinfMat <- sapply(parameterValues$Linf,function(x) {x>lowScBin})

  ##################################################################
  ##################################################################
  ##################################################################
  # Now calculate all parts required to run model
  ##################################################################
  ##################################################################
  ##################################################################

  ##################################################################
  # calculate the proportion leaving each size class per time step
  phi <- calc_phi(nSizeClass,nSpecies,uppScBin,lowScBin)
  ##################################################################
  # calculate the ration.
  ration <- calc_ration(nSizeClass,nSpecies,uppScBin,lowScBin,midScBin,phi$phiMin)
  ##################################################################
  # calculate maturity
  mature <- calc_maturity(nSizeClass,nSpecies,midScBin,scLinfMat)
  ##################################################################
  # calculate recruitment. Assumes all has a ricker form. All are scaled.
  # see Hall et al paper
  recruitAlphas <- exp(alphaInt - abs(alphaExp*log(parameterValues$Linf)))
  recruitBetas <- exp(betaInt - betaExp*log(parameterValues$Smax*SmaxScale))
  recruitAlphas[1] <- 400 # trial for forage fish
  ##################################################################
  # calculate M1 (residual natural mortality)
  M1 <- calc_M1(nSizeClass,nSpecies,lowScBin,midScBin,alphaM1,betaM1,cM1,scLinfMat,ration$scLinf,phi$phiMin)
  ##################################################################
  # calculated the size preference and the suitabilities
  M2PrefSuit <- calc_sizePref_suitability(nSizeClass,nSpecies,midScBin,spMu,spSigma,ration$wgt,ration$scLinf,FW)
  ##################################################################
  # calculates the fishing mortalities

  ##################################################################
  # calculates the predation mortalities
  ##################################################################
  ##################################################################
  ##################################################################
  ##################################################################


  ##################################################################
  ##################################################################
  ##################################################################
  # Now run the model
  ##################################################################
  ##################################################################
  ##################################################################


  return(M2PrefSuit)

}
