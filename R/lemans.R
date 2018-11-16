#'LeMANS program
#'
#'Runs the LeMANS model given the users specification for data need
#'
#'
#'
#'
#'
#'



# Note: all data fines are read in to memory when package is loaded
# initialValues, parameterValues, species, foodweb
lemans <- function() {
  # initial set up
  numSizeClass <- dim(initialValues)[2]
  numSpecies <- dim(initialValues)[1]
  otherFood <- 55000000 # (grams)
  convertCatch <- 3046527 # catch per tow to total numbers
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
  lowerSizeClassInterval <- seq(from=0,to=maxFishSize-maxFishSize/numSizeClass,length.out = numSizeClass)
  upperSizeClassInterval <- lowerSizeClassInterval + maxFishSize/numSizeClass
  midSizeClassInterval <- lowerSizeClassInterval + (upperSizeClassInterval-lowerSizeClassInterval)/2






}
