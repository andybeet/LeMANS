#'Key_Run (Length-based Multispecies Analysis by numerical Simulation)
#'
#'Runs the LeMANS model. currently loads Rochet et al data when package loads.
#'The code was adapted directly from MATLAB code used for the
#'Hall et al and Rochet et al papers.
#'
#'@param Ffull Fishing mortaliy rate for a fully recruited fish
#'@param nYrs Number of years to simulate.
#'@param modelSetupData list of parameters required for model (See \code{\link{data_modelSetup}})
#'@param parameterValues matrix of species specific parameters (See \code{\link{data_parameterValues}})
#'@param initialValues matrix of initial abundance estimates (See \code{\link{data_initialValues}})
#'@param foodweb predator prey food web (See \code{\link{data_foodweb}})
#'@param species matrix of species names and guild membership (See \code{\link{data_species}})
#'
#'@seealso \code{\link{data_foodweb}}, \code{\link{data_initialValues}}, \code{\link{data_parameterValues}}, \code{\link{data_species}}
#'@export



# Note: all data fines are read in to memory when package is loaded
# initialValues, parameterValues, species, foodweb
key_run <- function(Ffull,nYrs,modelSetupData,parameterValues,initialValues,foodweb,species) {
  start <- Sys.time()
  # initial set up
  nSizeClass <- dim(initialValues)[2]
  nSpecies <- dim(initialValues)[1]
  ##################################################################
  # calculate recruitment. Assumes all has a ricker form. All are scaled.
  # see Hall et al paper
  alphaInt <- modelSetupData$alphaInt # scaling for Ricker alpha wrt L_inf
  alphaExp <- modelSetupData$alphaExp # scaling for Ricker alpha wrt L_inf
  betaInt <- modelSetupData$betaInt # scaling for Ricker beta wrt Smax
  betaExp <- modelSetupData$betaExp # scaling for Ricker beta wrt Smax
  SmaxScale <-modelSetupData$SmaxScale # Catchability scaling for S.max estimate from survery
  recruitAlphas <- exp(alphaInt - abs(alphaExp*log(parameterValues$Linf)))
  recruitBetas <- exp(betaInt - betaExp*log(parameterValues$Smax*SmaxScale))
  recruitAlphas[1] <- modelSetupData$forageFishAlpha # trial for forage fish
  ##################################################################
  # M1 -parametes of beta function used to model M1
  alphaM1 <- modelSetupData$alphaM1
  betaM1 <- modelSetupData$betaM1
  cM1 <- modelSetupData$cM1 #scaling of final M1
  ##################################################################
  # Calculate upper and lower size class bins
  maxFishSize <- max(parameterValues$Linf) * 1.001
  lowScBin <- seq(from=0,to=maxFishSize-maxFishSize/nSizeClass,length.out = nSizeClass)
  uppScBin <- lowScBin + maxFishSize/nSizeClass
  midScBin <- lowScBin + (uppScBin-lowScBin)/2
  ##################################################################
  FW <- foodweb*modelSetupData$predationFlag
  initN <- t(initialValues)
  initN <- initN
  ##################################################################
  # size prefernce function Parameters
  spMu <- modelSetupData$spMu
  spSigma <- modelSetupData$spSigma

  #logical matrix reperesnteing size class bins applicable for each species
  scLinfMat <- sapply(parameterValues$Linf,function(x) {x>lowScBin})

  ##################################################################
  # Now calculate all parts required to run model
  ##################################################################
  ##################################################################
  # calculate the proportion leaving each size class per time step, and time step
  phi <- calc_phi(nSizeClass,nSpecies,uppScBin,lowScBin,parameterValues)

  ##################################################################
  # calculate the ration, weight, efficiency and largest size class
  ration <- calc_ration(nSizeClass,nSpecies,uppScBin,lowScBin,midScBin,phi$phiMin,parameterValues)
  ##################################################################
  # calculate maturity
  mature <- calc_maturity(nSizeClass,nSpecies,midScBin,scLinfMat,ration$scLinf,parameterValues)
  ##################################################################
  # calculate M1 (residual natural mortality)
  M1 <- calc_M1(nSizeClass,nSpecies,lowScBin,midScBin,alphaM1,betaM1,cM1,scLinfMat,ration$scLinf,phi$phiMin,parameterValues)
  ################################################################
  # calculated the size preference and the suitabilities
  M2PrefSuit <- calc_sizepref_suitability(nSizeClass,nSpecies,midScBin,spMu,spSigma,ration$wgt,ration$scLinf,FW)
  ##################################################################
  # calculates the predation mortalities
  #M2calcs <- calc_M2_r(nSizeClass,nSpecies,initN,ration,M2PrefSuit$suitability,phi$phiMin,otherFood)
  M2calcs <- calc_M2_c(nSizeClass,nSpecies,initN,ration$scLinf,ration$ration,ration$wgt,M2PrefSuit$suitability,phi$phiMin,modelSetupData$otherFood)
  ##################################################################
  # calculates the fishing mortalities
  # this will be extended to deal with more than one fleet
  eF <- calc_F(nSizeClass,nSpecies,midScBin,lowScBin,Ffull,modelSetupData$Falpha,modelSetupData$FL50,ration$scLinf,scLinfMat,phi$phiMin,parameterValues)
  ##################################################################
  # calculate Recruits and SSB
  recruits <- calc_recruits(initN,mature,ration$wgt,recruitAlphas,recruitBetas)

  ##################################################################
  # preallocate variables
  nTimeStepsPerYear <- round(1/phi$phiMin)
  nTimeSteps <- nYrs*nTimeStepsPerYear

  N <- array(data=NA,dim=c(nSizeClass,nSpecies,nTimeSteps))
  catch <- array(data=NA,dim=c(nSizeClass,nSpecies,nTimeSteps))
  M2 <- array(data=NA,dim=c(nSizeClass,nSpecies,nTimeSteps))
  R <- array(data=NA,dim=c(nSpecies,nYrs))
  SSB <- array(data=NA,dim=c(nSpecies,nYrs))

  ##################################################################
  # Now run the model
  ##################################################################
  # initialize time step 1
  iyear <- 1
  N[,,1] <- initN
  M2[,,1] <- M2calcs
  numberDead <- N[,,1]*(1-exp(-(eF+M1+M2[,,1])))
  numberRemain <- N[,,1]-numberDead
  # proportion dead due to fishing
  catch[,,1] <- (eF/(eF+M1+M2[,,1])) * numberDead
  #population growth
  updatedN <- calc_population_growth(nSizeClass,nSpecies,numberRemain,phi$probGrowOut)
  N[,,1] <- updatedN
  # move forward in time
  for (istep in 2:nTimeSteps){
    # calculate the growth of individuals
    N[,,istep] <- N[,,istep-1]

    # if we have entered the last step of the year (end of year) then a Recruitment event occurs
    if ((istep %% nTimeStepsPerYear) == 0) {
      recruits <- calc_recruits(N[,,istep],mature,ration$wgt,recruitAlphas,recruitBetas)
      R[,iyear] <- recruits$recruits
      SSB[,iyear] <- recruits$SSB
      # recruits added to the 1st size class
      N[1,,istep] <- N[1,,istep] + R[,iyear]
      iyear <- iyear + 1
    }

    # calculate predation mortality, M2
    M2calcs <- calc_M2_c(nSizeClass,nSpecies,N[,,istep],ration$scLinf,ration$ration,ration$wgt,M2PrefSuit$suitability,phi$phiMin,modelSetupData$otherFood)
    #M2calcs <- calc_M2_r(nSizeClass,nSpecies,N[,,istep],ration,M2PrefSuit$suitability,phi$phiMin,otherFood)
    M2[,,istep] <- M2calcs
    # calculate catch
    numberDead <- N[,,istep]*(1-exp(-(eF+M1+M2[,,istep])))
    # proportion number dead due to fishing
    catch[,,istep] <- (eF/(eF+M1+M2[,,istep])) * numberDead
    numberRemain <- N[,,istep]-numberDead
    updatedN <- calc_population_growth(nSizeClass,nSpecies,numberRemain,phi$probGrowOut)
    N[,,istep] <- updatedN
  }

  ##################################################################
  # end of run
  ##################################################################


  print(Sys.time()-start)
  return(list(N=N,M2=M2,catch=catch,SSB=SSB,Recruits=R))

}
