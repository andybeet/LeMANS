#'Key_Run (Length-based Multispecies Analysis by numerical Simulation)
#'
#'Runs the LeMANS model with the bundled data.
#'Bundled data conist of the 22 species used in the Rochet et al (2011) paper to represent Georges Bank. See \code{\link{rochet_GB_species}}.
#'The code was adapted directly from MATLAB code used in the Hall et al (2006) and Rochet et al (2011) papers.
#'
#'@param Ffull Fishing mortaliy rate for a fully recruited fish
#'@param nYrs Number of years to simulate.
#'@param modelSetup list of parameters required for model (See \code{\link{rochet_GB_modelSetup}})
#'@param parameterValues matrix of species specific parameters (See \code{\link{rochet_GB_parameterValues}})
#'@param initialValues matrix of initial abundance estimates (See \code{\link{rochet_GB_initialValues}})
#'@param foodweb predator prey food web (See \code{\link{rochet_GB_foodweb}})
#'@param species matrix of species names and guild membership (See \code{\link{rochet_GB_species}})
#'
#'@return A list containing:
#'
#'    \item{N}{3D array of abundance (numbers of animals). nsizeClass x nSpecies x nTimeSteps}
#'
#'    \item{M1}{Matrix of M1 mortality ("natural"). nsizeClass x nSpecies}
#'
#'    \item{M2}{3D array M2 mortality (predation). nsizeClass x nSpecies x nTimeSteps}
#'
#'    \item{F}{3D array M2 mortality (predation). nsizeClass x nSpecies x nTimeSteps}
#'
#'    \item{catch}{3D array of catch (numbers of animals). nsizeClass x nSpecies x nTimeSteps}
#'
#'    \item{SSB}{Matrix of spawning stock biomass (SSB). nSpecies x nYears}
#'
#'    \item{recruits    }{Matrix of recruits (numbers of animals). nSpecies x nYears}
#'
#'    \item{suitability}{3D array of prey suitability by predator. nsizeClass x nSpecies x (num Pred.size class combinations). See \code{\link{calc_sizepref_suitability}}}
#'
#'    \item{sizePreference}{3D array of prey preference by predator. nsizeClass x nSpecies x (num Pred.size class combinations). See \code{\link{calc_sizepref_suitability}}}
#'
#'    \item{ration}{Matrix. Amount consumed to account for growth. See \code{\link{calc_phi}}}
#'
#'    \item{growthEfficiency}{Matrix of growth efficiencies. See \code{\link{calc_phi}}}
#'
#'    \item{growthProportions}{Matrix of proportions. Proportion of individuals that leave each size class in each time step. See \code{\link{calc_phi}}}
#'
#'    \item{maturityProportions}{Matrix of maturity proportions. See \code{\link{calc_maturity}}}
#'
#'    \item{modelTimeStep}{scalar representing the fraction of a year each time step represents. See \code{\link{calc_phi}}}
#'
#'
#'
#'
#'@section Using the Rochet et al data set:
  #'The number of size classes and the width of the size class was decided upon a priori.
  #'This decision was based on the maximum L_inf among all species. max(Linf) = 148 cm
  #'The \code{rochet_GB_initialValues} should be set up to represent the number of size classes and width.
  #'
  #'The unit of output for Rochet et al data:
  #'
  #'    catch, N:  number of individuals
  #'
  #'    M1, M2, eF: rates
  #'
  #'    recruits and SSB: were scaled to common recruit and spawning stock size units (individuals × 10^6 and tonnes × 10^3, respectively)
  #'
  #'
  #'
#'@seealso \code{\link{plot_key_run}},  \code{\link{rochet_GB_foodweb}},  \code{\link{rochet_GB_initialValues}},  \code{\link{rochet_GB_parameterValues}},  \code{\link{rochet_GB_species}}
#'
#'
#'@examples
#'\dontrun{
#'#'# runs the model with bundled data from Rochet et al (2011).
#' output <- key_run(Ffull=.4,nYrs=50,rochet_GB_modelSetup,rochet_GB_parameterValues,rochet_GB_initialValues,rochet_GB_foodweb,rochet_GB_species)
#'}
#'@export



# Note: all data fines are read in to memory when package is loaded
# initialValues, parameterValues, species, foodweb
key_run <- function(Ffull,nYrs,modelSetup,parameterValues,initialValues,foodweb,species) {
  # initial set up
  nSizeClass <- dim(initialValues)[2]
  nSpecies <- dim(initialValues)[1]
  ##################################################################
  # calculate recruitment. Assumes all has a ricker form. All are scaled.
  # see Hall et al paper
  alphaInt <- modelSetup$alphaInt # scaling for Ricker alpha wrt L_inf
  alphaExp <- modelSetup$alphaExp # scaling for Ricker alpha wrt L_inf
  betaInt <- modelSetup$betaInt # scaling for Ricker beta wrt Smax
  betaExp <- modelSetup$betaExp # scaling for Ricker beta wrt Smax
  SmaxScale <-modelSetup$SmaxScale # Catchability scaling for S.max estimate from survery
  recruitAlphas <- exp(alphaInt - abs(alphaExp*log(parameterValues$Linf)))
  recruitBetas <- exp(betaInt - betaExp*log(parameterValues$Smax*SmaxScale))
  recruitAlphas[1] <- modelSetup$forageFishAlpha # trial for forage fish
  ##################################################################
  # Calculate upper and lower size class bins
  maxFishSize <- max(parameterValues$Linf) * 1.001
  lowScBin <- seq(from=0,to=maxFishSize-maxFishSize/nSizeClass,length.out = nSizeClass)
  uppScBin <- lowScBin + maxFishSize/nSizeClass
  midScBin <- lowScBin + (uppScBin-lowScBin)/2
  ##################################################################
  FW <- foodweb*modelSetup$predationFlag
  initN <- t(initialValues)
  initN <- initN

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
  M1 <- calc_M1(nSizeClass,nSpecies,lowScBin,midScBin,modelSetup$alphaM1,modelSetup$betaM1,modelSetup$cM1,scLinfMat,ration$scLinf,phi$phiMin,parameterValues)
  ################################################################
  # calculated the size preference and the suitabilities
  M2PrefSuit <- calc_sizepref_suitability(nSizeClass,nSpecies,midScBin,modelSetup$spMu,modelSetup$spSigma,ration$wgt,ration$scLinf,FW)
  ##################################################################
  # calculates the predation mortalities
  #M2calcs <- calc_M2_r(nSizeClass,nSpecies,initN,ration,M2PrefSuit$suitability,phi$phiMin,otherFood)
  M2calcs <- calc_M2_c(nSizeClass,nSpecies,initN,ration$scLinf,ration$ration,ration$wgt,M2PrefSuit$suitability,phi$phiMin,modelSetup$otherFood)
  ##################################################################
  # calculates the fishing mortalities
  # this will be extended to deal with more than one fleet
  eF <- calc_F(nSizeClass,nSpecies,midScBin,lowScBin,Ffull,modelSetup$Falpha,modelSetup$FL50,ration$scLinf,scLinfMat,phi$phiMin,parameterValues)
  ##################################################################
  # calculate Recruits and SSB
  #recruits <- calc_recruits(initN,mature,ration$wgt,recruitAlphas,recruitBetas)

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
      # calculate Recruits and SSB
      recruits <- calc_recruits(N[,,istep],mature,ration$wgt,recruitAlphas,recruitBetas)
      R[,iyear] <- recruits$recruits
      SSB[,iyear] <- recruits$SSB
      # recruits added to the 1st size class
      N[1,,istep] <- N[1,,istep] + R[,iyear]
      iyear <- iyear + 1
    }

    # calculate predation mortality, M2
    M2calcs <- calc_M2_c(nSizeClass,nSpecies,N[,,istep],ration$scLinf,ration$ration,ration$wgt,M2PrefSuit$suitability,phi$phiMin,modelSetup$otherFood)
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


  return(list(N=N,M1=M1,M2=M2,eF=eF,catch=catch,SSB=SSB,recruits=R,
              modelTimeStep=phi$phiMin,growthProportions=phi$probGrowOut,
              ration=ration$ration,growthEfficiency=ration$gEff,maturityProportions=mature,
              suitability=M2PrefSuit$suitability,sizePreference=M2PrefSuit$sizePref))

}
