# converts Jeremies files to rdata files

convert_to_rdata <- function() {
  path <- paste0("data-raw/")

  #reads in species list
  rochet_GB_species <- read.csv(paste0(path,"speciesID.csv"),header=TRUE)
  # reads in Foodweb structure
  rochet_GB_foodweb <- read.table(paste0(path,"GB_KeyRun_FW.dat"))
  names(rochet_GB_foodweb) <- rochet_GB_species$commonName
  row.names(rochet_GB_foodweb) <- rochet_GB_species$commonName

  # read in initial values
  rochet_GB_initialValues <- read.table(paste0(path,"GB_KeyRun_InitVal.dat"))
  rochet_GB_initialValues <- rochet_GB_initialValues * 3046527 # number of tows: converts catch per tow to total numbers
  row.names(rochet_GB_initialValues) <- rochet_GB_species$commonName

  # read in parameter values
  rochet_GB_parameterValues <- read.table(paste0(path,"GB_KeyRun_Params_SMAX2011_Linf.dat"))
  row.names(rochet_GB_parameterValues) <- rochet_GB_species$commonName
  names(rochet_GB_parameterValues) <- c("k","Linf","Lmat","kappa","wa","wb","Smax","IsFished")


  save(rochet_GB_species,file="data/rochet_GB_species.RData")
  save(rochet_GB_foodweb,file="data/rochet_GB_foodweb.RData")
  save(rochet_GB_initialValues,file="data/rochet_GB_initialValues.RData")
  save(rochet_GB_parameterValues,file="data/rochet_GB_parameterValues.RData")

  # set all default values in a list for the user to change
  rochet_GB_modelSetup <- list()
  rochet_GB_modelSetup$otherFood <- 55000000 #(grams)
  rochet_GB_modelSetup$predationFlag <- 1 # predation on or off in model
  rochet_GB_modelSetup$Falpha <- .25 # steepness of selectivity curve
  rochet_GB_modelSetup$FL50 <- 25 # length at 50% selection
  rochet_GB_modelSetup$alphaInt <- 11 # scaling for Ricker alpha wrt L_inf
  rochet_GB_modelSetup$alphaExp <- -2.298 # scaling for Ricker alpha wrt L_inf
  rochet_GB_modelSetup$betaInt <- .1513 # scaling for Ricker beta wrt Smax
  rochet_GB_modelSetup$betaExp <- .9484 # scaling for Ricker beta wrt Smax
  rochet_GB_modelSetup$SmaxScale <- 1 # Catchability scaling for Smax estimate from survery
  rochet_GB_modelSetup$forageFishAlpha <- 400 # trial for forage fish
  rochet_GB_modelSetup$alphaM1 <- .8 # beta func parameter to model M1
  rochet_GB_modelSetup$betaM1 <- .4 # beta func parameter to model M1
  rochet_GB_modelSetup$cM1 <- .35 # scaling of beta func to model M1
  rochet_GB_modelSetup$spMu <- 0.5   # size prefernce function mean
  rochet_GB_modelSetup$spSigma <- 2  # size prefernce function standard deviation


  save(rochet_GB_modelSetup,file="data/rochet_GB_modelSetup.RData")




  return()

}
