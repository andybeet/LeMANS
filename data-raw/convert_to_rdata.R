# converts Jeremies files to rdata files

convert_to_rdata <- function() {
  path <- paste0("data-raw/")

  #reads in species list
  data_species <- read.csv(paste0(path,"speciesID.csv"),header=TRUE)
  # reads in Foodweb structure
  data_foodweb <- read.table(paste0(path,"GB_KeyRun_FW.dat"))
  names(data_foodweb) <- data_species$commonName
  row.names(data_foodweb) <- data_species$commonName

  # read in initial values
  data_initialValues <- read.table(paste0(path,"GB_KeyRun_InitVal.dat"))
  data_initialValues <- data_initialValues * 3046527 # number of tows: converts catch per tow to total numbers
  row.names(data_initialValues) <- data_species$commonName

  # read in parameter values
  data_parameterValues <- read.table(paste0(path,"GB_KeyRun_Params_SMAX2011_Linf.dat"))
  row.names(data_parameterValues) <- data_species$commonName
  names(data_parameterValues) <- c("k","Linf","Lmat","kappa","wa","wb","Smax","IsFished")


  save(data_species,file="data/data_species.RData")
  save(data_foodweb,file="data/data_foodweb.RData")
  save(data_initialValues,file="data/data_initialValues.RData")
  save(data_parameterValues,file="data/data_parameterValues.RData")

  # set all default values in a list for the user to change
  data_modelSetup <- list()
  data_modelSetup$otherFood <- 55000000 #(grams)
  data_modelSetup$predationFlag <- 1 # predation on or off in model
  data_modelSetup$Falpha <- .25 # steepness of selectivity curve
  data_modelSetup$FL50 <- 25 # length at 50% selection
  data_modelSetup$alphaInt <- 11 # scaling for Ricker alpha wrt L_inf
  data_modelSetup$alphaExp <- -2.298 # scaling for Ricker alpha wrt L_inf
  data_modelSetup$betaInt <- .1513 # scaling for Ricker beta wrt Smax
  data_modelSetup$betaExp <- .9484 # scaling for Ricker beta wrt Smax
  data_modelSetup$SmaxScale <- 1 # Catchability scaling for Smax estimate from survery
  data_modelSetup$forageFishAlpha <- 400 # trial for forage fish
  data_modelSetup$alphaM1 <- .8 # beta func parameter to model M1
  data_modelSetup$betaM1 <- .4 # beta func parameter to model M1
  data_modelSetup$cM1 <- .35 # scaling of beta func to model M1
  data_modelSetup$spMu <- 0.5   # size prefernce function mean
  data_modelSetup$spSigma <- 2  # size prefernce function standard deviation


  save(data_modelSetup,file="data/data_modelSetup.RData")




  return()

}
