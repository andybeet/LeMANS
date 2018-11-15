# converts Jeremies files to rdata files

convert_to_rdata <- function() {
  path <- paste0("data-raw/")

  #reads in species list
  species <- read.csv(paste0(path,"speciesID.csv"),header=TRUE)
  # reads in Foodweb structure
  foodweb <- read.table(paste0(path,"GB_KeyRun_FW.dat"))
  names(foodweb) <- species$speciesName
  row.names(foodweb) <- species$speciesName

  # read in initial values
  initialValues <- read.table(paste0(path,"GB_KeyRun_InitVal.dat"))
  row.names(initialValues) <- species$speciesName

  # read in parameter values
  parameterValues <- read.table(paste0(path,"GB_KeyRun_Params_SMAX2011_Linf.dat"))
  row.names(parameterValues) <- species$speciesName
  names(parameterValues) <- c("k","Linf","Lmat","kappa","wa","wb","Smax","IsFished")


  save(species,file="data/species.RData")
  save(foodweb,file="data/foodweb.RData")
  save(initialValues,file="data/initialValues.RData")
  save(parameterValues,file="data/parameterValues.RData")



  return()

}
