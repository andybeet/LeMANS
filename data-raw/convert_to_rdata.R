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
  row.names(data_initialValues) <- data_species$commonName

  # read in parameter values
  data_parameterValues <- read.table(paste0(path,"GB_KeyRun_Params_SMAX2011_Linf.dat"))
  row.names(data_parameterValues) <- data_species$commonName
  names(data_parameterValues) <- c("k","Linf","Lmat","kappa","wa","wb","Smax","IsFished")


  save(data_species,file="data/data_species.RData")
  save(data_foodweb,file="data/data_foodweb.RData")
  save(data_initialValues,file="data/data_initialValues.RData")
  save(data_parameterValues,file="data/data_parameterValues.RData")



  return()

}
