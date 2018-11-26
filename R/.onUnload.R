# R packages p 101. H.Wickham

.onUnload <- function (libpath) {
  library.dynam.unload("LeMANS", libpath)
}
