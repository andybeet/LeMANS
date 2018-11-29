#' LeMANS: Length-based Multispecies Analysis by Numerical Simulation
#'
#'An implementation of the LeMANS model as described in the Hall et al (2006) and Rochet et al (2011) publications.
#'
#'\itemize{
#'\item Each component of the model (M1, M2, F, suitability, size preference etc) has been modularized
#'    such that a users code block can be substitued.
#'\item original data files (used in Rochet et al. (2011)) are bundled with the package.
#'\item The code is substantially faster (~9 times) than the original MATLAB implementation
#'\item It serves as a framework for extending the model
#'}
#'
#'To learn more about using \code{LeMANS}, start with the vignette: \code{browseVignettes(package="LeMANS")}
#'
#'
#' @docType package
#' @name LeMANS
#' @useDynLib LeMANS
#' @importFrom Rcpp sourceCpp
NULL
