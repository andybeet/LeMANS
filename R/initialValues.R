#' Ititial abundance values
#'
#' starting conditions for all species in each of the species size classes
#'
#' @format A data frame of size 22 x 15
#' \describe{abundance (# of individuals) in each species size class. Size classes are of equal lengt, the width determined by Linf (cm).
#' Species in rows, size classes in columns.
#' The list of species can be found in \code{\link{species}} and in the references below
#'}
#'
#'@seealso \code{\link{foodweb}}, \code{\link{initialValues}}, \code{\link{parameterValues}}, \code{\link{species}}
#'
#'@source Halle et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63: 1344-1359
#'@source Rochet et al. (2011). Does selective fishing conserve community biodiversiy? Predictions from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68: 469-486
"initialValues"
