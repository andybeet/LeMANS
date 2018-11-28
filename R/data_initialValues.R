#' Ititial abundance values
#'
#' starting conditions for all species in each of the species size classes
#'
#' @format A data frame of size 22 x 15
#' \describe{abundance (#catch per tow) in each species size class. Size classes are of equal length; the width (cm) determined by Linf.
#' Species in rows, size classes in columns.
#' The list of species can be found in \code{\link{data_species}} and in the references below
#'}
#'
#'@seealso \code{\link{data_foodweb}}, \code{\link{data_initialValues}}, \code{\link{data_parameterValues}}, \code{\link{data_species}}, \code{\link{data_modelSetup}}
#'
#'@source Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63: 1344-1359
#'@source Rochet et al. (2011). Does selective fishing conserve community biodiversiy? Predictions from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68: 469-486
"data_initialValues"
