#' Parameter values used in the model
#'
#' Species specific parameter values
#'
#' @format A data frame of size nSpecies x 8. The columns are defined as:
#' \describe{
#'     \item{k}{Von Bertalanfy growth parameter (year^-1)}
#'     \item{Linf}{Maximum length species can grow (cm)}
#'     \item{Lmat}{Length at maturity (cm)}
#'     \item{kappa}{Curvature parameter for the maturity ogive}
#'     \item{w_a}{Parameter of lenth-weight conversion}
#'     \item{w_b}{Parameter of lenth-weight conversion}
#'     \item{Smax}{Maximum observed spawning stock biomass (kg)}
#'     \item{IsFished}{Flag indicating if a species is fished or not}
#'}
#' The list of species can be found in \code{\link{rochet_GB_species}} and in the references below.
#'
#'@seealso \code{\link{rochet_GB_foodweb}}, \code{\link{rochet_GB_initialValues}}, \code{\link{rochet_GB_parameterValues}}, \code{\link{rochet_GB_species}}, \code{\link{rochet_GB_modelSetup}}
#'
#'@source Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63: 1344-1359
#'@source Rochet et al. (2011). Does selective fishing conserve community biodiversiy? Predictions from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68: 469-486
"rochet_GB_parameterValues"
