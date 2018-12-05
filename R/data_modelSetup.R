#' List of species independent parameters
#'
#' A list of parameter values that are not species specific that are required for the model to run.
#' All values are Scalar values.
#'
#' @format A list of species independent parameters:
#' \describe{
#'
#'      \item{otherFood}{Amount of other food available not accounted for in the model (g)}
#'      \item{predationFlag}{A flag denoting if predation is on 1 or off 0 in model}
#'      \item{Falpha}{Steepness of selectivity curve}
#'      \item{FL50}{Length at 50\% selectivity}
#'      \item{alphaInt}{Scaling for Ricker alpha wrt L_inf}
#'      \item{alphaExp}{Scaling for Ricker alpha wrt L_inf}
#'      \item{betaInt}{Scaling for Ricker beta wrt Smax}
#'      \item{betaExp}{Scaling for Ricker beta wrt Smax}
#'      \item{SmaxScale}{Catchability scaling for Smax estimate from survery}
#'      \item{forageFishAlpha}{Trial for forage fish}
#'      \item{alphaM1}{beta pdf parameter to model M1}
#'      \item{betaM1}{beta pdf parameter to model M1}
#'      \item{cM1}{Scaling of pdf to model M1}
#'      \item{spMu}{Mean of size preference function }
#'      \item{spSigma}{Standard deviation of size preference function}
#'
#' }
#'
#'@seealso \code{\link{rochet_GB_foodweb}}, \code{\link{rochet_GB_initialValues}}, \code{\link{rochet_GB_parameterValues}}, \code{\link{rochet_GB_species}}, \code{\link{rochet_GB_modelSetup}}
#'
#'@source Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63: 1344-1359
#'@source Rochet et al. (2011). Does selective fishing conserve community biodiversiy? Predictions from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68: 469-486
"rochet_GB_modelSetup"
