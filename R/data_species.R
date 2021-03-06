#' Species descriptions
#'
#' A table of species information
#'
#' @format A data frame of size nSpecies x 4. Species in rows. The columns are defined as:
#' \describe{
#'     \item{speciesNumber}{Ordered sequence of species numbers }
#'     \item{commonName}{Common name for species}
#'     \item{scientificName}{Latin name for species}
#'     \item{guild}{Associated guild name for species}
#'}
#'
#'@section Species in model:
#'Forage Fish,
#'Spiny dogfish,
#'Winter skate,
#'Little skate,
#'Atlantic herring,
#'Silver hake,
#'Atlantic cod,
#'Haddock,
#'Pollock,
#'White hake,
#'Red hake,
#'Fourspot flounder,
#'Summer flounder,
#'Yellowtail founder,
#'Winter flounder,
#'Witch flounder,
#'Windowpane flounder,
#'Atlantic mackerel,
#'Longhorn sculpin,
#'Sea raven,
#'Goosefish,
#'Sandlance
#'
#'@seealso \code{\link{rochet_GB_foodweb}}, \code{\link{rochet_GB_initialValues}}, \code{\link{rochet_GB_parameterValues}}, \code{\link{rochet_GB_species}}, \code{\link{rochet_GB_modelSetup}}
#'
#'@source Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63: 1344-1359
#'@source Rochet et al. (2011). Does selective fishing conserve community biodiversiy? Predictions from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68: 469-486
"rochet_GB_species"
