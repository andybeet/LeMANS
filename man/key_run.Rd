% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/key_run.R
\name{key_run}
\alias{key_run}
\title{Key_Run (Length-based Multispecies Analysis by numerical Simulation)}
\usage{
key_run(Ffull, nYrs, modelSetup, parameterValues, initialValues, foodweb,
  species)
}
\arguments{
\item{Ffull}{Fishing mortaliy rate for a fully recruited fish}

\item{nYrs}{Number of years to simulate.}

\item{modelSetup}{list of parameters required for model (See \code{\link{rochet_GB_modelSetup}})}

\item{parameterValues}{matrix of species specific parameters (See \code{\link{rochet_GB_parameterValues}})}

\item{initialValues}{matrix of initial abundance estimates (See \code{\link{rochet_GB_initialValues}})}

\item{foodweb}{predator prey food web (See \code{\link{rochet_GB_foodweb}})}

\item{species}{matrix of species names and guild membership (See \code{\link{rochet_GB_species}})}
}
\value{
A list containing:

   \item{N}{3D array of abundance (numbers of animals). nsizeClass x nSpecies x nTimeSteps}

   \item{M1}{Matrix of M1 mortality ("natural"). nsizeClass x nSpecies}

   \item{M2}{3D array M2 mortality (predation). nsizeClass x nSpecies x nTimeSteps}

   \item{F}{3D array M2 mortality (predation). nsizeClass x nSpecies x nTimeSteps}

   \item{catch}{3D array of catch (numbers of animals). nsizeClass x nSpecies x nTimeSteps}

   \item{SSB}{Matrix of spawning stock biomass (SSB). nSpecies x nYears}

   \item{recruits    }{Matrix of recruits (numbers of animals). nSpecies x nYears}

   \item{suitability}{3D array of prey suitability by predator. nsizeClass x nSpecies x (num Pred.size class combinations). See \code{\link{calc_sizepref_suitability}}}

   \item{sizePreference}{3D array of prey preference by predator. nsizeClass x nSpecies x (num Pred.size class combinations). See \code{\link{calc_sizepref_suitability}}}

   \item{ration}{Matrix. Amount consumed to account for growth. See \code{\link{calc_phi}}}

   \item{growthEfficiency}{Matrix of growth efficiencies. See \code{\link{calc_phi}}}

   \item{growthProportions}{Matrix of proportions. Proportion of individuals that leave each size class in each time step. See \code{\link{calc_phi}}}

   \item{maturityProportions}{Matrix of maturity proportions. See \code{\link{calc_maturity}}}

   \item{modelTimeStep}{scalar representing the fraction of a year each time step represents. See \code{\link{calc_phi}}}
}
\description{
Runs the LeMANS model with the bundled data.
Bundled data conist of the 22 species used in the Rochet et al (2011) paper to represent Georges Bank. See \code{\link{rochet_GB_species}}.
The code was adapted directly from MATLAB code used in the Hall et al (2006) and Rochet et al (2011) papers.
}
\section{Using the Rochet et al data set}{

The number of size classes and the width of the size class was decided upon a priori.
This decision was based on the maximum L_inf among all species. max(Linf) = 148 cm
The \code{rochet_GB_initialValues} should be set up to represent the number of size classes and width.

The unit of output for Rochet et al data:

   catch, N:  number of individuals

   M1, M2, eF: rates

   recruits and SSB: were scaled to common recruit and spawning stock size units (individuals ?? 10^6 and tonnes ?? 10^3, respectively)
}

\examples{
\dontrun{
#'# runs the model with bundled data from Rochet et al (2011).
output <- key_run(Ffull=.4,nYrs=50,rochet_GB_modelSetup,rochet_GB_parameterValues,rochet_GB_initialValues,rochet_GB_foodweb,rochet_GB_species)
}
}
\seealso{
\code{\link{plot_key_run}},  \code{\link{rochet_GB_foodweb}},  \code{\link{rochet_GB_initialValues}},  \code{\link{rochet_GB_parameterValues}},  \code{\link{rochet_GB_species}}
}
