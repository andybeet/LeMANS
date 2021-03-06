% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calc_M2_c}
\alias{calc_M2_c}
\title{Calculates predation mortality, M2 (C++)}
\usage{
calc_M2_c(nSize, nSpecies, N, scLinf, ration, wgt, suitability, phiMin,
  otherFood)
}
\arguments{
\item{nSize}{Number of size class intervals species can grow through}

\item{nSpecies}{Number of species in the model}

\item{N}{A matrix (nSize x nSpecies) of abundance (number of individuals)}

\item{scLinf}{A vector (nSpecies) of the size class at which each species reaches L_inf (maximum length)}

\item{ration}{A matrix (nSize x nSpecies) of (growth in time interval)/growth efficiency values. See \code{\link{calc_ration}}}

\item{wgt}{A matrix (nSize x nSpecies) of species weight at the mid point of each size class (Units: grams). See \code{\link{calc_ration}}}

\item{suitability}{3D array of predator size preference for prey size. See \code{\link{calc_sizepref_suitability}}}

\item{phiMin}{Scalar. Model timestep (years). See \code{\link{calc_phi}}}

\item{otherFood}{Scalar. Amount of other food available not accounted for in the model (g)}
}
\value{
A matrix is returned

  \item{M2}{A matrix (nSize x nSpecies) of M2 (predation mortality) values.
   Note: M2 = 0 for size classes in which a species is not preyed upon.}
}
\description{
Calculates the predation mortality for each species in each size class
}
\section{References}{

Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63:1344-1359.

Rochet et al. (2011). Does selective fishing conserve community biodiversity? Prediction from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68:469-486
alc_
}

\seealso{
\code{\link{calc_sizepref_suitability}}, \code{\link{calc_ration}} \code{\link{calc_phi}}
}
