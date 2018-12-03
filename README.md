# LeMANS (Length-based Multispecies Analysis by Numerical Simulation)

Details of the model can be found in the follwing publications:

* Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63:1344-1359.

* Rochet et al. (2011). Does selective fishing conserve community biodiversity? Prediction from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68:469-486

Originally coded in MATLAB

## Features

* Each component of the model (M1, M2, F, suitability, size preference etc) has been modularized such that a users code block can be substitued.

* original data files (used in Rochet et al. (2011)) are bundled with the package.

* The code is substantially faster (~9 times) than the original MATLAB implementation

* It serves as a framework for extending the model

## Usage

### Installation

devtools::install_github("andybeet/LeMANS",build_vignettes = TRUE)

### Help

browseVignettes("LeMANS")

The code was developed with a similar structure to the MATLAB code. Supporting data files are lazily loaded and ready to use.

* rochet_GB_foodweb  - Predator-prey matrix (binary). Predator (column), prey (row)
* rochet_GB_initialValues  - number of individuals per tow for each species/sizeclass combination
* rochet_GB_parameterValues - carrying capacity (k), Maximum length (Linf) + for each species
* rochet_GB_species - lists common name, scientific name, and guild each species is a member of.
* rochet_GB_modelSetup - A list containing species independent parameter values (otherFood, Falpha, FL50)

To run the model with parameter values described in the Rochet et al paper simply type:

results <- key_run(Ffull=0.4, nYears=50, rochet_GB_modelSetup, rochet_GB_parameterValues, rochet_GB_initialValues, rochet_GB_foodweb, rochet_GB_species)

To plot the catch output aggregated over sizeClass:

plot_key_run(results$catch, ylabel = "catch (millions individuals)",is.aggregated=T,rochet_GB_species)


