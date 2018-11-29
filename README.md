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

## Usage - STILL IN DEVELOPMENT

### Installation

devtools::install_github("andybeet/LeMANS",build_vignettes = TRUE)

### Help

browseVignettes("LeMANS")

The code was developed with a similar structure to the MATLAB code. Supporting data files are lazily loaded and ready to use.

* data_foodweb  - Predator-prey matrix (binary). Predator (column), prey (row)
* data_initialValues  - number of individuals per tow for each species/sizeclass combination
* data_parameterValues - carrying capacity (k), Maximum length (Linf) + for each species
* data_species - lists common name, scientific name, and guild each species is a member of.
* data_modelSetup - A list containing species independent parameter values (otherFood, Falpha, FL50)

To run the model with parameter values described in the Rochet et al paper simply type:

results <- key_run(Ffull=0.4, nYears=50, data_modelSetup, data_parameterValues, data_initialValues, data_foodweb, data_species)

To view a list of the functions:


