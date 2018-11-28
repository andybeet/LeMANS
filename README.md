# LeMANS (Length-based Multispecies Analysis by Numerical Simulation)

Details of the model can be found in the follwing publications:

* Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63:1344-1359.

* Rochet et al. (2011). Does selective fishing conserve community biodiversity? Prediction from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68:469-486

Originally coded in MATLAB

## Usage

### Installation

devtools::install_github("andybeet/LeMANS",build_vignettes = TRUE)

### Help

browseVignettes("LeMANS")

The code was developed with a similar structure to the MATLAB code. Supporting data files are lazily loaded and ready to use.

* data_foodweb  - Predator-prey matrix (binary). Predator (column), prey (row)
* data_initialValues  - number of individuals per tow for each species/sizeclass combination
* data_parameterValues - carrying capacity (k), Maximum length (Linf) + for each species
* data_species - lists common name, scientific name, and guild each species is a member of.

To run the model with parameter values described in the Rochet et al paper simply type:

