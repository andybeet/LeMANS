# LeMANS (Length-based Multispecies Analysis by Numerical Simulation)

![deploy to github pages](https://github.com/andybeet/LeMANS/workflows/deploy%20to%20github%20pages/badge.svg)
![gitleaks](https://github.com/andybeet/LeMANS/workflows/gitleaks/badge.svg)

Details of the model can be found in the following publications:

* Hall et al. (2006). A length-based multispecies model for evaluating community responses to fishing. Can. J. Fish. Aquat. Sci. 63:1344-1359. doi: [https://doi.org/10.1139/f06-039](https://doi.org/10.1139/f06-039)

* Rochet et al. (2011). Does selective fishing conserve community biodiversity? Prediction from a length-based multispecies model. Can. J. Fish. Aquat. Sci. 68:469-486. doi: [https://doi.org/10.1139/F10-159](https://doi.org/10.1139/F10-159)

Originally coded in MATLAB

Commments/questions/bugs should be documented at [https://github.com/NOAA-EDAB/LeMANS/issues](https://github.com/NOAA-EDAB/LeMANS/issues)

## Features

* Each component of the model (M1, M2, F, suitability, size preference etc) has been modularized such that a users code block can be substitued.

* original data files (used in Rochet et al. (2011)) are bundled with the package.

* The code is substantially faster (~9 times) than the original MATLAB implementation

* It serves as a framework for extending the model

## Usage

### Installation

``` r
# build the package from github with vignettes
remotes::install_github("NOAA-EDAB/LeMANS",build_vignettes = TRUE)
```
If you don't have remotes installed you will see an error "there is no package called 'remotes'"; if that happens install devtools with `install.packages("remotes")`.

If installation of LeMANS fails, it could be due to the Rcpp package. Try to install Rcpp first.

### Help

```r
# view the documentation on usage
browseVignettes("LeMANS")
```

The code was developed with a similar structure to the MATLAB code. Supporting data files are lazily loaded and ready to use.

* `rochet_GB_foodweb`  - Predator-prey matrix (binary). Predator (column), prey (row)
* `rochet_GB_initialValues`  - number of individuals per tow for each species/sizeclass combination
* `rochet_GB_parameterValues` - carrying capacity (k), Maximum length (Linf) + for each species
* `rochet_GB_species` - lists common name, scientific name, and guild each species is a member of.
* `rochet_GB_modelSetup` - A list containing species independent parameter values (otherFood, Falpha, FL50)

To run the model with parameter values described in the Rochet et al paper simply type:

`results <- key_run(Ffull=0.4, nYrs=50, rochet_GB_modelSetup, rochet_GB_parameterValues, rochet_GB_initialValues, rochet_GB_foodweb, rochet_GB_species)`

To plot the catch output aggregated over sizeClass:

`plot_key_run(results$catch/1E6, ylabel = "catch (millions individuals)", is.aggregated=T, speciesNames=rochet_GB_species, scale="free")`

## Contact

| [andybeet](https://github.com/andybeet)        
| ----------------------------------------------------------------------------------------------- 
| [![](https://avatars1.githubusercontent.com/u/22455149?s=100&v=4)](https://github.com/andybeet) | 



#### Legal disclaimer

*This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.*
