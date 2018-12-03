## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,   comment = "#>" ,fig.width = 8, fig.height = 8,fig.align="center")
library(LeMANS)

## ----run1, echo=TRUE-----------------------------------------------------
    results <- key_run(Ffull=0.4, nYrs=50, rochet_GB_modelSetup, rochet_GB_parameterValues, rochet_GB_initialValues, rochet_GB_foodweb, rochet_GB_species)


## ----run2, echo=TRUE-----------------------------------------------------
    str(results)

## ----plot_c, echo = TRUE-------------------------------------------------
    plot_key_run(results$catch/1E6, ylabel = "catch (millions individuals)", is.aggregated=F, rochet_GB_species, scales="free")


## ----plot_n, echo = TRUE-------------------------------------------------
    plot_key_run(results$N/1E6, ylabel = "Abundance (millions individuals)", is.aggregated=T, rochet_GB_species, scales="free")


## ----plot_m2, echo = TRUE------------------------------------------------
    plot_key_run(results$M2, ylabel = "Predation mortality", is.aggregated=F, rochet_GB_species, scales="fixed")


