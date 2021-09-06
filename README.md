# An updated guide for the perplexed: cytometry in the high-dimensional era 

This repository contains scripts and data to generate Figures 2D-F. The scripts here utilise the [Spectre](https://immunedynamics.io/spectre) toolkit in R. See [this page](https://immunedynamics.io/spectre/getting-started/) for installation instructions. 

Briefly, Spectre can be installed via `devtools`. In R/RStudio:

```
# Install devtools
if(!require('devtools')) {install.packages('devtools')}
```

```
# Install Spectre
library('devtools')
options(timeout=6000)
devtools::install_github("immunedynamics/spectre")
```
