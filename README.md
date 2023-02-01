# bremory (Version: 1.0.0) [beta]
### Modeling Memory Retention in Relational Event Data
### Table of contents
- [About the package](#about-the-package)
- [Programming Languages](#programming-languages)
- [Installing the package](#installing-the-package)
- [Vignettes](#vignettes)
- [Funding](#funding)
- [Author](#author)

### About the package
The `bremory` package offers several methods for inquiring about the presence of memory decay in Relational Event Networks:
 * bayesian semi-parametric approach for modeling memory decay in relational event networks  ([full-text here](https://doi.org/10.1177/00491241221113875), _this is the only method available in the current version, which is still in a beta version_)
 * parametric approach for modeling memory decay in one-type relational event networks (under review, _the method will be available in the development branch_)
 * parametric approach for modeling type-related memory decays in relational event networks (in progress, _the method will be available in the development branch_)

### Programming Languages
The package contains code written in:
* R (>= 4.0.0)
* Rcpp (>= 1.0.8.3) and RcppArmadillo (>= 0.11)
* C++14
	
### Installing the package
To install the package in R using `devtools`:

```
library(devtools)
devtools::install_github(repo = "TilburgNetworkGroup/bremory", build_vignettes = TRUE)

# load the package
library(bremory)
```

### Vignettes
In order to provide a thorough explanation over the functions inside the package, the user can read the vignette with examples and explanation

```
vignette(package = "bremory") # returns all the vignettes available with the current version of the package 
```

### Funding
The funder of this work is the ERC and the ERC project number is 758791.

### Author
Giuseppe Arena, Tilburg University (Tilburg, The Netherlands). (g.arena@tilburguniversity.edu)
