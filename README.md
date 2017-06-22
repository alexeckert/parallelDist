## parallelDist [![Build Status](https://travis-ci.org/alexeckert/parallelDist.svg?branch=master)](https://travis-ci.org/alexeckert/parallelDist) [![CRAN](http://www.r-pkg.org/badges/version/parallelDist)](https://CRAN.R-project.org/package=parallelDist)

### Introduction

The [parallelDist](https://CRAN.R-project.org/package=parallelDist) package provides a fast parallelized alternative to R's native 'dist' function to calculate distance matrices for continuous, binary, and multi-dimensional input matrices and offers a broad variety of distance functions from the 'stats', 'proxy' and 'dtw' R packages. For ease of use, the 'parDist' function extends the signature of the 'dist' function and uses the same parameter naming conventions as distance methods of existing R packages. Currently 39 different distance methods are supported.

The package is mainly implemented in C++ and leverages the '[Rcpp](https://CRAN.R-project.org/package=Rcpp)' and '[RcppParallel](https://CRAN.R-project.org/package=RcppParallel)' package to parallelize the distance computations with the help of the 'TinyThread' library. Furthermore, the Armadillo linear algebra library is used via '[RcppArmadillo](https://CRAN.R-project.org/package=RcppArmadillo)' for optimized matrix operations for distance calculations. The curiously recurring template pattern (CRTP) technique is applied to avoid virtual functions, which improves the Dynamic Time Warping calculations while keeping the implementation flexible enough to support different step patterns and normalization methods.

### Documentation and Usage Examples

Usage examples and performance benchmarks can be found in the included vignette.

Details about the 39 supported distance methods and their parameters are described on the help page of the 'parDist' function. The help page can be displayed with the following command:

```R
?parDist
```

### Installation

parallelDist is available on [CRAN](https://CRAN.R-project.org/package=parallelDist) and can be installed with the following command: 

```R
install.packages("parallelDist")
```

The current version from github can be installed using the 'devtools' package:

```R
library(devtools)
install_github("alexeckert/parallelDist")
```

### Authors

Alexander Eckert

### License

GPL (>= 2)
