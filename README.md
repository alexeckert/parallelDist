## parallelDist [![CRAN](http://www.r-pkg.org/badges/version/parallelDist)](https://CRAN.R-project.org/package=parallelDist) [![Build Status](https://travis-ci.org/alexeckert/parallelDist.svg?branch=master)](https://travis-ci.org/alexeckert/parallelDist) [![Build status](https://ci.appveyor.com/api/projects/status/d6o2d529gdf7qjyu/branch/master?svg=true)](https://ci.appveyor.com/project/alexeckert/paralleldist/branch/master) [![codecov](https://codecov.io/gh/alexeckert/parallelDist/branch/master/graph/badge.svg)](https://codecov.io/gh/alexeckert/parallelDist) [![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2Falexeckert%2FparallelDist.svg?type=shield)](https://app.fossa.io/projects/git%2Bgithub.com%2Falexeckert%2FparallelDist?ref=badge_shield)


### Introduction

The [parallelDist](https://CRAN.R-project.org/package=parallelDist) package provides a fast parallelized alternative to R's native 'dist' function to calculate distance matrices for continuous, binary, and multi-dimensional input matrices and offers a broad variety of predefined distance functions from the 'stats', 'proxy' and 'dtw' R packages, as well as support for user-defined distance functions written in C++. For ease of use, the 'parDist' function extends the signature of the 'dist' function and uses the same parameter naming conventions as distance methods of existing R packages. Currently 41 different distance methods are supported.

The package is mainly implemented in C++ and leverages the '[Rcpp](https://CRAN.R-project.org/package=Rcpp)' and '[RcppParallel](https://CRAN.R-project.org/package=RcppParallel)' package to parallelize the distance computations with the help of the 'TinyThread' library. Furthermore, the Armadillo linear algebra library is used via '[RcppArmadillo](https://CRAN.R-project.org/package=RcppArmadillo)' for optimized matrix operations for distance calculations. The curiously recurring template pattern (CRTP) technique is applied to avoid virtual functions, which improves the Dynamic Time Warping calculations while keeping the implementation flexible enough to support different step patterns and normalization methods.

### Documentation and Usage Examples

Usage examples and performance benchmarks can be found in the included vignette.

Details about the 41 supported distance methods and their parameters are described on the help page of the 'parDist' function. The help page can be displayed with the following command:

```R
?parDist
```

#### User-defined distance functions

Since version 0.2.0, parallelDist supports fast parallel distance matrix computations for user-defined distance functions written in C++.

A user-defined function needs to have the following signature (also see the [Armadillo documentation](http://arma.sourceforge.net/docs.html)):

```Cpp
double customDist(const arma::mat &A, const arma::mat &B)
```

Defining and compiling the function, as well as creating an external pointer to the user-defined function can easily be achieved with the *cppXPtr* function of the '[RcppXPtrUtils](https://CRAN.R-project.org/package=RcppXPtrUtils)' package. The following code shows a full example of defining and using a user-defined euclidean distance function:

```R
# RcppArmadillo is used as dependency
library(RcppArmadillo)
# RcppXPtrUtils is used for simple handling of C++ external pointers
library(RcppXPtrUtils)

# compile user-defined function and return pointer (RcppArmadillo is used as dependency)
euclideanFuncPtr <- cppXPtr("double customDist(const arma::mat &A, const arma::mat &B) { return sqrt(arma::accu(arma::square(A - B))); }",
                            depends = c("RcppArmadillo"))

# distance matrix for user-defined euclidean distance function
# (note that method is set to "custom")
parDist(matrix(1:16, ncol=2), method="custom", func = euclideanFuncPtr)
```

More information can be found in the vignette and the help pages.

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


[![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2Falexeckert%2FparallelDist.svg?type=large)](https://app.fossa.io/projects/git%2Bgithub.com%2Falexeckert%2FparallelDist?ref=badge_large)