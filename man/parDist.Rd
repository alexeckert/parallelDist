\name{parDist}
\alias{parDist}
\alias{parallelDist}
\title{Parallel Distance Matrix Computation using multiple Threads}
\usage{
parDist(x, method = "euclidean", diag = FALSE, upper = FALSE, threads = NULL, ...)
parallelDist(x, method = "euclidean", diag = FALSE, upper = FALSE, threads = NULL, ...)
}
\arguments{
\item{x}{a numeric matrix (each row is one series) or list of numeric matrices for multidimensional series (each matrix is one series, a row is a dimension of a series)}

\item{method}{the distance measure to be used. A list of all available distance methods can be found in the details section below.}

\item{diag}{logical value indicating whether the diagonal of the distance matrix should be printed by print.dist.}

\item{upper}{logical value indicating whether the upper triangle of the distance matrix should be printed by print.dist}

\item{threads}{number of cpu threads for calculating a distance matrix. Default is the maximum amount of cpu threads available on the system.}

\item{...}{additional parameters which will be passed to the distance methods. See details section below.}

}
\description{
Calculates distance matrices in parallel using multiple threads. Supports 39 predefined distance measures and user-defined distance functions.
}

\details{
\subsection{User-defined distance functions}{
  \describe{
      \item{\code{custom}}{
      Defining and compiling a user-defined C++ distance function, as well as creating an external pointer to the function can easily be achieved with the \code{\link[RcppXPtrUtils]{cppXPtr}} function of the \pkg{RcppXPtrUtils} package. The resulting \verb{Xptr} external pointer object needs to be passed to \code{\link[parallelDist]{parDist}} using the \code{func} parameter.

Parameters:
          \itemize{
            \item{
              \describe{
                \item{\code{func} (Xptr)}{External pointer to a user-defined distance function with the following signature: \cr \code{double customDist(const arma::mat &A, const arma::mat &B)} \cr Note that the return value must be a \verb{double} and the two parameters must be of type \verb{const arma::mat &param}.  \cr \cr
                More information about the Armadillo library can be found at \url{http://arma.sourceforge.net/docs.html} or as part of the documentation of the \pkg{RcppArmadillo} package.}
              }
            }
          }
An exemplary definition and usage of an user-defined euclidean distance function can be found in the examples section below.
    }
  }
}

\subsection{Available predefined distance measures (written for two vectors \eqn{x} and \eqn{y})}{

  \bold{Distance methods for continuous input variables}

  \describe{
    \item{\code{bhjattacharyya}}{
      The Bhjattacharyya distance.\cr Type: continuous\cr Formula: \eqn{sqrt(sum_i (sqrt(x_i) - sqrt(y_i))^2))}.\cr Details: See \command{pr_DB$get_entry("bhjattacharyya")} in \pkg{proxy}.
    }
    \item{\code{bray}}{
      The Bray/Curtis dissimilarity.\cr Type: continuous\cr Formula: \eqn{sum_i |x_i - y_i| / sum_i (x_i + y_i)}.\cr Details: See \command{pr_DB$get_entry("bray")} in \pkg{proxy}.
    }
    \item{\code{canberra}}{
      The Canberra distance (with compensation for excluded components). Terms with zero numerator and denominator are omitted from the sum and treated as if the values were missing. \cr Type: continuous\cr Formula: \eqn{sum_i |x_i - y_i| / |x_i + y_i|}.\cr Details: See \command{pr_DB$get_entry("canberra")} in \pkg{proxy}.
    }
    \item{\code{chord}}{
      The Chord distance.\cr Type: continuous\cr Formula: \eqn{sqrt(2 * (1 - xy / sqrt(xx * yy)))}.\cr Details: See \command{pr_DB$get_entry("chord")} in \pkg{proxy}.
    }
    \item{\code{divergence}}{
      The Divergence distance.\cr Type: continuous\cr Formula: \eqn{sum_i (x_i - y_i)^2 / (x_i + y_i)^2}.\cr Details: See \command{pr_DB$get_entry("divergence")} in \pkg{proxy}.
    }

    \item{\code{dtw}}{Implementation of a multi-dimensional Dynamic Time Warping algorithm.\cr Type: continuous\cr Formula: Euclidean distance \eqn{sqrt(sum_i (x_i - y_i)^2)}.\cr
    Parameters:
          \itemize{
            \item{
              \describe{
                \item{\code{window.size} (integer, optional)}{Size of the window of the Sakoe-Chiba band. If the absolute length difference of two series x and y is larger than the window.size, the window.size is set to the length difference.}
              }
            }
            \item{
              \describe{
                \item{\code{norm.method} (character, optional)}{Normalization method for DTW distances.
                  \itemize{
                    \item{\code{path.length} Normalization with the length of the warping path.}
                    \item{\code{n} Normalization with n. n is the length of series x.}
                    \item{\code{n+m} Normalization with n + m. n is the length of series x, m is the length of series y.}
                  }
                }
              }
            }
            \item{
              \describe{
                \item{\code{step.pattern} (character or stepPattern of \pkg{dtw} package, default: \code{symmetric1})}{
                The following step patterns of the \pkg{dtw} package are supported:

                   \itemize{
                      \item{
                        \code{asymmetric} (Normalization hint: n)
                      }
                      \item{
                        \code{asymmetricP0} (Normalization hint: n)
                      }
                      \item{
                        \code{asymmetricP05} (Normalization hint: n)
                      }
                      \item{
                        \code{asymmetricP1} (Normalization hint: n)
                      }
                      \item{
                        \code{asymmetricP2} (Normalization hint: n)
                      }
                      \item{
                        \code{symmetric1} (Normalization hint: path.length)
                      }
                      \item{
                        \code{symmetric2} or \code{symmetricP0} (Normalization hint: n+m)
                      }
                      \item{
                        \code{symmetricP05} (Normalization hint: n+m)
                      }
                      \item{
                        \code{symmetricP1} (Normalization hint: n+m)
                      }
                      \item{
                        \code{symmetricP2} (Normalization hint: n+m)
                      }
                  }
              For a detailed description see \code{\link[dtw]{stepPattern}} of the \pkg{dtw} package.
              }
            }
        }
      }
    }

    \item{\code{euclidean}}{
      The Euclidean distance/L_2-norm (with compensation for excluded components).\cr Type: continuous\cr Formula: \eqn{sqrt(sum_i (x_i - y_i)^2))}.\cr Details: See \command{pr_DB$get_entry("euclidean")} in \pkg{proxy}.
    }
    \item{\code{fJaccard}}{
      The fuzzy Jaccard distance.\cr Type: binary\cr Formula: \eqn{sum_i (min{x_i, y_i}) / sum_i(max{x_i, y_i})}.\cr Details: See \command{pr_DB$get_entry("fJaccard")} in \pkg{proxy}.
    }
    \item{\code{geodesic}}{
      The geoedesic distance, i.e. the angle between x and y.\cr Type: continuous\cr Formula: \eqn{arccos(xy / sqrt(xx * yy))}.\cr Details: See \command{pr_DB$get_entry("geodesic")} in \pkg{proxy}.
    }
    \item{\code{hellinger}}{
      The Hellinger distance.\cr Type: continuous\cr Formula: \eqn{sqrt(sum_i (sqrt(x_i / sum_i x) - sqrt(y_i / sum_i y)) ^ 2)}.\cr Details: See \command{pr_DB$get_entry("hellinger")} in \pkg{proxy}.
    }
    \item{\code{kullback}}{
      The Kullback-Leibler distance.\cr Type: continuous\cr Formula: \eqn{sum_i [x_i * log((x_i / sum_j x_j) / (y_i / sum_j y_j)) / sum_j x_j)]}.\cr Details: See \command{pr_DB$get_entry("kullback")} in \pkg{proxy}.
    }
    \item{\code{mahalanobis}}{
      The Mahalanobis distance. The Variance-Covariance-Matrix is estimated from the input data if unspecified.\cr Type: continuous\cr Formula: \eqn{sqrt((x - y) Sigma^(-1) (x - y))}.\cr Parameters:
          \itemize{
            \item{
              \describe{
                \item{\code{cov} (numeric matrix, optional)}{The covariance matrix  (p x p) of the distribution.}
              }
            }
            \item{
              \describe{
                \item{\code{inverted} (logical, optional)}{If TRUE, cov is supposed to contain the inverse of the covariance matrix.}
              }
            }
          }
          Details: See \command{pr_DB$get_entry("mahalanobis")} in \pkg{proxy} or \command{mahalanobis} in \pkg{stats}.
    }
    \item{\code{manhattan}}{
      The Manhattan/City-Block/Taxi/L_1-norm distance (with compensation for excluded components).\cr Type: continuous\cr Formula: \eqn{sum_i |x_i - y_i|}.\cr Details: See \command{pr_DB$get_entry("manhattan")} in \pkg{proxy}.
    }
    \item{\code{maximum}}{
      The Maximum/Supremum/Chebyshev distance.\cr Type: continuous\cr Formula: \eqn{max_i |x_i - y_i|}.\cr Details: See \command{pr_DB$get_entry("maximum")} in \pkg{proxy}.
    }
    \item{\code{minkowski}}{
      The Minkowski distance/p-norm (with compensation for excluded components). \cr Type: continuous\cr Formula: \eqn{(sum_i (x_i - y_i)^p)^(1/p)}.\cr Parameters:
          \itemize{
            \item{
              \describe{
                \item{\code{p} (integer, optional)}{The \eqn{p}th root of the sum of the \eqn{p}th powers of the differences of the components.}
              }
            }
          }
          Details: See \command{pr_DB$get_entry("minkowski")} in \pkg{proxy}.
    }

    \item{\code{podani}}{
      The Podany measure of discordance is defined on ranks with ties. In the formula, for two given objects x and y, n is the number of variables, a is is the number of pairs of variables ordered identically, b the number of pairs reversely ordered, c the number of pairs tied in both x and y (corresponding to either joint presence or absence), and d the number of all pairs of variables tied at least for one of the objects compared such that one, two, or thee scores are zero.\cr Type: continuous\cr Formula: \eqn{1 - 2 * (a - b + c - d) / (n * (n - 1))}.\cr Details: See \command{pr_DB$get_entry("podani")} in \pkg{proxy}.
    }
    \item{\code{soergel}}{
      The Soergel distance.\cr Type: continuous\cr Formula: \eqn{sum_i |x_i - y_i| / sum_i max{x_i, y_i}}.\cr Details: See \command{pr_DB$get_entry("soergel")} in \pkg{proxy}.
    }
    \item{\code{wave}}{
      The Wave/Hedges distance.\cr Type: continuous\cr Formula: \eqn{sum_i (1 - min(x_i, y_i) / max(x_i, y_i))}.\cr Details: See \command{pr_DB$get_entry("wave")} in \pkg{proxy}.
    }
    \item{\code{whittaker}}{
      The Whittaker distance.\cr Type: continuous\cr Formula: \eqn{sum_i |x_i / sum_i x - y_i / sum_i y| / 2}.\cr Details: See \command{pr_DB$get_entry("whittaker")} in \pkg{proxy}.
    }
}

\bold{Distance methods for binary input variables}

\emph{Notation:}
\itemize{
  \item a: number of (TRUE, TRUE) pairs
  \item b: number of (FALSE, TRUE) pairs
  \item c: number of (TRUE, FALSE) pairs
  \item d: number of (FALSE, FALSE) pairs
}

\emph{Note:} Similarities are converted to distances.

\describe{
    \item{\code{binary}}{
      The Jaccard Similarity for binary data. It is the proportion of (TRUE, TRUE) pairs, but not considering (FALSE, FALSE) pairs.\cr Type: binary\cr Formula: \eqn{a / (a + b + c)}.\cr Details: See \command{pr_DB$get_entry("binary")} in \pkg{proxy}.
    }
    \item{\code{braun-blanquet}}{
      The Braun-Blanquet similarity.\cr Type: binary\cr Formula: \eqn{a / max{(a + b), (a + c)}}.\cr Details: See \command{pr_DB$get_entry("braun-blanquet")} in \pkg{proxy}.
    }
    \item{\code{dice}}{
      The Dice similarity.\cr Type: binary\cr Formula: \eqn{2a / (2a + b + c)}.\cr Details: See \command{pr_DB$get_entry("dice")} in \pkg{proxy}.
    }
    \item{\code{fager}}{
      The Fager / McGowan distance.\cr Type: binary\cr Formula: \eqn{a / sqrt((a + b)(a + c)) - sqrt(a + c) / 2}.\cr Details: See \command{pr_DB$get_entry("fager")} in \pkg{proxy}.
    }
    \item{\code{faith}}{
      The Faith similarity.\cr Type: binary\cr Formula: \eqn{(a + d/2) / n}.\cr Details: See \command{pr_DB$get_entry("faith")} in \pkg{proxy}.
    }
    \item{\code{hamman}}{
      The Hamman Matching similarity for binary data. It is the proportion difference of the concordant and discordant pairs.\cr Type: binary\cr Formula: \eqn{([a + d] - [b + c]) / n}.\cr Details: See \command{pr_DB$get_entry("hamman")} in \pkg{proxy}.
    }
    \item{\code{kulczynski1}}{
      Kulczynski similarity for binary data. Relates the (TRUE, TRUE) pairs to discordant pairs.\cr Type: binary\cr Formula: \eqn{a / (b + c)}.\cr Details: See \command{pr_DB$get_entry("kulczynski1")} in \pkg{proxy}.
    }
    \item{\code{kulczynski2}}{
      Kulczynski similarity for binary data. Relates the (TRUE, TRUE) pairs to the discordant pairs.\cr Type: binary\cr Formula: \eqn{[a / (a + b) + a / (a + c)] / 2}.\cr Details: See \command{pr_DB$get_entry("kulczynski2")} in \pkg{proxy}.
    }
    \item{\code{michael}}{
      The Michael similarity.\cr Type: binary\cr Formula: \eqn{4(ad - bc) / [(a + d)^2 + (b + c)^2]}.\cr Details: See \command{pr_DB$get_entry("michael")} in \pkg{proxy}.
    }
    \item{\code{mountford}}{
      The Mountford similarity for binary data.\cr Type: binary\cr Formula: \eqn{2a / (ab + ac + 2bc)}.\cr Details: See \command{pr_DB$get_entry("mountford")} in \pkg{proxy}.
    }
    \item{\code{mozley}}{
      The Mozley/Margalef similarity.\cr Type: binary\cr Formula: \eqn{an / (a + b)(a + c)}.\cr Details: See \command{pr_DB$get_entry("mozley")} in \pkg{proxy}.
    }
    \item{\code{ochiai}}{
      The Ochiai similarity.\cr Type: binary\cr Formula: \eqn{a / sqrt[(a + b)(a + c)]}.\cr Details: See \command{pr_DB$get_entry("ochiai")} in \pkg{proxy}.
    }
    \item{\code{phi}}{
      The Phi similarity (= Product-Moment-Correlation for binary variables).\cr Type: binary\cr Formula: \eqn{(ad - bc) / sqrt[(a + b)(c + d)(a + c)(b + d)]}.\cr Details: See \command{pr_DB$get_entry("phi")} in \pkg{proxy}.
    }
    \item{\code{russel}}{
      The Russel/Raosimilarity for binary data. It is just the proportion of (TRUE, TRUE) pairs.\cr Type: binary\cr Formula: \eqn{a / n}.\cr Details: See \command{pr_DB$get_entry("russel")} in \pkg{proxy}.
    }
    \item{\code{simple matching}}{
      The Simple Matching similarity for binary data. It is the proportion of concordant pairs.\cr Type: binary\cr Formula: \eqn{(a + d) / n}.\cr Details: See \command{pr_DB$get_entry("simple matching")} in \pkg{proxy}.
    }
    \item{\code{simpson}}{
      The Simpson similarity.\cr Type: binary\cr Formula: \eqn{a / min{(a + b), (a + c)}}.\cr Details: See \command{pr_DB$get_entry("simpson")} in \pkg{proxy}.
    }
    \item{\code{stiles}}{
      The Stiles similarity. Identical to the logarithm of Krylov's distance.\cr Type: binary\cr Formula: \eqn{log(n(|ad-bc| - 0.5n)^2 / [(a + b)(c + d)(a + c)(b + d)])}.\cr Details: See \command{pr_DB$get_entry("stiles")} in \pkg{proxy}.
    }
    \item{\code{tanimoto}}{
      The Rogers/Tanimoto similarity for binary data. Similar to the simple matching coefficient, but putting double weight on the discordant pairs.\cr Type: binary\cr Formula: \eqn{(a + d) / (a + 2b + 2c + d)}.\cr Details: See \command{pr_DB$get_entry("tanimoto")} in \pkg{proxy}.
    }
    \item{\code{yule}}{
     The Yule similarity.\cr Type: binary\cr Formula: \eqn{(ad - bc) / (ad + bc)}.\cr Details: See \command{pr_DB$get_entry("yule")} in \pkg{proxy}.
    }
    \item{\code{yule2}}{
     The Yule similarity.\cr Type: binary\cr Formula: \eqn{(sqrt(ad) - sqrt(bc)) / (sqrt(ad) + sqrt(bc))}.\cr Details: See \command{pr_DB$get_entry("yule2")} in \pkg{proxy}.
    }
  }
  }
}


\value{
  \code{parDist} returns an object of class \code{"dist"}.

  The lower triangle of the distance matrix stored by columns in a
  vector, say \code{do}. If \code{n} is the number of
  observations, i.e., \code{n <- attr(do, "Size")}, then
  for \eqn{i < j \le n}, the dissimilarity between (row) i and j is
  \code{do[n*(i-1) - i*(i-1)/2 + j-i]}.
  The length of the vector is \eqn{n*(n-1)/2}, i.e., of order \eqn{n^2}.

  The object has the following attributes (besides \code{"class"} equal
  to \code{"dist"}):
  \item{Size}{integer, the number of observations in the dataset.}
  \item{Labels}{optionally, contains the labels, if any, of the
    observations of the dataset.}
  \item{Diag, Upper}{logicals corresponding to the arguments \code{diag}
    and \code{upper} above, specifying how the object should be printed.}
  \item{call}{optionally, the \code{\link{call}} used to create the
    object.}
  \item{method}{optionally, the distance method used; resulting from
    \code{\link{parDist}()}, the (\code{\link{match.arg}()}ed) \code{method}
    argument.}
}

\examples{
\dontrun{
## predefined distance functions
# defining a matrix, where each row corresponds to one series
sample.matrix <- matrix(c(1:100), ncol = 10)

# euclidean distance
parDist(x = sample.matrix, method = "euclidean")
# minkowski distance with parameter p=2
parDist(x = sample.matrix, method = "minkowski", p=2)
# dynamic time warping distance
parDist(x = sample.matrix, method = "dtw")
# dynamic time warping distance normalized with warping path length
parDist(x = sample.matrix, method = "dtw", norm.method="path.length")
# dynamic time warping with different step pattern
parDist(x = sample.matrix, method = "dtw", step.pattern="symmetric2")
# dynamic time warping with window size constraint
parDist(x = sample.matrix, method = "dtw", step.pattern="symmetric2", window.size=1)

## multi-dimensional distance functions using list of matrices
# defining a list of matrices, where each list entry row corresponds to a two dimensional series
tmp.mat <- matrix(c(1:40), ncol = 10)
sample.matrix.list <- list(tmp.mat[1:2,], tmp.mat[3:4,])

# multi-dimensional euclidean distance
parDist(x = sample.matrix.list, method = "euclidean")
# multi-dimensional dynamic time warping
parDist(x = sample.matrix.list, method = "dtw")

## user-defined distance function
library(RcppArmadillo)
# Use RcppXPtrUtils for simple usage of C++ external pointers
library(RcppXPtrUtils)

# compile user-defined function and return pointer (RcppArmadillo is used as dependency)
euclideanFuncPtr <- cppXPtr(
"double customDist(const arma::mat &A, const arma::mat &B) {
  return sqrt(arma::accu(arma::square(A - B)));
}", depends = c("RcppArmadillo"))

# distance matrix for user-defined euclidean distance function (note that method is set to "custom")
parDist(matrix(1:16, ncol=2), method="custom", func = euclideanFuncPtr)}
}
