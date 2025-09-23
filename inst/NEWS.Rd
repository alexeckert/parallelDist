\name{NEWS}
\title{NEWS for Package \pkg{parallelDist}}
\section{Changes in parallelDist version 0.2.7}{
  \itemize{
    \item Removed C++11 dependency to support newer Armadillo versions
  }
}
\section{Changes in parallelDist version 0.2.6}{
  \itemize{
    \item Remove LazyData from DESCRIPTION
    \item Switch vignette font
  }
}
\section{Changes in parallelDist version 0.2.5}{
  \itemize{
    \item Reflect changes in proxy 0.4-26: Change coercion function for cosine similarity to distance to 1-x instead of 1-abs(x)
  }
}
\section{Changes in parallelDist version 0.2.4}{
  \itemize{
    \item Fixed build on Solaris using explicit casts.
  }
}
\section{Changes in parallelDist version 0.2.3}{
  \itemize{
    \item Preserve label attribute of dist class.
    \item Allow double values for Minkowski distance.
    \item Added hamming distance (contributed by Srikanth KS).
    \item Added cosine similarity (contributed by Lucas Godoy).
  }
}
\section{Changes in parallelDist version 0.2.2}{
  \itemize{
    \item Fixed DTW calculation errors for matrices of different length when using multiple threads.
    \item Reduced number of tests executed on CRAN.
  }
}
\section{Changes in parallelDist version 0.2.1}{
  \itemize{
    \item Fixed vignette generation errors for upcoming R releases.
    \item Removed highlight as vignette builder and dependency.
  }
}
\section{Changes in parallelDist version 0.2.0}{
  \itemize{
    \item Added support for user-defined distance functions.
	  \item Added R 3.0.2 as requirement.
	  \item Fixed clang UBSAN downcast warning.
  }
}
\section{Changes in parallelDist version 0.1.1}{
  \itemize{
    \item Missing include statement of the string library in Utility.h
    caused compile error when clang++ was used.
  }
}
