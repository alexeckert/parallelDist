\name{NEWS}
\title{NEWS for Package \pkg{parallelDist}}
\section{Changes in parallelDist version 0.2.2}{
  \itemize{
    \item Fixed DTW calculation errors for matrices of different length when using multiple threads.
    \item Reduced number of tests executed on CRAN
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
