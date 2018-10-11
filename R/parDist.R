## parDist.R
##
## Copyright (C)  2017  Alexander Eckert
##
## This file is part of parallelDist.
##
## parallelDist is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## parallelDist is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with parallelDist. If not, see <http://www.gnu.org/licenses/>.

#
# Calculates distance matrices in parallel
#
parDist <- parallelDist <- function (x, method = "euclidean", diag = FALSE, upper = FALSE, threads = NULL, ...) {
  METHODS <- c("bhjattacharyya", "bray", "canberra", "chord", "divergence",
               "dtw", "euclidean", "fJaccard", "geodesic", "hellinger",
               "kullback", "mahalanobis", "manhattan", "maximum", "minkowski",
               "podani", "soergel", "wave", "whittaker",
               "binary", "braun-blanquet", "dice", "fager", "faith",
               "hamman", "kulczynski1", "kulczynski2", "michael", "mountford",
               "mozley", "ochiai", "phi", "russel", "simple matching",
               "simpson", "stiles", "tanimoto", "yule", "yule2", "cosine", "custom") # w/o "levenshtein"
  methodIdx <- pmatch(method, METHODS)
  if (is.na(methodIdx))
    stop("Invalid distance method")
  method <- METHODS[methodIdx]

  arguments <- list(...)
  # set step pattern (for dtw distances)
  step.pattern.name <- getStepPatternName(arguments)
  if (!any(is.na(step.pattern.name))) {
    arguments[["step.pattern"]] <- step.pattern.name[1]
  }

  # check funct argument for custom distance measure
  if (method == "custom") {
    funcPtr = arguments[["func"]]
    if (is.null(funcPtr)) {
      stop("Parameter 'func' is missing.")
    }
    checkPtr(funcPtr)
  }

  # set number of threads
  if (!is.null(threads)) {
    RcppParallel::setThreadOptions(numThreads = threads)
  }

  N <- ifelse(is.list(x), length(x), nrow(x))
  attrs <- list(Size = N, Labels = names(x), Diag = diag, Upper = upper,
                method = METHODS[methodIdx], call = match.call(), class = "dist")

  # check data type
  if (is.list(x) && inherits(x, "list")) {
    methods.first.row.only <- c("chord", "geodesic", "podani")
    if (method %in% methods.first.row.only) {
      warning("Only first row of each matrix is used for distance calculation.")
    }
    return(.Call('_parallelDist_cpp_parallelDistVec', PACKAGE = 'parallelDist', x, attrs, arguments = arguments))
  } else {
    if (is.matrix(x)) {
      return(.Call('_parallelDist_cpp_parallelDistMatrixVec', PACKAGE = 'parallelDist', x, attrs, arguments = arguments))
    } else {
      stop("x must be a matrix or a list of matrices.")
    }
  }
}

getType <- function (code) {
  tokenize <- strsplit(code, "[[:space:]]*(\\(|\\)){1}[[:space:]]*")[[1]]
  tokens <- strsplit(tokenize[[1]], "[[:space:]]+")[[1]]
  tokens <- tokens[seq_len(length(tokens) - 1)]
  paste(tokens[tokens != ""], collapse = " ")
}

getStepPatternName <- function(arguments) {
  step.pattern.name <- NA
  supported.patterns.names <- c("asymmetric", "asymmetricP0", "asymmetricP05", "asymmetricP1", "asymmetricP2",
                                "symmetric1", "symmetric2", "symmetricP0", "symmetricP05", "symmetricP1", "symmetricP2")
  sp.candidate <- arguments[["step.pattern"]]

  if (!is.null(sp.candidate)) {
    if (class(sp.candidate) == "stepPattern" && requireNamespace("dtw", quietly = TRUE))  {
      supported.patterns <- list(dtw::asymmetric, dtw::asymmetricP0, dtw::asymmetricP05, dtw::asymmetricP1, dtw::asymmetricP2,
                                 dtw::symmetric1, dtw::symmetric2, dtw::symmetricP0, dtw::symmetricP05, dtw::symmetricP1, dtw::symmetricP2)
      # check if step pattern is supported (using object)
      sp.found <- sapply(supported.patterns, FUN = function(x){
        if (identical(dim(x), dim(sp.candidate))) {
          all(x == sp.candidate)
        } else {
          FALSE
        }
      })
    } else {
      if (is.character(sp.candidate)) {
        # check if step pattern is supported (using name)
        sp.found <- supported.patterns.names == sp.candidate
      }
    }
    if (any(sp.found)) {
      step.pattern.name <- supported.patterns.names[which(sp.found)]
    } else {
      stop("Step pattern is not supported.")
    }
  }
  step.pattern.name
}

checkPtr <- function(ptr) {
  stopifnot(inherits(ptr, "XPtr"))

  actualReturnType <- attr(ptr, "type")
  expectedReturnType <- "double"

  actualArgTypes <- sapply(attr(ptr, "args"), getType, USE.NAMES = FALSE)
  expectedArgTypes <- rep("const arma::mat&", 2)

  msg <- character()
  # check return type
  if (actualReturnType != expectedReturnType) {
    msg <- paste(c(msg, paste0("  Wrong return type '", actualReturnType,  "', should be '", expectedReturnType, "'.")), collapse = "\n")
  }
  # check number of arguments
  if (length(actualArgTypes) != 2) {
    msg <- paste(c(msg, paste0("  Wrong number of arguments ('", length(actualArgTypes) ,"'), should be 2.")), collapse = "\n")
  } else {
    # check argument types
    for (i in which(!(expectedArgTypes == actualArgTypes))) {
      msg <- paste(c(msg, paste0("  Wrong argument type '", actualArgTypes[[i]], "', should be '", expectedArgTypes[[i]], "'.")), collapse = "\n")
    }
  }

  if (length(msg))
    stop("Bad XPtr signature:\n", msg)
}
