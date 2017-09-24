## testMatrixCustomDistances.R
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

context("Custom distance methods using matrix as input")

mat.sample1 <- matrix(c(0,1,0,1,0,0,1,0), nrow = 2)
mat.sample2 <- matrix(c(0,1,0,1,0,0,1,0,1,1), nrow = 2)
mat.sample3 <- matrix(c(1:500), ncol = 5)
mat.sample4 <- matrix(rep(0,100), ncol = 5)
mat.sample5 <- matrix(c(-500:499), ncol = 5)
mat.sample6 <- matrix(c(1:2), ncol = 1)
mat.sample7 <- matrix(c(0.5,1,0,1,0,0,1,0.3,1,1), nrow = 2)

mat.list <- list(mat.sample1, mat.sample2, mat.sample3, mat.sample4, mat.sample5, mat.sample6, mat.sample7)

testMatrixEquality <- function(matrix, comparisonMethod, func, ...) {
  expect_equal(as.matrix(parDist(matrix, method = "custom", func=func, ...)), as.matrix(parDist(matrix, method = comparisonMethod, ...)))
}

testMatrixListEquality <- function(matlist, comparisonMethod, func, ...) {
  invisible(sapply(matlist, function(x) { testMatrixEquality(x, comparisonMethod, func, ...) }))
}

library(RcppXPtrUtils)
test_that("custom euclidean distance method produces same outputs as native method", {
  ptr <- cppXPtr("double customDist(const arma::mat &A, const arma::mat &B) { return sqrt(arma::accu(arma::square(A - B))); }", depends = c("RcppArmadillo"))
  testMatrixListEquality(mat.list, "euclidean", func=ptr)
})

test_that("error for missing func parameter shows up", {
  expect_error(parDist(mat.sample1, method = "custom"), "Parameter 'func' is missing.")
})

test_that("error for wrong return value of func pointer shows up", {
  ptr <- cppXPtr("int customDist(const arma::mat &A, const arma::mat &B) { return 0; }", depends = c("RcppArmadillo"))
  expect_error(parDist(mat.sample1, method = "custom", func=ptr), "Wrong return type 'int', should be 'double'.")
})

test_that("error for wrong argument type of func pointer shows up", {
  ptr <- cppXPtr("double customDist(int a, const arma::mat &B) { return 0; }", depends = c("RcppArmadillo"))
  expect_error(parDist(mat.sample1, method = "custom", func=ptr), "Wrong argument type 'int', should be 'const arma::mat&'.")
})

test_that("error for wrong number of arguments of func pointer shows up", {
#  ptr <- ""
#  class(ptr) <- "XPtr"
#  attr(ptr, "type")="double"
#  attr(ptr, "fname")="customDist"
#  attr(ptr, "args")="const arma::mat& A"
  ptr <- cppXPtr("double customDist(const arma::mat &A) { return 0; }", depends = c("RcppArmadillo"))
  expect_error(parDist(mat.sample1, method = "custom", func=ptr), "Wrong number of arguments \\('1'\\), should be 2.")
})
