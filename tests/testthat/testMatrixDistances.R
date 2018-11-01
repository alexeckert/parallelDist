## testMatrixDistances.R
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

context("Distance methods using matrix as input")

mat.sample1 <- matrix(c(0, 1, 0, 1, 0, 0, 1, 0), nrow = 2)
mat.sample2 <- matrix(c(0, 1, 0, 1, 0, 0, 1, 0, 1, 1), nrow = 2)
mat.sample3 <- matrix(c(1:500), ncol = 5)
mat.sample4 <- matrix(rep(0, 100), ncol = 5)
mat.sample5 <- matrix(c(-500:499), ncol = 5)
mat.sample6 <- matrix(c(1:2), ncol = 1)
mat.sample7 <- matrix(c(0.5, 1, 0, 1, 0, 0, 1, 0.3, 1, 1), nrow = 2)

mat.list <- list(mat.sample1, mat.sample2, mat.sample3, mat.sample4, mat.sample5, mat.sample6, mat.sample7)

if (isCran()) {
  mat.list <- mat.list[1:4]
}

testMatrixEquality <- function(matrix, method, ...) {
  expect_equal(as.matrix(parDist(matrix, method = method, ...)), as.matrix(dist(matrix, method = method, ...)))
}

library(proxy)
testMatrixListEquality <- function(matlist, method, ...) {
  invisible(sapply(matlist, function(x) {
    testMatrixEquality(x, method, ...)
  }))
}

test_that("error for invalid distance method shows up", {
  expect_error(parDist(mat.sample1, method = "unknown"), "Invalid distance method")
})

test_that("error for invalid input type", {
  expect_error(parDist(data.frame(c(1:2))), "x must be a matrix or a list of matrices.")
})

# works
test_that("bhjattacharyya method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "bhjattacharyya")
})
# works
test_that("bray method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "bray")
})
# works
test_that("canberra method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "canberra")
})
# works
test_that("chord method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "chord")
})
# works
test_that("divergence method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "divergence")
})
# works
test_that("euclidean method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "euclidean")
})
# works (only vals [0,1])
test_that("fJaccard method produces same outputs as dist", {
  testMatrixListEquality(list(mat.sample1, mat.sample2, mat.sample7), "fJaccard")
})
# works
test_that("geodesic method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "geodesic")
})
# works
test_that("hellinger method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "hellinger")
})
# works
test_that("kullback method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "kullback")
})
#works
test_that("mahalanobis method produces same outputs as dist", {
  mat.mahalanobis <- cbind(1:6, 1:3)
  testMatrixEquality(mat.mahalanobis, "mahalanobis")
  testMatrixEquality(mat.mahalanobis, "mahalanobis", cov = cov(mat.mahalanobis))
  expect_equal(as.matrix(parDist(mat.mahalanobis, "mahalanobis", cov = solve(cov(mat.mahalanobis)), inverted = TRUE)),
               as.matrix(dist(mat.mahalanobis, method = "mahalanobis", cov = cov(mat.mahalanobis))))
})
test_that("mahalanobis method throws error for list of matrices input", {
  mat.mahalanobis <- cbind(1:6, 1:3)
  mat.list <- list(mat.mahalanobis, mat.mahalanobis)
  expect_error(parDist(mat.list, "mahalanobis"),
               "Calculation of inverted covariance matrix is only supported for input data in matrix format.")
})
# works
test_that("manhattan method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "manhattan")
})
# works
test_that("maximum method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "maximum")
})
# works
test_that("minkowski method produces same outputs as dist", {
  for (p_param in c(1:10)) {
    testMatrixListEquality(mat.list, "minkowski", p = as.numeric(p_param))
  }
})
# works
test_that("podani method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "podani")
})
# possible wrong implementation in proxy?
test_that("soergel method produces expected output", {
  expect_equal(as.matrix(parDist(mat.sample1, method = "soergel"))[1, 2], 1)
  expect_equal(as.matrix(parDist(mat.sample2, method = "soergel"))[1, 2], 0.75)
  expect_equal(as.matrix(parDist(mat.sample3, method = "soergel"))[1, 2], 0.0049504950495049506)
  expect_equal(as.matrix(parDist(mat.sample4, method = "soergel"))[1, 2], NaN)
  expect_equal(as.matrix(parDist(mat.sample5, method = "soergel"))[1, 2], -0.010101010101010102)
  expect_equal(as.matrix(parDist(mat.sample6, method = "soergel"))[1, 2], 0.5)
  expect_equal(as.matrix(parDist(mat.sample7, method = "soergel"))[1, 2], 0.55)
})
# wrong implementation in proxy?
test_that("wave method produces same outputs as dist", {
  expect_equal(as.matrix(parDist(mat.sample1, method = "wave"))[1, 2], NaN)
  expect_equal(as.matrix(parDist(mat.sample2, method = "wave"))[1, 2], NaN)
  expect_equal(as.matrix(parDist(mat.sample3, method = "wave"))[1, 2], 0.5205532370853329)
  expect_equal(as.matrix(parDist(mat.sample4, method = "wave"))[1, 2], NaN)
  expect_equal(as.matrix(parDist(mat.sample5, method = "wave"))[1, 2], -0.0022262504871707334)
  expect_equal(as.matrix(parDist(mat.sample6, method = "wave"))[1, 2], 0.5)
  expect_equal(as.matrix(parDist(mat.sample7, method = "wave"))[1, 2], NaN)
})
# works
test_that("whittaker method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "whittaker")
})
# binary distances
#works
test_that("binary method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "binary")
})
# works
test_that("braunblanquet method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "braun-blanquet")
})
# works
test_that("dice method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "dice")
})
# works
test_that("fager method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "fager")
})
# works
test_that("faith method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "faith")
})
# works
test_that("hamman method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "hamman")
})
# works
test_that("kulczynski1 method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "kulczynski1")
})
# works
test_that("kulczynski2 method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "kulczynski2")
})
# works
test_that("michael method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "michael")
})
# works
test_that("mountford method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "mountford")
})
# works
test_that("mozley method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "mozley")
})
# works
test_that("ochiai method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "ochiai")
})
# works
test_that("phi method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "phi")
})
# works
test_that("russel method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "russel")
})
# works
test_that("simplematching method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "simple matching")
})
# works
test_that("simpson method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "simpson")
})
# works
test_that("stiles method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "stiles")
})
# works
test_that("tanimoto method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "tanimoto")
})
# works
test_that("yule method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "yule")
})
# works
test_that("yule2 method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "yule2")
})
# works
test_that("cosine method produces same outputs as dist", {
  testMatrixListEquality(mat.list[-4], "cosine") # dividing by zero
})
