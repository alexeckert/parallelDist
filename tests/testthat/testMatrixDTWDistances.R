## testMatrixDTWDistances.R
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

mat.sample1 <- matrix(c(0,1,0,1,0,0,1,0), nrow = 2)
mat.sample2 <- matrix(c(0,1,0,1,0,0,1,0,1,1), nrow = 2)
mat.sample3 <- matrix(c(0,1,0,1,0,0,1,0,1,1,0,1), nrow = 3)
mat.sample4 <- matrix(c(1:100), ncol = 5)
mat.sample5 <- matrix(rep(0,100), ncol = 5)
mat.sample6 <- matrix(c(-50:49), ncol = 5)
mat.sample7 <- matrix(c(1:8), ncol = 2)
mat.sample8 <- matrix(c(1:15), ncol = 5)
mat.sample9 <- matrix(c(1,88,55,72,3,33,44,11,3,42, 52,29,3,45,34,21,34,59,35,84), nrow = 2, byrow = TRUE)
tolerance <- 1e-8

mat.list <- list(mat.sample1, mat.sample2, mat.sample3, mat.sample4, mat.sample5, mat.sample6, mat.sample7, mat.sample8, mat.sample9)

if (isCran()) {
  mat.list <- mat.list[1:5]
}

library(dtw)
testMatrixEquality <- function(matrix, method, ...) {
  expect_equal(as.matrix(parDist(matrix, method = method, ...)), as.matrix(dist(matrix, method = method, ...)))
}

testMatrixListEquality <- function(matlist, method, ...) {
  invisible(sapply(matlist, function(x) { testMatrixEquality(x, method, ...) }))
}

test_that("error for unsupported step pattern shows up", {
  expect_error(parDist(mat.sample1, method = "dtw", step.pattern="unknown"), "Step pattern is not supported.")
})

# symmetric1 / symmetric
testMatrixListEquality(mat.list, method="dtw", window.type="none", step.pattern=symmetric1)
invisible(sapply(mat.list, function(x) {
  expect_equal(as.matrix(parDist(x, method = "dtw",  window.type="none", step.pattern="symmetric1")),
               as.matrix(dist(x, method = "dtw",  step.pattern=symmetric1)))
})) # step.pattern as string
testMatrixListEquality(mat.list, method="dtw", window.size = 5, window.type=sakoeChibaWindow, step.pattern=symmetric1)

# symmetric2 / symmetricP0
testMatrixListEquality(mat.list, method="dtw", window.type="none", step.pattern=symmetric2)
testMatrixListEquality(mat.list, method="dtw", window.size = 5,  window.type=sakoeChibaWindow, step.pattern=symmetric2)

# symmetricP05
testMatrixListEquality(mat.list, method="dtw",  window.type="none", step.pattern=symmetricP05)
testMatrixListEquality(mat.list, method="dtw", window.size = 5,  window.type=sakoeChibaWindow, step.pattern=symmetricP05)

# symmetricP1
testMatrixListEquality(mat.list, method="dtw", window.type="none", step.pattern=symmetricP1)
testMatrixListEquality(mat.list, method="dtw", window.size = 5, window.type=sakoeChibaWindow, step.pattern=symmetricP1)

# symmetricP2
testMatrixListEquality(mat.list, method="dtw", window.type="none", step.pattern=symmetricP2)
testMatrixListEquality(mat.list, method="dtw", window.size = 5, window.type=sakoeChibaWindow, step.pattern=symmetricP2)

# asymmetric
testMatrixListEquality(mat.list, method="dtw", window.type="none", step.pattern=asymmetric)
testMatrixListEquality(mat.list, method="dtw", window.size = 5, window.type=sakoeChibaWindow, step.pattern=asymmetric)

# asymmetricP0
testMatrixListEquality(mat.list, method="dtw", window.type="none", step.pattern=asymmetricP0)
testMatrixListEquality(mat.list, method="dtw", window.size = 5, window.type=sakoeChibaWindow, step.pattern=asymmetricP0)

# asymmetricP05
testMatrixListEquality(mat.list, method="dtw", window.type="none", step.pattern=asymmetricP05)
testMatrixListEquality(mat.list, method="dtw", window.size = 5, window.type=sakoeChibaWindow, step.pattern=asymmetricP05)

# asymmetricP1
testMatrixListEquality(mat.list, method="dtw", window.type="none", step.pattern=asymmetricP1)
testMatrixListEquality(mat.list, method="dtw", window.size = 5, window.type=sakoeChibaWindow, step.pattern=asymmetricP1)

# asymmetricP2
testMatrixListEquality(mat.list, method="dtw", window.type="none", step.pattern=asymmetricP2)
testMatrixListEquality(mat.list, method="dtw", window.size = 5, window.type=sakoeChibaWindow, step.pattern=asymmetricP2)

# warping path normalization
expect_equal(as.matrix(parDist(mat.sample9[1:2, 1:3], method = "dtw", norm.method="path.length", threads=1)),
             as.matrix(dist(mat.sample9[1:2, 1:3], method = "dtw", step.pattern=symmetric1)) / 3,
             tolerance=tolerance)
expect_equal(as.matrix(parDist(mat.sample3, method = "dtw", norm.method="path.length", threads=1)),
             as.matrix(dist(mat.sample3, method = "dtw", step.pattern=symmetric1)) / 5,
             tolerance=tolerance)
expect_equal(as.matrix(parDist(mat.sample4[c(3, 20),], method = "dtw", norm.method="path.length", threads=1)),
             as.matrix(dist(mat.sample4[c(3, 20),], method = "dtw", step.pattern=symmetric1)) / 7, tolerance=tolerance)
expect_equal(as.matrix(parDist(mat.sample9, method = "dtw", norm.method="n")),
             as.matrix(dist(mat.sample9, method = "dtw", step.pattern=symmetric1)) / dim(mat.sample9)[2],
             tolerance=tolerance)
expect_equal(as.matrix(parDist(mat.sample9, method = "dtw", norm.method="n+m")),
             as.matrix(dist(mat.sample9, method = "dtw", step.pattern=symmetric1)) / (dim(mat.sample9)[2] * 2),
             tolerance=tolerance)
