## testMatrixListDistances.R
##
## Copyright (C)  2017, 2021  Alexander Eckert
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

context("Distance methods using matrix list as input")

matToList <- function(matrix, vert=TRUE) {
  mat <- if (vert == TRUE) {
    lapply(as.list(data.frame(t(matrix))), function(x) { matrix(x)})
  } else {
    lapply(as.list(data.frame(t(matrix))), function(x) {matrix(x, nrow=1)})
  }
  names(mat) <- NULL
  mat
}

mat.sample1 <- matrix(c(0,1,0,1,0,0,1,0), nrow = 2)
mat.sample2 <- matrix(c(0,1,0,1,0,0,1,0,1,1), nrow = 2)
mat.sample3 <- matrix(c(1:500), ncol = 5)
mat.sample4 <- matrix(rep(0,100), ncol = 5)
mat.sample5 <- matrix(c(-500:499), ncol = 5)
mat.sample6 <- matrix(c(1:2), ncol = 1)
mat.sample7 <- matrix(c(0.5,1,0,1,0,0,1,0.3,1,1), nrow = 2)

mat.list <- list(mat.sample1, mat.sample2, mat.sample3, mat.sample4, mat.sample5, mat.sample6, mat.sample7)
matlist.list.h <- lapply(mat.list, function(x) matToList(x, vert=FALSE))

if (isCran()) {
  mat.list <- mat.list[1:4]
  matlist.list.h <- lapply(mat.list, function(x) matToList(x, vert=FALSE))
}

testMatrixListMatrixEquality <- function(matList, matrix, method, ...) {
  suppressNaNWarnings(expect_equal(as.matrix(parDist(matList, method = method, ...)), as.matrix(dist(matrix, method = method, ...))))
}

testMatrixListEquality <- function(matListList, matList, method, ...) {
  mapply(function(X,Y) {
    testMatrixListMatrixEquality(X, Y, method, ...)
  }, X=matListList, Y=matList)
}

library(proxy)

test_that("bhjattacharyya method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "bhjattacharyya")
})

test_that("bray method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "bray")
})

test_that("canberra method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "canberra")
})
# first row only
test_that("chord method produces same outputs as dist", {
  expect_warning(testMatrixListEquality(matlist.list.h, mat.list, "chord"))
})

test_that("divergence method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "divergence")
})

test_that("euclidean method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "euclidean")
})

test_that("fJaccard method produces same outputs as dist", {
  testMatrixListEquality(list(matToList(mat.sample1, vert=FALSE), matToList(mat.sample2, vert=FALSE), matToList(mat.sample7, vert=FALSE)),
                         list(mat.sample1, mat.sample2, mat.sample7),
                         "fJaccard")
})
# first row only
test_that("geodesic method produces same outputs as dist", {
  expect_warning(testMatrixListEquality(matlist.list.h, mat.list, "geodesic"))
})

test_that("hellinger method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "hellinger")
})

test_that("kullback method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "kullback")
})

test_that("manhattan method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "manhattan")
})

test_that("maximum method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "maximum")
})

test_that("minkowski method produces same outputs as dist", {
  for (p_param in c(1:10)) {
    testMatrixListEquality(matlist.list.h, mat.list, "minkowski", p=as.numeric(p_param))
  }
})
# first row only
test_that("podani method produces same outputs as dist", {
  expect_warning(testMatrixListEquality(matlist.list.h, mat.list, "podani"))
})
# wrong implementation in proxy?
#test_that("soergel method produces same outputs as dist", {
#  testMatrixListEquality(matlist.list.h, mat.list, "soergel")
#})
# wrong implementation in proxy?
#test_that("wave method produces same outputs as dist", {
#  testMatrixListEquality(matlist.list.h, mat.list, "wave")
#})

test_that("whittaker method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "whittaker")
})

# binary distances
test_that("binary method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "binary")
})

test_that("braunblanquet method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "braun-blanquet")
})

test_that("dice method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "dice")
})

test_that("fager method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "fager")
})

test_that("faith method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "faith")
})

test_that("hamman method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "hamman")
})

test_that("kulczynski1 method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "kulczynski1")
})

test_that("kulczynski2 method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "kulczynski2")
})

test_that("michael method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "michael")
})

test_that("mountford method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "mountford")
})

test_that("mozley method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "mozley")
})

test_that("ochiai method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "ochiai")
})

test_that("phi method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "phi")
})

test_that("russel method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "russel")
})

test_that("simplematching method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "simple matching")
})

test_that("simpson method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "simpson")
})

test_that("stiles method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "stiles")
})

test_that("tanimoto method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "tanimoto")
})

test_that("yule method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "yule")
})

test_that("yule2 method produces same outputs as dist", {
  testMatrixListEquality(matlist.list.h, mat.list, "yule2")
})
