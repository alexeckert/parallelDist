context("Distance methods using matrix as input")

set.seed(1)
mat.sample1 <- matrix(c(NA, 2, NA, 4, NA, 6, 7, 8 ,9, NA, 11, 12), nrow = 4)
mat.sample2 <- matrix(c(0, 1, 0, 1, NA, 0, 1, 0, 1, 1), nrow = 2)
mat.sample3 <- matrix(c(1:500), ncol = 5)
p <- 0.3
mat.sample3 <- apply(mat.sample3, 1:2, function(x) sample(c(x, NA), 1, prob=c((1 - p), p)))
mat.sample4 <- matrix(rep(0, 100), ncol = 5)
p <- 0.4
mat.sample4 <- apply(mat.sample4, 1:2, function(x) sample(c(x, NA), 1, prob=c((1 - p), p)))
mat.sample5 <- matrix(c(-500:499), ncol = 5)
p <- 0.2
mat.sample5 <- apply(mat.sample5, 1:2, function(x) sample(c(x, NA), 1, prob=c((1 - p), p)))
mat.sample6 <- matrix(c(1, NA), ncol = 1)
mat.sample7 <- matrix(NA, ncol = 1)

mat.list <- list(mat.sample1, mat.sample2, mat.sample3, mat.sample4, mat.sample5, mat.sample6, mat.sample7)

if (isCran()) {
  mat.list <- mat.list[1:4]
}

testMatrixEquality <- function(matrix, method, ...) {
  suppressNaNWarnings(expect_equal(as.matrix(parDist(matrix, method = method, ...)), as.matrix(dist(matrix, method = method, ...))))
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

test_that("dist class label attribute keeps preserved", {
  namedMatrix <- matrix(1:12, 4)
  colnames(namedMatrix) <- c("A", "B", "C")
  rownames(namedMatrix) <- c("a", "b", "c", "d")

  attributes1 <- attributes(parDist(namedMatrix))
  attributes2 <- attributes(dist(namedMatrix))

  expect_equal(attributes1$Labels, attributes2$Labels)
})

test_that("canberra method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "canberra")
})

test_that("euclidean method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "euclidean")
})


test_that("manhattan method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "manhattan")
})

test_that("maximum method produces same outputs as dist", {
  testMatrixListEquality(mat.list, "maximum")
})

test_that("minkowski method produces same outputs as dist", {
  for (p_param in c(1:10)) {
    testMatrixListEquality(mat.list, "minkowski", p = as.numeric(p_param))
  }
})


