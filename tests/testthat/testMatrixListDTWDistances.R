## testMatrixListDTWDistances.R
##
## Copyright (C)  2018, 2021  Alexander Eckert
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

context("DTW distance methods using matrix as input")

# Sample data

matList.sample1 <- list(structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1), .Dim = c(1L, 40L)),
           structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), .Dim = c(1L, 26L)),
           structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1), .Dim = c(1L, 18L)), structure(1, .Dim = c(1L, 1L)))

matList.sample2 <- list(structure(c(0, 0.666221290024121, 0.992138579223845, 0.837675741488087,
                       0.298093968979853, -0.36598027656586, -0.861963396802543, -0.992827766724864,
                       -0.725930890382089, -0.185671194528465, 0.413176851770962, 0.855075968430664,
                       0.999129796868121, 0.816424749017454, 0.384112677013609, -0.153156594492941,
                       -0.635636580058955, -0.935956802538981, -0.988865502668731, -0.798723238550414,
                       -0.427611028978791, 0.0285324490339255, 0.467078394773374, 0.801915060528429,
                       0.978294933963238, 0.978050536506311, 0.816424749017454, 0.533244498772231,
                       0.18159436371588, -0.183237680976828, -0.512615351168829, -0.770254350539228,
                       -0.934664583522371, -0.998831629222526, -0.967951368250845, -0.856101340405466,
                       -0.682627878063188, -0.46882893202083, -0.235286772144856, -7.34763812293426e-16
), .Dim = c(1L, 40L)), structure(c(0, 0.637294298455617, 0.98319606956513,
                                   0.868860740689089, 0.339292108143587, -0.35771316353353, -0.883391258820254,
                                   -0.972575815035646, -0.571496467513928, 0.124331294650235, 0.759722381112957,
                                   0.999367384948098, 0.707879926506996, 0.032504191469351, -0.664423623249325,
                                   -0.99683424992037, -0.771156544793101, -0.105476926857648, 0.624552817559215,
                                   0.993872910172706, 0.777587097879622, 0.0946947548310227, -0.649096030626008,
                                   -0.998323195253652, -0.729083830715588, -7.34763812293426e-16
), .Dim = c(1L, 26L)), structure(c(0, 0.246969768376824, 0.508644727112303,
                                   0.752713514852865, 0.933926624849675, 0.999757468870668, 0.902994631700102,
                                   0.621064698874997, 0.178051104865969, -0.339697691026722, -0.787475688668867,
                                   -0.99758720536247, -0.847446820461593, -0.340578183393838, 0.343461298353224,
                                   0.883479231147892, 0.961169928913234, 0.470238319859503, -0.339697691026722,
                                   -0.938153498032876, -0.851671382567693, -0.0794079161883605,
                                   0.780014632342609, 0.945203224149745, 0.193789516277419, -0.772441004193791,
                                   -0.911415695442887, -7.34763812293426e-16), .Dim = c(1L, 28L)),
structure(c(0, 0.498706765350665, 0.875807932165434, 0.997939444221616,
            0.793049199327792, 0.297537649041793, -0.328267139499277,
            -0.838183694383601, -0.994438048115168, -0.688929896390012,
            -0.0318455120202145, 0.658786647043114, 0.99698746548441,
            0.752329817086689, 0.0318455120202155, -0.72482797811622,
            -0.994438048115168, -0.545388021934346, 0.328267139499275,
            0.954710085524754, 0.793049199327793, -0.0641628059443414,
            -0.875807932165433, -0.866770766808318, -7.34763812293426e-16
), .Dim = c(1L, 25L)))


matListList <- list(matList.sample1, matList.sample2)

# Test methods

# Method for calculating a dtw distance matrix with vectors of different length
calcDTWDistMatForMatList <- function(matrixList, step.pattern) {
  matrixListLength <- length(matrixList)
  mat <- matrix(NA, nrow=matrixListLength, ncol=matrixListLength)
  for (i in c(1:(matrixListLength -1))) {
    for (j in (i+1):matrixListLength) {
      mat[j,i] <- tryCatch({
        dtw(as.vector(matrixList[[j]]), as.vector(matrixList[[i]]), step.pattern = step.pattern, distance.only = TRUE)$distance
      }, error=function(e){
        Inf
      })
    }
  }
  as.dist(mat)
}

library(dtw)
testDtwMatrixListEquality <- function(matList, threads = 2, ...) {
  expect_equal(as.matrix(calcDTWDistMatForMatList(matList, ...)),
               as.matrix(parDist(matList, method = "dtw", threads = threads, window.type="none", ...)))
}

# Tests

test_that("parDist produces same distance for list of matrices with different length for different step patterns", {
  threadsForTest = 2
  for (matList in matListList) {
    testDtwMatrixListEquality(matList, threads = threadsForTest, step.pattern = symmetric1)
    testDtwMatrixListEquality(matList, threads = threadsForTest, step.pattern = symmetric2)
    testDtwMatrixListEquality(matList, threads = threadsForTest, step.pattern = symmetricP05)
    testDtwMatrixListEquality(matList, threads = threadsForTest, step.pattern = symmetricP1)
    testDtwMatrixListEquality(matList, threads = threadsForTest, step.pattern = symmetricP2)
    testDtwMatrixListEquality(matList, threads = threadsForTest, step.pattern = asymmetric)
    testDtwMatrixListEquality(matList, threads = threadsForTest, step.pattern = asymmetricP0)
    testDtwMatrixListEquality(matList, threads = threadsForTest, step.pattern = asymmetricP05)
    testDtwMatrixListEquality(matList, threads = threadsForTest, step.pattern = asymmetricP1)
    testDtwMatrixListEquality(matList, threads = threadsForTest, step.pattern = asymmetricP2)
  }
})
