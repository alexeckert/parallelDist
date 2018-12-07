// DistanceDTWGeneric.h
//
// Copyright (C)  2017, 2018  Alexander Eckert
//
// This file is part of parallelDist.
//
// parallelDist is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// parallelDist is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with parallelDist. If not, see <http://www.gnu.org/licenses/>.

#ifndef DISTANCEDTWGENERIC_H_
#define DISTANCEDTWGENERIC_H_

#include "IDistance.h"
#include <algorithm>
#include <utility>

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) < (y) ? (y) : (x))
#define min3(x, y, z) (min(x, min(y, z)))
#define unused(x) ((void) x)

enum class NormMethod { NoNorm, PathLength, ALength, ABLength };

//==============================
// Dynamic Time Warping distance
//==============================
// Generic implementation
template <typename Implementation>
class DistanceDTWGeneric : public IDistance {
private:
    unsigned int windowSize;
    bool warpingWindow;
    NormMethod normalizationMethod;

    const unsigned int getPatternOffset() {
        return Implementation::patternOffset;
    }

    inline Implementation& impl() {
        return *static_cast<Implementation*>(this);
    }

    /**
     Calculate costs for two entries of input matrices A and B
     @param pen penality matrix
     @param bSizeOffset B.ncol + patternOffset
     @param A matrix A
     @param B matrix B
     @param i index i
     @param j index j
     @return costs for two entries of input matrices A and B
     */
    std::pair<double, int> getCost(double *pen, unsigned int bSizeOffset, const arma::mat &A, const arma::mat &B,
      unsigned int i, unsigned int j) {
        return impl().getCost(pen, bSizeOffset, A, B, i, j);
    }

protected:
  // calculates minimum and argmin of a double array
  std::pair<double, int> argmin(double arr[], unsigned int len) {
    double min = arr[0];
    unsigned int argMin = 0;

    for (unsigned int i = 1; i < len; i++) {
      if (arr[i] < min) {
        min = arr[i];
        argMin = i;
      }
    }
    return std::make_pair(min, argMin);
  }
  // returns value of cell of matrix
  double getCell(double *matrix, unsigned int bMatrixSizeOffset, unsigned int i, unsigned int j) {
    return matrix[i * bMatrixSizeOffset + j];
  }
  // calculates euclidean distance between two matrices
  double getDistance(const arma::mat &A, const arma::mat &B, unsigned int i, unsigned int j) {
    if (i < getPatternOffset() || j < getPatternOffset()) {
      return INFINITY;
    } else {
      return sqrt(arma::accu(arma::square(A.col(i - getPatternOffset()) - B.col(j - getPatternOffset()))));
    }
  }

public:
    DistanceDTWGeneric(bool warpingWindow = false, unsigned int windowSize = 0,
      NormMethod normalizationMethod = NormMethod::NoNorm) {
        this->warpingWindow = warpingWindow;
        this->windowSize = windowSize;
        this->normalizationMethod = normalizationMethod;
    }

    ~DistanceDTWGeneric() {
    }

    double calcDistance(const arma::mat &A, const arma::mat &B) {
        const unsigned int patternOffset = getPatternOffset();
        // vector sizes for convenience
        const unsigned int Asize = A.n_cols, Bsize = B.n_cols;
        const unsigned int aSizeOffset = Asize + patternOffset;
        const unsigned int bSizeOffset = Bsize + patternOffset;

        // size penality matrix according to the possible offset of the steppattern
        double *pen = new double[aSizeOffset * bSizeOffset];

        char *pre = 0;
        if (normalizationMethod == NormMethod::PathLength) {
            pre = new char[aSizeOffset * bSizeOffset];
        } else {
            unused(pre);
        }

        // initialize fully costalty matrix
        for (unsigned int i = 0; i < aSizeOffset; ++i) {
            unsigned int currIdx = i * bSizeOffset;
            for (unsigned int j = 0; j < bSizeOffset; ++j) {
                pen[currIdx + j] = INFINITY;
            }
        }

        // adjust window if needed
        unsigned int effectiveWindowSize;
        if (warpingWindow) {
            effectiveWindowSize = max(windowSize, Asize > Bsize ? Asize - Bsize : Bsize - Asize);
        } else {
            effectiveWindowSize = max(Asize, Bsize);
        }

        for (unsigned int i = patternOffset; i < aSizeOffset; ++i) {
            unsigned int lower = patternOffset, upper = bSizeOffset;

            lower = i > effectiveWindowSize + patternOffset ? i - effectiveWindowSize : patternOffset;
            upper = min(bSizeOffset, i + effectiveWindowSize + 1);

            unsigned int currIdx = i * bSizeOffset;
            for (unsigned int j = lower; j < upper; ++j) {
                if (i == patternOffset && j == patternOffset) {
                    pen[currIdx + j] = getDistance(A, B, i, j);
                } else {
                    std::pair<double, int> cost = getCost(pen, bSizeOffset, A, B, i, j);
                    pen[currIdx + j] = cost.first;

                    if (normalizationMethod == NormMethod::PathLength) {
                        pre[currIdx + j] = cost.second;
                    }
                }
            }
        }

        // remember the optimal distance measure
        double dist = pen[aSizeOffset * bSizeOffset - 1];

        // free memory
        delete [] pen;

        // calc warp path for normalization
        if (normalizationMethod == NormMethod::PathLength) {
            unsigned int warpPathLength = 0;
            // start at the end of the warping path
            unsigned int i = aSizeOffset - 1, j = bSizeOffset - 1;
            unsigned int patternOffsetPlusOne = patternOffset + 1;
            // until starting node reached
            while (i != patternOffset && j != patternOffset) {
                // add node to warping path
                ++warpPathLength;
                if (i == patternOffsetPlusOne) {
                    j -= 1;
                } else if (j == patternOffsetPlusOne) {
                    i -= 1;
                } else if (pre[i * bSizeOffset + j] == 0) {
                    i -= 1;
                } else if (pre[i * bSizeOffset + j] == 1) {
                    i -= 1;
                    j -= 1;
                } else if (pre[i * bSizeOffset + j] == 2) {
                    j -= 1;
                }
            }
            // free memory
            delete [] pre;
            dist /= warpPathLength;
        } else if (normalizationMethod == NormMethod::ABLength) {
            dist /= (Asize + Bsize);
        } else if (normalizationMethod == NormMethod::ALength) {
            dist /= Asize;
        }

        return dist;
    }
};

#endif  // DISTANCEDTWGENERIC_H_
