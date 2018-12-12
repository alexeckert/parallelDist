// StepPattern.h
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

#ifndef STEPPATTERN_H_
#define STEPPATTERN_H_

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include "DistanceDTWGeneric.h"

// StepPatternSymmetric1
//   g[i,j] = min(
//     g[i-1,j-1] +     d[i  ,j  ] ,
//     g[i  ,j-1] +     d[i  ,j  ] ,
//     g[i-1,j  ] +     d[i  ,j  ] ,
//   )
class StepPatternSymmetric1 : public DistanceDTWGeneric<StepPatternSymmetric1> {
public:
    using DistanceDTWGeneric::DistanceDTWGeneric;
    static constexpr unsigned int patternOffset = 1;
    std::pair<double, int> getCost(double *pen, unsigned int bSizeOffset, const arma::mat &A, const arma::mat &B,
      unsigned int i, unsigned int j) {
        double distance = getDistance(A, B, i, j);
        // bound checking for this step pattern
        double minArray[3] = {
            getCell(pen, bSizeOffset, i-1, j-1),
            getCell(pen, bSizeOffset, i, j-1),
            getCell(pen, bSizeOffset, i-1, j)
        };

        std::pair<double, int> result = argmin(minArray, 3);
        result.first += distance;
        return result;
    }
};

// StepPatternSymmetric2
// g[i,j] = min(
//   g[i-1,j-1] + 2 * d[i  ,j  ] ,
//   g[i  ,j-1] +     d[i  ,j  ] ,
//   g[i-1,j  ] +     d[i  ,j  ] ,
// )
class StepPatternSymmetric2 : public DistanceDTWGeneric<StepPatternSymmetric2> {
public:
    using DistanceDTWGeneric::DistanceDTWGeneric;
    static constexpr unsigned int patternOffset = 1;
    std::pair<double, int> getCost(double *pen, unsigned int bSizeOffset, const arma::mat &A, const arma::mat &B,
      unsigned int i, unsigned int j) {
        double distance = getDistance(A, B, i, j);
        // bound checking for this step pattern
        double minArray[3] = {
            2 * distance + getCell(pen, bSizeOffset, i-1, j-1),
            distance + getCell(pen, bSizeOffset, i, j-1),
            distance + getCell(pen, bSizeOffset, i-1, j)
        };
        return argmin(minArray, 3);
    }
};

// StepPatternAsymmetric
// g[i,j] = min(
//   g[i-1,j  ] +     d[i  ,j  ] ,
//   g[i-1,j-1] +     d[i  ,j  ] ,
//   g[i-1,j-2] +     d[i  ,j  ] ,
// )
class StepPatternAsymmetric : public DistanceDTWGeneric<StepPatternAsymmetric> {
public:
    using DistanceDTWGeneric::DistanceDTWGeneric;
    static constexpr unsigned int patternOffset = 2;
    std::pair<double, int> getCost(double *pen, unsigned int bSizeOffset, const arma::mat &A, const arma::mat &B,
      unsigned int i, unsigned int j) {
        double distance = getDistance(A, B, i, j);
        // bound checking for this step pattern
        double minArray[3] = {
            distance + getCell(pen, bSizeOffset, i-1, j),
            distance + getCell(pen, bSizeOffset, i-1, j-1),
            distance + getCell(pen, bSizeOffset, i-1, j-2)
        };

        return argmin(minArray, 3);
    }
};

// StepPatternAsymmetricP0
// g[i,j] = min(
//   g[i  ,j-1] + 0 * d[i  ,j  ] ,
//   g[i-1,j-1] +     d[i  ,j  ] ,
//   g[i-1,j  ] +     d[i  ,j  ] ,
// )
class StepPatternAsymmetricP0 : public DistanceDTWGeneric<StepPatternAsymmetricP0> {
public:
    using DistanceDTWGeneric::DistanceDTWGeneric;
    static constexpr unsigned int patternOffset = 1;
    std::pair<double, int> getCost(double *pen, unsigned int bSizeOffset, const arma::mat &A, const arma::mat &B,
      unsigned int i, unsigned int j) {
        double distance = getDistance(A, B, i, j);
        // bound checking for this step pattern
        double minArray[3] = {
            getCell(pen, bSizeOffset, i, j-1),
            distance + getCell(pen, bSizeOffset, i-1, j-1),
            distance + getCell(pen, bSizeOffset, i-1, j)
        };

        return argmin(minArray, 3);
    }
};

// StepPatternAsymmetricP05
// g[i,j] = min(
//      g[i-1,j-3] +0.33 * d[i  ,j-2] +0.33 * d[i  ,j-1] +0.33 * d[i  ,j  ] ,
//      g[i-1,j-2] +0.5 * d[i  ,j-1] +0.5 * d[i  ,j  ] ,
//      g[i-1,j-1] +     d[i  ,j  ] ,
//      g[i-2,j-1] +     d[i-1,j  ] +     d[i  ,j  ] ,
//      g[i-3,j-1] +     d[i-2,j  ] +     d[i-1,j  ] +     d[i  ,j  ] ,
//   )
class StepPatternAsymmetricP05 : public DistanceDTWGeneric<StepPatternAsymmetricP05> {
public:
    using DistanceDTWGeneric::DistanceDTWGeneric;
    static constexpr unsigned int patternOffset = 3;
    std::pair<double, int> getCost(double *pen, unsigned int bSizeOffset, const arma::mat &A, const arma::mat &B,
      unsigned int i, unsigned int j) {
        double multiplier = 1.0 / 3;
        double minArray[5] = {
            getCell(pen, bSizeOffset, i-1, j-3) + multiplier * getDistance(A, B, i, j-2) +
              multiplier * getDistance(A, B, i, j-1) + multiplier * getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-1, j-2) + 0.5 * getDistance(A, B, i, j-1) + 0.5 * getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-1, j-1) + getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-2, j-1) + getDistance(A, B, i-1, j) + getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-3, j-1) + getDistance(A, B, i-2, j) + getDistance(A, B, i-1, j) +
              getDistance(A, B, i, j)
        };
        return argmin(minArray, 5);
    }
};

//
// g[i,j] = min(
//   g[i-1,j-3] + 2 * d[i  ,j-2] +     d[i  ,j-1] +     d[i  ,j  ] ,
//   g[i-1,j-2] + 2 * d[i  ,j-1] +     d[i  ,j  ] ,
//   g[i-1,j-1] + 2 * d[i  ,j  ] ,
//   g[i-2,j-1] + 2 * d[i-1,j  ] +     d[i  ,j  ] ,
//   g[i-3,j-1] + 2 * d[i-2,j  ] +     d[i-1,j  ] +     d[i  ,j  ] ,
// )
class StepPatternSymmetricP05 : public DistanceDTWGeneric<StepPatternSymmetricP05> {
public:
    using DistanceDTWGeneric::DistanceDTWGeneric;
    static constexpr unsigned int patternOffset = 3;
    std::pair<double, int> getCost(double *pen, unsigned int bSizeOffset, const arma::mat &A, const arma::mat &B,
      unsigned int i, unsigned int j) {
        double minArray[5] = {
            getCell(pen, bSizeOffset, i-1, j-3) + 2 * getDistance(A, B, i, j-2) + getDistance(A, B, i, j-1) +
              getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-1, j-2) + 2 * getDistance(A, B, i, j-1) + getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-1, j-1) + 2 * getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-2, j-1) + 2 * getDistance(A, B, i-1, j) + getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-3, j-1) + 2 * getDistance(A, B, i-2, j) + getDistance(A, B, i-1, j) +
              getDistance(A, B, i, j)
        };
        return argmin(minArray, 5);
    }
};

class StepPatternSymmetricP1 : public DistanceDTWGeneric<StepPatternSymmetricP1> {
public:
    using DistanceDTWGeneric::DistanceDTWGeneric;
    static constexpr unsigned int patternOffset = 2;
    std::pair<double, int> getCost(double *pen, unsigned int bSizeOffset, const arma::mat &A, const arma::mat &B,
      unsigned int i, unsigned int j) {
        double minArray[3] = {
            getCell(pen, bSizeOffset, i-1, j-2) + 2 * getDistance(A, B, i, j-1) + getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-1, j-1) + 2 * getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-2, j-1) + 2 * getDistance(A, B, i-1, j) + getDistance(A, B, i, j)
        };
        return argmin(minArray, 3);
    }
};

class StepPatternAsymmetricP1 : public DistanceDTWGeneric<StepPatternAsymmetricP1> {
public:
    using DistanceDTWGeneric::DistanceDTWGeneric;
    static constexpr unsigned int patternOffset = 2;
    std::pair<double, int> getCost(double *pen, unsigned int bSizeOffset, const arma::mat &A, const arma::mat &B,
      unsigned int i, unsigned int j) {
        double minArray[3] = {
            getCell(pen, bSizeOffset, i-1, j-2) + 0.5 * getDistance(A, B, i, j-1) + 0.5 * getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-1, j-1) + getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-2, j-1) + getDistance(A, B, i-1, j) + getDistance(A, B, i, j)
        };
        return argmin(minArray, 3);
    }
};


// StepPatternAsymmetricP2
//   g[i,j] = min(
//     g[i-2,j-3] +0.67 * d[i-1,j-2] +0.67 * d[i  ,j-1] +0.67 * d[i  ,j  ] ,
//     g[i-1,j-1] +     d[i  ,j  ] ,
//     g[i-3,j-2] +     d[i-2,j-1] +     d[i-1,j  ] +     d[i  ,j  ] ,
//   )
class StepPatternAsymmetricP2 : public DistanceDTWGeneric<StepPatternAsymmetricP2> {
public:
    using DistanceDTWGeneric::DistanceDTWGeneric;
    static constexpr unsigned int patternOffset = 3;
    std::pair<double, int> getCost(double *pen, unsigned int bSizeOffset, const arma::mat &A, const arma::mat &B,
      unsigned int i, unsigned int j) {
        double multiplier = 2.0 / 3;
        double minArray[3] = {
            getCell(pen, bSizeOffset, i-2, j-3) + multiplier * getDistance(A, B, i-1, j-2) +
              multiplier * getDistance(A, B, i, j-1) + multiplier * getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-1, j-1) + getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-3, j-2) + getDistance(A, B, i-2, j-1) + getDistance(A, B, i-1, j) +
              getDistance(A, B, i, j)
        };
        return argmin(minArray, 3);
    }
};

// StepPatternSymmetricP2
// g[i,j] = min(
//   g[i-2,j-3] + 2 * d[i-1,j-2] + 2 * d[i  ,j-1] +     d[i  ,j  ] ,
//   g[i-1,j-1] + 2 * d[i  ,j  ] ,
//   g[i-3,j-2] + 2 * d[i-2,j-1] + 2 * d[i-1,j  ] +     d[i  ,j  ] ,
// )
class StepPatternSymmetricP2 : public DistanceDTWGeneric<StepPatternSymmetricP2> {
public:
    using DistanceDTWGeneric::DistanceDTWGeneric;
    static constexpr unsigned int patternOffset = 3;
    std::pair<double, int> getCost(double *pen, unsigned int bSizeOffset, const arma::mat &A, const arma::mat &B,
      unsigned int i, unsigned int j) {
        double minArray[3] = {
            getCell(pen, bSizeOffset, i-2, j-3) + 2 * getDistance(A, B, i-1, j-2) + 2 * getDistance(A, B, i, j-1) +
              getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-1, j-1) + 2 * getDistance(A, B, i, j),
            getCell(pen, bSizeOffset, i-3, j-2) + 2 * getDistance(A, B, i-2, j-1) + 2* getDistance(A, B, i-1, j) +
              getDistance(A, B, i, j)
        };
        return argmin(minArray, 3);
    }
};

#endif  // STEPPATTERN_H_
