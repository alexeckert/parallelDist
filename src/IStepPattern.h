// IStepPattern.h
//
// Copyright (C)  2017  Alexander Eckert
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

#ifndef ISTEPPATTERN_H_
#define ISTEPPATTERN_H_

#include <iostream>
#include <string>
#include <math.h>
#include <RcppArmadillo.h>

template <typename Implementation>
class IStepPattern {
private:
    unsigned int matrixSize;
    unsigned int idxOffset;
protected:
    // returns value of cell of matrix
    double getCell(double *matrix, unsigned int i, unsigned int j) {
        return matrix[i * idxOffset + j];
    }
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
    inline Implementation& impl() {
        return *static_cast<Implementation*>(this);
    }
public:
    void setMatrixSize(unsigned int matrixSize) {
        this->matrixSize = matrixSize;
        this->idxOffset = matrixSize + getPatternOffset();
    }
    double getPatternOffset() {
        return impl().patternOffset;
    }
    double getDistance(const arma::mat &A, const arma::mat &B, unsigned int i, unsigned int j) {

        if (i < getPatternOffset() || j < getPatternOffset()) {
            //  std::cout << "distI " << i << " distJ " << j << std::endl;
            //  std::cout << "inf" << std::endl;
            return INFINITY;
        } else {
            //return sqrt(arma::accu(arma::pow(A.col(i - patternOffset) - B.col(j - patternOffset), 2)));
            return arma::accu(arma::abs(A.col(i - getPatternOffset()) - B.col(j - getPatternOffset())));
        }
    }
    // constructor
    IStepPattern() : matrixSize(0), idxOffset(0) {}
    ~IStepPattern() {}

    std::pair<double, int> getCost(double *pen, const arma::mat &A, const arma::mat &B, unsigned int i, unsigned int j) {
        return impl().getCost(pen, A, B, i, j);
    };
};


// method for printing a 2dArray
// template <size_t rows, size_t cols>
// void print2dArray(double (&array)[rows][cols])
// {
//   std::cout << __func__ << std::endl;
//   for (size_t i = 0; i < rows; ++i)
//   {
//     std::cout << i << ": ";
//     for (size_t j = 0; j < cols; ++j)
//       std::cout << array[i][j] << '\t';
//     std::cout << std::endl;
//   }
// }

#endif
