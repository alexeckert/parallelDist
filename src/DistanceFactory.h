// DistanceFactory.h
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

#ifndef DISTANCEFACTORY_H_
#define DISTANCEFACTORY_H_

#include "IDistance.h"
#include <memory>
#include <vector>

//==============================
// Distance Factory
//==============================
class DistanceFactory {
private:
    // Store reference to data objects to enable distance method parameter precalculations
    arma::mat dataMatrix;
    std::vector<arma::mat> dataMatrixList;
    bool isDataMatrix;
public:
    explicit DistanceFactory(arma::mat &dataMatrix) : dataMatrix(dataMatrix), isDataMatrix(true) {}
    explicit DistanceFactory(std::vector<arma::mat> &dataMatrixList) : dataMatrixList(dataMatrixList),
    isDataMatrix(false) {}
    std::shared_ptr<IDistance> createDistanceFunction(const Rcpp::List& attrs, const Rcpp::List& arguments);
};

#endif  // DISTANCEFACTORY_H_
