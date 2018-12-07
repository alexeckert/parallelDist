// parallelDist.cpp
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

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <cmath>
#include <list>
#include <vector>
#include <algorithm>
#include <iostream>

#include "IDistance.h"
#include "DistanceFactory.h"


uint64_t sumForm(const uint64_t n) {
    return (pow(n, 2) + n) / 2;
}

uint64_t matToVecIdx(const uint64_t i, const uint64_t j, const uint64_t N) {
    return i * N - i - sumForm(i) - 1 + j;
}

inline bool isInteger(const std::string& s) {
    if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false;
    char * p;
    strtol(s.c_str(), &p, 10);
    return (*p == 0);
}

struct DistanceVec : public RcppParallel::Worker {
    // input vector of matrices
    const std::vector<arma::mat>& seriesVec;

    int vecSize = 0;

    // output vector by reference
    Rcpp::NumericVector& rvec;

    // distance function
    std::shared_ptr<IDistance> distance;

    // initialize from Rcpp input and output matrixes (the RMatrix class
    // can be automatically converted to from the Rcpp matrix type)
    DistanceVec(const std::vector<arma::mat>& seriesVec, Rcpp::NumericVector& rvec,
      const std::shared_ptr<IDistance>& distance)
        : seriesVec(seriesVec), rvec(rvec), distance(distance) {
        vecSize = seriesVec.size();
    }

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            for (std::size_t j = 0; j < i; j++) {
                rvec[matToVecIdx(j, i, vecSize)] = distance->calcDistance(seriesVec.at(i), seriesVec.at(j));
            }
        }
    }
};

// uses not a list but the matrix
struct DistanceMatrixVec : public RcppParallel::Worker {
    // input vector of matrices
    const arma::mat& seriesVec;

    int vecSize = 0;

    // output vector by reference
    Rcpp::NumericVector& rvec;

    // distance function
    std::shared_ptr<IDistance> distance;

    // initialize from Rcpp input and output matrixes (the RMatrix class
    // can be automatically converted to from the Rcpp matrix type)
    DistanceMatrixVec(const arma::mat& seriesVec, Rcpp::NumericVector& rvec, const std::shared_ptr<IDistance>& distance)
        : seriesVec(seriesVec), rvec(rvec), distance(distance) {
        vecSize = seriesVec.n_rows;
    }

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            for (std::size_t j = 0; j < i; j++) {
                rvec[matToVecIdx(j, i, vecSize)] = distance->calcDistance(seriesVec.row(i), seriesVec.row(j));
            }
        }
    }
};

void setVectorAttributes(Rcpp::NumericVector &rvec, const Rcpp::List &attrs) {
    rvec.attr("Size") = attrs["Size"];
    rvec.attr("Labels") = attrs["Labels"];
    rvec.attr("Diag") = Rcpp::as<bool >(attrs["Diag"]);
    rvec.attr("Upper") = Rcpp::as<bool >(attrs["Upper"]);
    rvec.attr("method") = attrs["method"];
    rvec.attr("call") = attrs["call"];
    rvec.attr("class") = "dist";
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_parallelDistVec(Rcpp::List dataList, Rcpp::List attrs, Rcpp::List arguments) {
    uint64_t n = dataList.size();
    // result matrix
    Rcpp::NumericVector rvec(sumForm(n) - n);

    setVectorAttributes(rvec, attrs);

    // Convert list to vector of double matrices
    std::vector<arma::mat> listVec;
    for (Rcpp::List::iterator it = dataList.begin(); it != dataList.end(); ++it) {
        listVec.push_back(Rcpp::as<arma::mat >(*it));
    }
    std::shared_ptr<IDistance> distanceFunction = DistanceFactory(listVec).createDistanceFunction(attrs, arguments);

    DistanceVec* distanceWorker = new DistanceVec(listVec, rvec, distanceFunction);
    // call it with parallelFor
    RcppParallel::parallelFor(0, n, (*distanceWorker));
    delete distanceWorker;
    distanceWorker = NULL;

    return rvec;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_parallelDistMatrixVec(arma::mat dataMatrix, Rcpp::List attrs, Rcpp::List arguments) {
    uint64_t n = dataMatrix.n_rows;

    // result matrix
    Rcpp::NumericVector rvec(sumForm(n) - n);

    setVectorAttributes(rvec, attrs);

    std::shared_ptr<IDistance> distanceFunction = DistanceFactory(dataMatrix).createDistanceFunction(attrs, arguments);
    DistanceMatrixVec* distanceWorker = new DistanceMatrixVec(dataMatrix, rvec, distanceFunction);
    // call it with parallelFor
    RcppParallel::parallelFor(0, n, (*distanceWorker));
    delete distanceWorker;
    distanceWorker = NULL;

    return rvec;
}
