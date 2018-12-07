// DistanceFactory.cpp
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

#include "DistanceFactory.h"
#include "DistanceDTWFactory.h"
#include "Util.h"
#include "DistanceDist.h"
#include "DistanceBinary.h"

std::shared_ptr<IDistance> DistanceFactory::createDistanceFunction(const Rcpp::List& attrs,
  const Rcpp::List& arguments) {
    using util::isEqualStr;
    std::string distName = attrs["method"];
    std::shared_ptr<IDistance> distanceFunction = NULL;

    if (isEqualStr(distName, "bhjattacharyya")) {
        distanceFunction = std::make_shared<DistanceBhjattacharyya>();
    } else if (isEqualStr(distName, "bray")) {
        distanceFunction = std::make_shared<DistanceBray>();
    } else if (isEqualStr(distName, "canberra")) {
        distanceFunction = std::make_shared<DistanceCanberra>();
    } else if (isEqualStr(distName, "chord")) {
        distanceFunction = std::make_shared<DistanceChord>();
    } else if (isEqualStr(distName, "divergence")) {
        distanceFunction = std::make_shared<DistanceDivergence>();
    } else if (isEqualStr(distName, "dtw")) {
        distanceFunction = DistanceDTWFactory().createDistanceFunction(distName, arguments);
    } else if (isEqualStr(distName, "fJaccard")) {
        distanceFunction = std::make_shared<DistanceFJaccard>();
    } else if (isEqualStr(distName, "geodesic")) {
        distanceFunction = std::make_shared<DistanceGeodesic>();
    } else if (isEqualStr(distName, "hellinger")) {
        distanceFunction = std::make_shared<DistanceHellinger>();
    } else if (isEqualStr(distName, "kullback")) {
        distanceFunction = std::make_shared<DistanceKullback>();
    } else if (isEqualStr(distName, "cosine")) {
        distanceFunction = std::make_shared<DistanceCosine>();
    } else if (isEqualStr(distName, "mahalanobis")) {
        bool isInvertedCov = false;
        bool isCov = arguments.containsElementNamed("cov");
        arma::Mat<double> cov;
        if (isCov) {
          cov = Rcpp::as<arma::Mat<double>>(arguments["cov"]);
        } else {
          // if data was provided as matrix
          if (this->isDataMatrix) {
            // calc covariance matrix if input data is in matrix format
            cov = arma::cov(dataMatrix);
          } else {
            Rcpp::stop("Calculation of inverted covariance matrix is only supported for input data in matrix format.");
          }
        }
        if (arguments.containsElementNamed("inverted")) {
          isInvertedCov = Rcpp::as<bool >(arguments["inverted"]);
        }
        if (!isInvertedCov) {
          cov = arma::inv(cov);
        }
        distanceFunction = std::make_shared<DistanceMahalanobis>(cov);
    } else if (isEqualStr(distName, "manhattan")) {
        distanceFunction = std::make_shared<DistanceManhattan>();
    } else if (isEqualStr(distName, "maximum")) {
        distanceFunction = std::make_shared<DistanceMaximum>();
    } else if (isEqualStr(distName, "minkowski")) {
        double p = 2;
        if (arguments.containsElementNamed("p")) {
          p = Rcpp::as<double >(arguments["p"]);
        }
        distanceFunction = std::make_shared<DistanceMinkowski>(p);
    } else if (isEqualStr(distName, "podani")) {
        distanceFunction = std::make_shared<DistancePodani>();
    } else if (isEqualStr(distName, "soergel")) {
        distanceFunction = std::make_shared<DistanceSoergel>();
    } else if (isEqualStr(distName, "wave")) {
        distanceFunction = std::make_shared<DistanceWave>();
    } else if (isEqualStr(distName, "whittaker")) {
        distanceFunction = std::make_shared<DistanceWhittaker>();
    } else if (isEqualStr(distName, "binary")) {
        distanceFunction = std::make_shared<DistanceBinary>();
    } else if (isEqualStr(distName, "braun-blanquet")) {
        distanceFunction = std::make_shared<DistanceBraunblanquet>();
    } else if (isEqualStr(distName, "dice")) {
        distanceFunction = std::make_shared<DistanceDice>();
    } else if (isEqualStr(distName, "fager")) {
        distanceFunction = std::make_shared<DistanceFager>();
    } else if (isEqualStr(distName, "faith")) {
        distanceFunction = std::make_shared<DistanceFaith>();
    } else if (isEqualStr(distName, "hamman")) {
        distanceFunction = std::make_shared<DistanceHamman>();
    } else if (isEqualStr(distName, "kulczynski1")) {
        distanceFunction = std::make_shared<DistanceKulczynski1>();
    } else if (isEqualStr(distName, "kulczynski2")) {
        distanceFunction = std::make_shared<DistanceKulczynski2>();
    } else if (isEqualStr(distName, "michael")) {
        distanceFunction = std::make_shared<DistanceMichael>();
    } else if (isEqualStr(distName, "mountford")) {
        distanceFunction = std::make_shared<DistanceMountford>();
    } else if (isEqualStr(distName, "mozley")) {
        distanceFunction = std::make_shared<DistanceMozley>();
    } else if (isEqualStr(distName, "ochiai")) {
        distanceFunction = std::make_shared<DistanceOchiai>();
    } else if (isEqualStr(distName, "phi")) {
        distanceFunction = std::make_shared<DistancePhi>();
    } else if (isEqualStr(distName, "russel")) {
        distanceFunction = std::make_shared<DistanceRussel>();
    } else if (isEqualStr(distName, "simple matching")) {
        distanceFunction = std::make_shared<DistanceSimplematching>();
    } else if (isEqualStr(distName, "simpson")) {
        distanceFunction = std::make_shared<DistanceSimpson>();
    } else if (isEqualStr(distName, "stiles")) {
        distanceFunction = std::make_shared<DistanceStiles>();
    } else if (isEqualStr(distName, "tanimoto")) {
        distanceFunction = std::make_shared<DistanceTanimoto>();
    } else if (isEqualStr(distName, "yule")) {
        distanceFunction = std::make_shared<DistanceYule>();
    } else if (isEqualStr(distName, "yule2")) {
        distanceFunction = std::make_shared<DistanceYule2>();
    } else if (isEqualStr(distName, "hamming")) {
      distanceFunction = std::make_shared<DistanceHamming>();
    } else if (isEqualStr(distName, "custom")) {
        SEXP func_ = arguments["func"];
        funcPtr func = *Rcpp::XPtr<funcPtr>(func_);
        distanceFunction = std::make_shared<DistanceCustom>(func);
    } else {
        distanceFunction = std::make_shared<DistanceEuclidean>();
    }
    return distanceFunction;
}
