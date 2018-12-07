// DistanceDTWFactory.cpp
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

#include "DistanceDTWFactory.h"
#include "StepPattern.h"
#include "Util.h"

std::shared_ptr<IDistance> DistanceDTWFactory::createDistanceFunction(
  const std::string& distName, const Rcpp::List& arguments) {
    using util::isEqualStr;

    std::shared_ptr<IDistance> distanceFunction = NULL;
    unsigned int windowSize = 0;
    NormMethod normMethod = NormMethod::NoNorm;
    bool warpingWindow = false;
    std::string stepPatternName = "symmetric1";

    warpingWindow = arguments.containsElementNamed("window.size");
    if (warpingWindow) {
        windowSize = Rcpp::as<unsigned int >(arguments["window.size"]);
    }
    if (arguments.containsElementNamed("norm.method")) {
        std::string normMethodStr = Rcpp::as<std::string >(arguments["norm.method"]);
        if (isEqualStr(normMethodStr, "n")) {
            normMethod = NormMethod::ALength;
        } else if (isEqualStr(normMethodStr, "n+m")) {
            normMethod = NormMethod::ABLength;
        } else if (isEqualStr(normMethodStr, "path.length")) {
            normMethod = NormMethod::PathLength;
        }
    }
    if (arguments.containsElementNamed("step.pattern")) {
        stepPatternName = Rcpp::as<std::string >(arguments["step.pattern"]);
    }

    if (isEqualStr(stepPatternName, "asymmetric")) {
        distanceFunction = std::make_shared<StepPatternAsymmetric>(warpingWindow, windowSize, normMethod);
    } else if (isEqualStr(stepPatternName, "asymmetricP0")) {
        distanceFunction = std::make_shared<StepPatternAsymmetricP0>(warpingWindow, windowSize, normMethod);
    } else if (isEqualStr(stepPatternName, "asymmetricP05")) {
        distanceFunction = std::make_shared<StepPatternAsymmetricP05>(warpingWindow, windowSize, normMethod);
    } else if (isEqualStr(stepPatternName, "asymmetricP1")) {
        distanceFunction = std::make_shared<StepPatternAsymmetricP1>(warpingWindow, windowSize, normMethod);
    } else if (isEqualStr(stepPatternName, "asymmetricP2")) {
        distanceFunction = std::make_shared<StepPatternAsymmetricP2>(warpingWindow, windowSize, normMethod);
    } else if (
      isEqualStr(stepPatternName, "symmetric2") ||
      isEqualStr(stepPatternName, "symmetricP0")
    ) {
        distanceFunction = std::make_shared<StepPatternSymmetric2>(warpingWindow, windowSize, normMethod);
    } else if (isEqualStr(stepPatternName, "symmetricP05")) {
        distanceFunction = std::make_shared<StepPatternSymmetricP05>(warpingWindow, windowSize, normMethod);
    } else if (isEqualStr(stepPatternName, "symmetricP1")) {
        distanceFunction = std::make_shared<StepPatternSymmetricP1>(warpingWindow, windowSize, normMethod);
    } else if (isEqualStr(stepPatternName, "symmetricP2")) {
        distanceFunction = std::make_shared<StepPatternSymmetricP2>(warpingWindow, windowSize, normMethod);
    } else {
        distanceFunction = std::make_shared<StepPatternSymmetric1>(warpingWindow, windowSize, normMethod);
    }
    return distanceFunction;
}
