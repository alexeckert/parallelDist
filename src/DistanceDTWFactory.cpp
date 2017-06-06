// DistanceDTWFactory.cpp
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

#include "DistanceDTWFactory.h"
#include "StepPattern.h"
#include "Utility.h"

std::shared_ptr<IDistance> DistanceDTWFactory::createDistanceFunction(std::string& distName, Rcpp::List& arguments) {
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
        if (utility::isEqualStr(normMethodStr, "n")) {
            normMethod = NormMethod::ALength;
        } else if (utility::isEqualStr(normMethodStr, "n+m")) {
            normMethod = NormMethod::ABLength;
        } else if (utility::isEqualStr(normMethodStr, "path.length")) {
            normMethod = NormMethod::PathLength;
        }
    }
    if (arguments.containsElementNamed("step.pattern")) {
        stepPatternName = Rcpp::as<std::string >(arguments["step.pattern"]);
    }

    if (utility::isEqualStr(stepPatternName, "asymmetric")) {
        distanceFunction = std::make_shared<DistanceDTWGeneric<StepPatternAsymmetric>>(warpingWindow, windowSize, normMethod);
    } else if (utility::isEqualStr(stepPatternName, "asymmetricP0")) {
        distanceFunction = std::make_shared<DistanceDTWGeneric<StepPatternAsymmetricP0>>(warpingWindow, windowSize, normMethod);
    } else if (utility::isEqualStr(stepPatternName, "asymmetricP05")) {
        distanceFunction = std::make_shared<DistanceDTWGeneric<StepPatternAsymmetricP05>>(warpingWindow, windowSize, normMethod);
    } else if (utility::isEqualStr(stepPatternName, "asymmetricP1")) {
        distanceFunction = std::make_shared<DistanceDTWGeneric<StepPatternAsymmetricP1>>(warpingWindow, windowSize, normMethod);
    } else if (utility::isEqualStr(stepPatternName, "asymmetricP2")) {
        distanceFunction = std::make_shared<DistanceDTWGeneric<StepPatternAsymmetricP2>>(warpingWindow, windowSize, normMethod);
    } else if (utility::isEqualStr(stepPatternName, "symmetric2") || utility::isEqualStr(stepPatternName, "symmetricP0")) {
        distanceFunction = std::make_shared<DistanceDTWGeneric<StepPatternSymmetric2>>(warpingWindow, windowSize, normMethod);
    } else if (utility::isEqualStr(stepPatternName, "symmetricP05")) {
        distanceFunction = std::make_shared<DistanceDTWGeneric<StepPatternSymmetricP05>>(warpingWindow, windowSize, normMethod);
    } else if (utility::isEqualStr(stepPatternName, "symmetricP1")) {
        distanceFunction = std::make_shared<DistanceDTWGeneric<StepPatternSymmetricP1>>(warpingWindow, windowSize, normMethod);
    } else if (utility::isEqualStr(stepPatternName, "symmetricP2")) {
        distanceFunction = std::make_shared<DistanceDTWGeneric<StepPatternSymmetricP2>>(warpingWindow, windowSize, normMethod);
    } else {
        distanceFunction = std::make_shared<DistanceDTWGeneric<StepPatternSymmetric1>>(warpingWindow, windowSize, normMethod);
    }
    return distanceFunction;
}
