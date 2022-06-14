// Util.h
//
// Copyright (C)  2017, 2021  Alexander Eckert
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

#ifndef UTIL_H_
#define UTIL_H_

#include <cmath>
#include <memory>
#include <string>
#include <RcppArmadillo.h>

namespace util {

bool isEqualStr(const std::string &str1, std::string str2);

double similarityToDistance(const double distance);

inline double proportion(const int countFinite, const int countCol) { return (double)countFinite/countCol; }

void remove_nan(arma::mat &res, int &countFinite, int &countCol);

} // namespace util

#endif // UTIL_H_
