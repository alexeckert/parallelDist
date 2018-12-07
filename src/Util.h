// Util.h
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

#ifndef UTIL_H_
#define UTIL_H_

#include <memory>
#include <cmath>
#include <string>

namespace util {

bool isEqualStr(const std::string& str1, std::string str2);

double similarityToDistance(const double distance);

}

#endif  // UTIL_H_
