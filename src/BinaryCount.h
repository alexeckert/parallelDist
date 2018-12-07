// BinaryCount.h
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

#ifndef BINARYCOUNT_H_
#define BINARYCOUNT_H_

#include <RcppArmadillo.h>

class BinaryCount {
 private:
    uint64_t a;
    uint64_t b;
    uint64_t c;
    uint64_t d;

 public:
    BinaryCount(uint64_t a, uint64_t b, uint64_t c, uint64_t d) : a(a), b(b), c(c), d(d) {}
    ~BinaryCount() {}
    static BinaryCount getBinaryCount(const arma::mat &A, const arma::mat &B) {
        uint64_t a = 0;
        uint64_t b = 0;
        uint64_t c = 0;
        uint64_t d = 0;

        for (arma::uword idx = 0; idx < A.size(); ++idx) {
            bool aZero = A.at(idx) == 0.0;
            bool bZero = B.at(idx) == 0.0;

            if (!aZero && !bZero) {
                ++a;
            } else if (!aZero && bZero) {
                ++b;
            } else if (aZero && !bZero) {
                ++c;
            } else if (aZero && bZero) {
                ++d;
            }
        }

        return BinaryCount(a, b, c, d);
    }
    uint64_t getA() {
        return a;
    }
    uint64_t getB() {
        return b;
    }
    uint64_t getC() {
        return c;
    }
    uint64_t getD() {
        return d;
    }
};

#endif  // BINARYCOUNT_H_
