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
    unsigned long a;
    unsigned long b;
    unsigned long c;
    unsigned long d;
public:
    BinaryCount (unsigned long a, unsigned long b, unsigned long c, unsigned long d) : a(a), b(b), c(c), d(d) {};
    ~BinaryCount () {};
    static BinaryCount getBinaryCount(const arma::mat &A, const arma::mat &B) {
        unsigned long a = 0;
        unsigned long b = 0;
        unsigned long c = 0;
        unsigned long d = 0;

        for(arma::uword idx=0; idx < A.size(); ++idx) {
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
    unsigned long getA() {
        return a;
    };
    unsigned long getB() {
        return b;
    };
    unsigned long getC() {
        return c;
    };
    unsigned long getD() {
        return d;
    };
};

#endif
