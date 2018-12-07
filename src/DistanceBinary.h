// DistanceBinary.h
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

#ifndef DISTANCEBINARY_H_
#define DISTANCEBINARY_H_

#include <RcppArmadillo.h>
#include <math.h>
#include "IDistance.h"
#include "Util.h"
#include "BinaryCount.h"

#undef max
#define minOfPair(x, y) ((x) < (y) ? (x) : (y))
#define maxOfPair(x, y) ((x) < (y) ? (y) : (x))

//=======================
// Binary distance
//=======================
class DistanceBinary : public IDistance {
public:
    double calcDistance(const arma::mat &A, const arma::mat &B) {
        BinaryCount bc = BinaryCount::getBinaryCount(A, B);
        uint64_t denominator = bc.getA() + bc.getB() + bc.getC();
        return ((denominator == 0) ? 0 : static_cast<double>(bc.getB() + bc.getC()) / denominator);
    }
};

//=======================
// Braun-Blanquet
//=======================
class DistanceBraunblanquet : public IDistance {
public:
  double calcDistance(const arma::mat &A, const arma::mat &B) {
    BinaryCount bc = BinaryCount::getBinaryCount(A, B);
    uint64_t denominator = maxOfPair((bc.getA() + bc.getB()), (bc.getA() + bc.getC()));
    return util::similarityToDistance(static_cast<double>(bc.getA()) / denominator);
  }
};

//=======================
// Dice distance
//=======================
class DistanceDice : public IDistance {
public:
  double calcDistance(const arma::mat &A, const arma::mat &B) {
    BinaryCount bc = BinaryCount::getBinaryCount(A, B);
    uint64_t denominator = 2 * bc.getA() + bc.getB() + bc.getC();
    return util::similarityToDistance(static_cast<double>(2 * bc.getA()) / denominator);
  }
};

//=======================
// Fager distance (like in proxy)
//=======================
class DistanceFager : public IDistance {
public:
  double calcDistance(const arma::mat &A, const arma::mat &B) {
    BinaryCount bc = BinaryCount::getBinaryCount(A, B);
    return util::similarityToDistance((static_cast<double>(bc.getA()) /
    sqrt((bc.getA() + bc.getB()) * (bc.getA() + bc.getC()))) - ((sqrt(bc.getA() + bc.getC())) / 2.0));
  }
};

//=======================
// Faith distance
//=======================
class DistanceFaith : public IDistance {
public:
  double calcDistance(const arma::mat &A, const arma::mat &B) {
    BinaryCount bc = BinaryCount::getBinaryCount(A, B);
    return util::similarityToDistance((bc.getA() + static_cast<double>(bc.getD()) / 2.0) / A.n_cols);
  }
};

//=======================
// Hamman distance
//=======================
class DistanceHamman : public IDistance {
public:
  double calcDistance(const arma::mat &A, const arma::mat &B) {
    BinaryCount bc = BinaryCount::getBinaryCount(A, B);
    return util::similarityToDistance(
      (static_cast<double>(bc.getA()) + bc.getD() - bc.getB() - bc.getC()) / A.n_cols);
  }
};

//=======================
// Kulczynski1 distance
//=======================
class DistanceKulczynski1 : public IDistance {
public:
    double calcDistance(const arma::mat &A, const arma::mat &B) {
        BinaryCount bc = BinaryCount::getBinaryCount(A, B);
        return util::similarityToDistance(static_cast<double>(bc.getA()) / (bc.getB() + bc.getC()));
    }
};

//=======================
// Kulczynski2 distance
//=======================
class DistanceKulczynski2 : public IDistance {
public:
    double calcDistance(const arma::mat &A, const arma::mat &B) {
        BinaryCount bc = BinaryCount::getBinaryCount(A, B);
        double div1 = static_cast<double>(bc.getA()) / (bc.getA() + bc.getB());
        double div2 = static_cast<double>(bc.getA()) / (bc.getA() + bc.getC());
        return util::similarityToDistance((div1 + div2) / 2.0);
    }
};

//=======================
// Michael distance
//=======================
class DistanceMichael : public IDistance {
public:
  double calcDistance(const arma::mat &A, const arma::mat &B) {
    BinaryCount bc = BinaryCount::getBinaryCount(A, B);
    double denominator = pow(bc.getA() + bc.getD(), 2) + pow(bc.getB() + bc.getC(), 2);
    return util::similarityToDistance((4.0 * (static_cast<double>(bc.getA() * bc.getD()) -
    (bc.getB() * bc.getC()))) / denominator);
  }
};

//=======================
// Mountford distance
//=======================
class DistanceMountford : public IDistance {
public:
    double calcDistance(const arma::mat &A, const arma::mat &B) {
        BinaryCount bc = BinaryCount::getBinaryCount(A, B);
        uint64_t denominator = bc.getA() * (bc.getB() + bc.getC()) + 2 * bc.getB() * bc.getC();
        return util::similarityToDistance(static_cast<double>(2 * bc.getA()) / denominator);
    }
};

//=======================
// Mozley distance
//=======================
class DistanceMozley : public IDistance {
public:
  double calcDistance(const arma::mat &A, const arma::mat &B) {
    BinaryCount bc = BinaryCount::getBinaryCount(A, B);
    uint64_t denominator = (bc.getA() + bc.getB()) * (bc.getA() + bc.getC());
    return util::similarityToDistance((static_cast<double>(bc.getA() * A.n_cols)) / denominator);
  }
};

//=======================
// Ochiai distance
//=======================
class DistanceOchiai : public IDistance {
public:
  double calcDistance(const arma::mat &A, const arma::mat &B) {
    BinaryCount bc = BinaryCount::getBinaryCount(A, B);
    double denominator = sqrt((bc.getA() + bc.getB()) * (bc.getA() + bc.getC()));
    return util::similarityToDistance(static_cast<double>(bc.getA()) / denominator);
  }
};

//=======================
// Phi distance
//=======================
class DistancePhi : public IDistance {
public:
  double calcDistance(const arma::mat &A, const arma::mat &B) {
    BinaryCount bc = BinaryCount::getBinaryCount(A, B);
    double denominator = (sqrt(bc.getA() + bc.getB()) * sqrt(bc.getC() + bc.getD()) *
      sqrt(bc.getA() + bc.getC()) * sqrt(bc.getB() + bc.getD()));
    return util::similarityToDistance(
      (static_cast<double>(bc.getA() * bc.getD()) - (bc.getB() * bc.getC()))/ denominator);
  }
};

//=======================
// Russel distance
//=======================
class DistanceRussel : public IDistance {
public:
    double calcDistance(const arma::mat &A, const arma::mat &B) {
        BinaryCount bc = BinaryCount::getBinaryCount(A, B);
        return util::similarityToDistance(static_cast<double>(bc.getA()) / A.n_cols);
    }
};

//=======================
// SimpleMatching distance
//=======================
class DistanceSimplematching : public IDistance {
public:
    double calcDistance(const arma::mat &A, const arma::mat &B) {
        BinaryCount bc = BinaryCount::getBinaryCount(A, B);
        return util::similarityToDistance(static_cast<double>(bc.getA() + bc.getD()) / A.n_cols);
    }
};

//=======================
// Simpson distance
//=======================
class DistanceSimpson : public IDistance {
public:
  double calcDistance(const arma::mat &A, const arma::mat &B) {
    BinaryCount bc = BinaryCount::getBinaryCount(A, B);
    uint64_t denominator = minOfPair((bc.getA() + bc.getB()), (bc.getA() + bc.getC()));
    return util::similarityToDistance(static_cast<double>(bc.getA()) / denominator);
  }
};

//=======================
// Stiles distance
//=======================
class DistanceStiles : public IDistance {
public:
  double calcDistance(const arma::mat &A, const arma::mat &B) {
    BinaryCount bc = BinaryCount::getBinaryCount(A, B);
    unsigned int n = A.n_cols;
    return util::similarityToDistance(
      (log(n) + 2 * log(fabs(static_cast<double>(bc.getA() * bc.getD()) - bc.getB() * bc.getC()) - n / 2.0) -
      log(bc.getA() + bc.getB()) - log(bc.getC() + bc.getD()) -
      log(bc.getA() + bc.getC()) - log(bc.getB() + bc.getD())));
  }
};

//=======================
// Tanimoto distance
//=======================
class DistanceTanimoto : public IDistance {
public:
    double calcDistance(const arma::mat &A, const arma::mat &B) {
        BinaryCount bc = BinaryCount::getBinaryCount(A, B);
        uint64_t denominator = bc.getA() + 2 * bc.getB() + 2 * bc.getC() + bc.getD();
        return util::similarityToDistance(static_cast<double>(bc.getA() + bc.getD()) / denominator);
    }
};

//=======================
// Yule distance
//=======================
class DistanceYule : public IDistance {
public:
    double calcDistance(const arma::mat &A, const arma::mat &B) {
        BinaryCount bc = BinaryCount::getBinaryCount(A, B);
        uint64_t denominator = (bc.getA() * bc.getD()) + (bc.getB() * bc.getC());
        return util::similarityToDistance((static_cast<double>(bc.getA() * bc.getD()) -
        (bc.getB() * bc.getC())) / denominator);
    }
};

//=======================
// Yule2 distance
//=======================
class DistanceYule2 : public IDistance {
public:
    double calcDistance(const arma::mat &A, const arma::mat &B) {
        BinaryCount bc = BinaryCount::getBinaryCount(A, B);
        double denominator = sqrt(bc.getA() * bc.getD()) + sqrt(bc.getB() * bc.getC());
        return util::similarityToDistance((static_cast<double>(sqrt(bc.getA() * bc.getD())) -
        sqrt(bc.getB() * bc.getC())) / denominator);
    }
};

#endif  // DISTANCEBINARY_H_
