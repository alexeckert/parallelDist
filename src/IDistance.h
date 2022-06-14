// IDistance.h
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

#ifndef IDISTANCE_H_
#define IDISTANCE_H_

#include <RcppArmadillo.h>

using arma::mat;
using arma::Mat;
using arma::uword;

typedef double (*funcPtr)(const mat &A, const mat &B);

template <typename T>
Mat<T> colwise_max_idx(const Mat<T> &A) {
    Mat<T> res = mat(1, A.n_cols);
    for (uword i = 0; i != A.n_cols; ++i) {
        res.at(0, i) = A.col(i).max();
    }
    return res;
}

template <typename T>
Mat<T> colwise_min_idx(const Mat<T> &A) {
    Mat<T> res = mat(1, A.n_cols);
    for (uword i = 0; i != A.n_cols; ++i) {
        res.at(0, i) = A.col(i).min();
    }
    return res;
}

class IDistance {
  public:
    virtual ~IDistance() {}
    virtual double calcDistance(const mat &A, const mat &B) = 0;

  protected:
    double res;
    int countFinite;
    int countCol;

    double proportion() { return (double)countFinite/countCol; }
    
    void remove_nan(mat &res) {        
        countFinite = arma::uvec(arma::find_finite( res )).n_elem;
        countCol = res.n_cols;
        res.elem( find_nonfinite(res) ).zeros();
    }

};

#endif // IDISTANCE_H_
