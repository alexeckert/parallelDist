// parallelDist_init.c
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

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _parallelDist_cpp_parallelDistMatrixVec(SEXP, SEXP, SEXP);
extern SEXP _parallelDist_cpp_parallelDistVec(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_parallelDist_cpp_parallelDistMatrixVec", (DL_FUNC) &_parallelDist_cpp_parallelDistMatrixVec, 3},
  {"_parallelDist_cpp_parallelDistVec",       (DL_FUNC) &_parallelDist_cpp_parallelDistVec,       3},
  {NULL, NULL, 0}
};

void R_init_parallelDist(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
