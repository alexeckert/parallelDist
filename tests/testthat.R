## testthat.R
##
## Copyright (C)  2017  Alexander Eckert
##
## This file is part of parallelDist.
##
## parallelDist is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## parallelDist is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with parallelDist. If not, see <http://www.gnu.org/licenses/>.

library(testthat)
library(parallelDist)

isCran <- function() {
  !identical(Sys.getenv("NOT_CRAN"), "true")
}

# workaround for error "cannot open file 'startup.Rs': No such file or directory"
Sys.setenv("R_TESTS" = "")
test_check("parallelDist")
