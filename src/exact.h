#ifndef _volmodels_exact
#define _volmodels_exact

#include <Rcpp.h>

RcppExport SEXP pintsqr_c(SEXP kappa_, SEXP theta_, SEXP sigma_, SEXP u_, SEXP vu_, SEXP t_, SEXP vt_, SEXP x_);

#endif
