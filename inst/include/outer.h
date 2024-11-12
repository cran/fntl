#ifndef FNTL_OUTER_H
#define FNTL_OUTER_H

#include <Rcpp.h>
#include "typedefs.h"
#include "typedefs-rcpp.h"
#include "util.h"

namespace fntl {

inline Rcpp::NumericMatrix outer(const Rcpp::NumericMatrix& X, const dfvv& f)
{
	unsigned int n = X.nrow();
	Rcpp::NumericMatrix out(n, n);

	for (unsigned int j = 0; j < n; j++) {
		for (unsigned int i = 0; i < j; i++) {
			out(i,j) = f(X.row(i), X.row(j));
			out(j,i) = out(i,j);
		}
	}

	for (unsigned int i = 0; i < n; i++) {
		out(i,i) = f(X.row(i), X.row(i));
	}

	return out;
}

inline Rcpp::NumericMatrix outer(const Rcpp::NumericMatrix& X,
	const Rcpp::NumericMatrix& Y, const dfvv& f)
{
	unsigned int m = X.nrow();
	unsigned int n = Y.nrow();
	Rcpp::NumericMatrix out(m, n);

	for (unsigned int j = 0; j < n; j++) {
		for (unsigned int i = 0; i < m; i++) {
			out(i,j) = f(X.row(i), Y.row(j));
		}
	}

	return out;
}

inline Rcpp::NumericVector outer_matvec(const Rcpp::NumericMatrix& X,
	const dfvv& f, const Rcpp::NumericVector& a)
{
	unsigned int n = X.nrow();
	if (n != a.size()) {
		Rcpp::stop("Dimension mismatch");
	}

	Rcpp::NumericVector out(n, 0);

	for (unsigned int j = 0; j < n; j++) {
		for (unsigned int i = 0; i < j; i++) {
			double f_ij = f(X.row(i), X.row(j));
			double f_ji = f_ij;
			out(i) += f_ij * a(j);
			out(j) += f_ji * a(i);
		}
	}

	for (unsigned int i = 0; i < n; i++) {
		double f_ii = f(X.row(i), X.row(i));
		out(i) += f_ii * a(i);
	}

	return out;
}

inline Rcpp::NumericVector outer_matvec(const Rcpp::NumericMatrix& X,
	const Rcpp::NumericMatrix& Y, const dfvv& f,
	const Rcpp::NumericVector& a)
{
	unsigned int m = X.nrow();
	unsigned int n = Y.nrow();

	if (n != a.size()) {
		Rcpp::stop("Dimension mismatch");
	}

	Rcpp::NumericVector out(m, 0);

	for (unsigned int j = 0; j < n; j++) {
		for (unsigned int i = 0; i < m; i++) {
			double f_ij = f(X.row(i), Y.row(j));
			out(i) += f_ij * a(j);
		}
	}

	return out;
}

}

#endif
