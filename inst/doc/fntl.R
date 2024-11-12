## -----------------------------------------------------------------------------
#| include: false
library(fntl)
library(tidyverse)

set.seed(1234)


## -----------------------------------------------------------------------------
mu_true = 0
x = rnorm(n = 200, mean = mu_true, sd = 1)
loglik = function(mu) { sum(dnorm(x, mu, 1, log = TRUE)) }
optimize(loglik, lower = -100, upper = 100, maximum = TRUE)


## auto f = [](double x, double y) -> double { return x*y; };

## auto f = [](double x, double y) { return x*y; };

## Rcpp::NumericVector x = Rcpp::rnorm(200);
## auto loglik = [&](double mu) {
##     double out = Rcpp::sum(Rcpp::dnorm(x, mu, 1, true));
## 	return out;
## };

## std::function<double(double)> loglik = [&](double mu) {
##     double out = Rcpp::sum(Rcpp::dnorm(x, mu, 1, true));
## 	return out;
## };

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/first.cpp")
out = first_ex(2, 3)
print(out$value)


## // Create args and export to a List
## fntl::integrate_args args0;
## Rcpp::List x = Rcpp::wrap(args0);
## 
## // Instantiate a second args struct from the list x
## x["stop_on_error"] = false;
## fntl::integrate_args args1 = Rcpp::as<fntl::integrate_args>(x);

## // This will cause an exception
## x["abcdefg"] = 0;
## fntl::integrate_args args1 = Rcpp::as<fntl::integrate_args>(x);

## fntl::integrate_result out = fntl::integrate(f, 0, 1);

## fntl::integrate_result out = fntl::integrate(f, 0, 1, args);
## Rcpp::List x = Rcpp::wrap(out);

## fntl::integrate_status status = fntl::integrate_status::OK;  // Define a status
## int err = to_underlying(status);                             // status to int
## status1 = fntl::integrate_status(err);                       // int to status

## typedef function<double(double)> dfd;
## typedef function<double(const NumericVector&)> dfv;
## typedef function<double(const NumericVector&, const NumericVector&)> dfvv;
## typedef function<NumericVector(const NumericVector&)> vfv;
## typedef function<NumericMatrix(const NumericVector&)> mfv;

## double mach_eps = std::numeric_limits<double>::epsilon();
## double mach_eps_2r = sqrt(mach_eps);
## double mach_eps_4r = std::pow(mach_eps, 0.25);
## unsigned int uint_max = std::numeric_limits<unsigned int>::max();

## enum class error_action : unsigned int {
## 	STOP = 3L,     // <1>
## 	WARNING = 2L,  // <2>
## 	MESSAGE = 1L,  // <3>
## 	NONE = 0L      // <4>
## };

## // [[Rcpp::depends(fntl)]]
## #include "fntl.h"
## 
## // [[Rcpp::export]]
## Rcpp::List crash_ex(Rcpp::NumericVector x0)
## {
## 	fntl::dfv f = [](Rcpp::NumericVector x) { return Rcpp::sum(x*x); };
## 	auto out = fntl::gradient(f, x0);
## 	return Rcpp::wrap(out);
## }

## fntl::dfv f = [](Rcpp::NumericVector x) -> double { return Rcpp::sum(x*x); };

## fntl::dfv f = [](Rcpp::NumericVector x) {
## 	double out = Rcpp::sum(x*x);
## 	return out;
## };

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/callr.cpp")
a = 2
b = 3
f = function(x) { x^(a - 1) * (1 - x)^(b - 1) }
out = callr_ex(f)
print(out$value)


## -----------------------------------------------------------------------------
args = integrate_args()
print(args)
a = 2; b = 3
f = function(x) { x^(a-1) * (1-x)^(b-1) }
out = integrate0(f, 0, 1, args)
print(out$value)


## -----------------------------------------------------------------------------
#| eval: false
# source("examples/timing1.R")
# Rcpp::sourceCpp("examples/timing2.cpp")
# Rcpp::sourceCpp("examples/timing3.cpp")
# Rcpp::sourceCpp("examples/timing4.cpp")
# n_levels = c(100, 200, 500, 1000, 10000)


## -----------------------------------------------------------------------------
#| eval: false
# set.seed(1234)
# start = Sys.time()
# timing1_ex(R = 200, n_levels)
# print(Sys.time() - start)


## -----------------------------------------------------------------------------
#| eval: false
# set.seed(1234)
# start = Sys.time()
# timing2_ex(R = 200, n_levels)
# print(Sys.time() - start)


## -----------------------------------------------------------------------------
#| eval: false
# set.seed(1234)
# start = Sys.time()
# timing3_ex(R = 200, n_levels)
# print(Sys.time() - start)


## -----------------------------------------------------------------------------
#| eval: false
# set.seed(1234)
# start = Sys.time()
# timing4_ex(R = 200, n_levels)
# print(Sys.time() - start)


## double sum1p(Rcpp::NumericMatrix x) { return Rcpp::sum(x + 1); }

## double sum1p_1(Rcpp::NumericMatrix& x) { return Rcpp::sum(x + 1); }

## double sum1p_2(const Rcpp::NumericMatrix& x) { return Rcpp::sum(x + 1); }

## integrate_result integrate(
## 	const dfd& f,        // <1>
## 	double lower,                // <2>
## 	double upper,                // <3>
##     const integrate_args& args   // <4>
## )
## 
## integrate_result integrate(
## 	const dfd& f,        // <1>
## 	double lower,                // <2>
## 	double upper                 // <3>
## )

## struct integrate_args {
## 	unsigned int subdivisions = 100L;  // <1>
## 	double rel_tol = mach_eps_4r;      // <2>
## 	double abs_tol = mach_eps_4r;      // <3>
## 	bool stop_on_error = true;         // <4>
## };

## struct integrate_result {
## 	double value;             // <1>
## 	double abs_error;         // <2>
## 	int subdivisions;         // <3>
## 	integrate_status status;  // <4>
## 	int n_eval;               // <5>
## 	std::string message;      // <6>
## 
## 	operator SEXP() const;    // <7>
## };

## enum class integrate_status : int {
## 	OK = 0L,                                  // <1>
## 	MAX_SUBDIVISIONS = 1L,                    // <2>
## 	ROUNDOFF_ERROR = 2L,                      // <3>
## 	BAD_INTEGRAND_BEHAVIOR = 3L,              // <4>
## 	ROUNDOFF_ERROR_EXTRAPOLATION_TABLE = 4L,  // <5>
## 	PROBABLY_DIVERGENT_INTEGRAL = 5L,         // <6>
## 	INVALID_INPUT = 6L                        // <7>
## };

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/integrate.cpp")
out = integrate_ex(0 , Inf)
print(out)


## double fd_deriv(
## 	const dfv& f,                                  // <1>
## 	const Rcpp::NumericVector& x,                  // <2>
## 	unsigned int i,                                // <3>
## 	double h,                                      // <5>
## 	const fd_types& fd_type = fd_types::SYMMETRIC  // <7>
## )
## 
## double fd_deriv2(
## 	const dfv& f,                                  // <1>
## 	const Rcpp::NumericVector& x,                  // <2>
## 	unsigned int i,                                // <3>
## 	unsigned int j,                                // <4>
## 	double h_i,                                    // <5>
## 	double h_j,                                    // <6>
## 	const fd_types& fd_type = fd_types::SYMMETRIC  // <7>
## )

## enum class fd_types : unsigned int {
## 	SYMMETRIC = 0L,
## 	FORWARD = 1L,
## 	BACKWARD = 2L
## };

## -----------------------------------------------------------------------------
n = 3
b = seq_len(n)
x0 = rep(0.5, n)
eps = 0.001  ## Fix a step size for finite differences
g = function(x) { b * cos(sum(b * x)) }
h = function(x) { -tcrossprod(b) * sin(sum(b * x)) }


## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/fd-deriv.cpp")
out3 = out2 = out1 = numeric(n)
for (i in 1:n) {
	out1[i] = fd_deriv_ex(x0, i-1, eps, type = 0)  ## Symmetric
	out2[i] = fd_deriv_ex(x0, i-1, eps, type = 1)  ## Forward
	out3[i] = fd_deriv_ex(x0, i-1, eps, type = 2)  ## Backward
}
print(out1)
print(out2)
print(out3)
print(g(x0))


## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/fd-deriv2.cpp")
out3 = out2 = out1 = matrix(NA, n, n)
for (i in 1:n) {
for (j in 1:n) {
	out1[i,j] = fd_deriv2_ex(x0, i-1, j-1, eps, eps, type = 0)  ## Symmetric
	out2[i,j] = fd_deriv2_ex(x0, i-1, j-1, eps, eps, type = 1)  ## Forward
	out3[i,j] = fd_deriv2_ex(x0, i-1, j-1, eps, eps, type = 2)  ## Backward
}
}
print(out1)
print(out2)
print(out3)
h(x0)


## richardson_result deriv(
## 	const dfv& f,                                  // <1>
## 	const Rcpp::NumericVector& x,                  // <2>
## 	unsigned int i,                                // <3>
## 	const richardson_args& args,                   // <5>
## 	const fd_types& fd_type = fd_types::SYMMETRIC  // <6>
## )
## 
## richardson_result deriv(
## 	const dfv& f,                                  // <1>
## 	const Rcpp::NumericVector& x,                  // <2>
## 	unsigned int i,                                // <3>
## 	const fd_types& fd_type = fd_types::SYMMETRIC  // <6>
## )
## 
## richardson_result deriv2(
## 	const dfv& f,                                  // <1>
## 	const Rcpp::NumericVector& x,                  // <2>
## 	unsigned int i,                                // <3>
## 	unsigned int j,                                // <4>
## 	const richardson_args& args,                   // <5>
## 	const fd_types& fd_type = fd_types::SYMMETRIC  // <6>
## )
## 
## richardson_result deriv2(
## 	const dfv& f,                                  // <1>
## 	const Rcpp::NumericVector& x,                  // <2>
## 	unsigned int i,                                // <3>
## 	unsigned int j,                                // <4>
## 	const fd_types& fd_type = fd_types::SYMMETRIC  // <6>
## )

## struct richardson_args
## {
## 	double delta = 0.5;                 // <1>
## 	unsigned int maxiter = 10;          // <2>
## 	double h = 1;                       // <3>
## 	double tol = mach_eps_4r;           // <4>
## 	double accuracy_factor = R_PosInf;  // <5>
## 
## 	richardson_args() { };              // <6>
## 	richardson_args(SEXP obj);          // <7>
## 	operator SEXP() const;              // <8>
## };
## 

## struct richardson_result
## {
## 	double value;              // <1>
## 	double err;                // <2>
## 	unsigned int iter;         // <3>
## 	richardson_status status;  // <4>
## 
## 	operator SEXP() const;     // <5>
## };

## enum class richardson_status : unsigned int {
## 	OK = 0L,                  // <1>
## 	NOT_CONVERGED = 1L,       // <2>
## 	NUMERICAL_PRECISION = 2L  // <3>
## };

## -----------------------------------------------------------------------------
n = 3
x0 = rep(0.5, n)


## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/deriv.cpp")
out3 = out2 = out1 = numeric(n)
for (i in 1:n) {
	out1[i] = deriv_ex(x0, i-1, type = 0)$value  ## Symmetric
	out2[i] = deriv_ex(x0, i-1, type = 1)$value  ## Forward
	out3[i] = deriv_ex(x0, i-1, type = 2)$value  ## Backward
}
print(out1)
print(out2)
print(out3)


## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/deriv2.cpp")
out3 = out2 = out1 = matrix(NA, n, n)
for (i in 1:n) {
for (j in 1:n) {
	out1[i,j] = deriv2_ex(x0, i-1, j-1, type = 0)$value  ## Symmetric
	out2[i,j] = deriv2_ex(x0, i-1, j-1, type = 1)$value  ## Forward
	out3[i,j] = deriv2_ex(x0, i-1, j-1, type = 2)$value  ## Backward
}
}
print(out1)
print(out2)
print(out3)


## gradient_result gradient(
## 	const dfv& f,                                  // <1>
## 	const Rcpp::NumericVector& x,                  // <2>
## 	const richardson_args& args,                   // <3>
## 	const fd_types& fd_type = fd_types::SYMMETRIC  // <4>
## )
## 
## gradient_result gradient(
## 	const dfv& f,                                  // <1>
## 	const Rcpp::NumericVector& x,                  // <2>
## 	const fd_types& fd_type = fd_types::SYMMETRIC  // <4>
## )

## struct gradient_result {
## 	std::vector<double> value;       // <1>
## 	std::vector<double> err;         // <2>
## 	std::vector<unsigned int> iter;  // <3>
## 
## 	operator SEXP() const;           // <4>
## };

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/gradient.cpp")
gradient_ex(x0 = 1:4)


## -----------------------------------------------------------------------------
f = function(x) { sum(x^2) }
numDeriv::grad(f, x = 1:4)


## jacobian_result jacobian(
## 	const vfv& f,                                 // <1>
## 	const Rcpp::NumericVector& x,                 // <2>
## 	const richardson_args& args                   // <3>
## 	const fd_types& fd_type = fd_types::SYMMETRIC // <4>
## )
## 
## jacobian_result jacobian(
## 	const vfv& f,                                 // <1>
## 	const Rcpp::NumericVector& x,                 // <2>
## 	const fd_types& fd_type = fd_types::SYMMETRIC // <4>
## )

## struct jacobian_result
## {
## 	std::vector<double> value;       // <1>
## 	std::vector<double> err;         // <2>
## 	std::vector<unsigned int> iter;  // <3>
## 	double rows;                     // <4>
## 	double cols;                     // <5>
## 
## 	operator SEXP() const;           // <6>
## };

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/jacobian.cpp")
out = jacobian_ex(x0 = 1:4)
names(out)
print(out$value)


## -----------------------------------------------------------------------------
f = function(x) { cumsum(sin(x)) }
numDeriv::jacobian(f, x = 1:4)


## hessian_result hessian(
## 	const dfv& f,                                  // <1>
## 	const Rcpp::NumericVector& x,                  // <2>
## 	const richardson_args& args,                   // <3>
## 	const fd_types& fd_type = fd_types::SYMMETRIC  // <4>
## )
## 
## hessian_result hessian(
## 	const dfv& f,                                  // <1>
## 	const Rcpp::NumericVector& x,                  // <2>
## 	const fd_types& fd_type = fd_types::SYMMETRIC  // <4>
## )

## struct hessian_result
## {
## 	std::vector<double> value;       // <1>
## 	std::vector<double> err;         // <2>
## 	std::vector<unsigned int> iter;  // <3>
## 	double dim;                      // <4>
## 
## 	operator SEXP() const;           // <5>
## };

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/hessian.cpp")
out = hessian_ex(x0 = c(1,2))
print(out)


## -----------------------------------------------------------------------------
f = function(x) { sum(sin(x)) }
numDeriv::hessian(f, x = c(1,2))


## struct findroot_args {
## 	double tol = mach_eps_4r;                      // <1>
## 	unsigned int maxiter = 1000;                   // <2>
## 	error_action action = error_action::STOP;      // <3>
## 	
## 	findroot_args() { };                           // <4>
## 	findroot_args(SEXP obj);                       // <5>
## 	operator SEXP() const;                         // <6>
## };

## struct findroot_result {
## 	double root;                 // <1>
## 	double f_root;               // <2>
## 	unsigned int iter;           // <3>
## 	double tol;                  // <4>
## 	findroot_status status;      // <5>
## 	std::string message;         // <6>
## 
## 	operator SEXP() const;       // <7>
## };

## enum class findroot_status : unsigned int {
## 	OK = 0L,                  // <1>
## 	NUMERICAL_OVERFLOW = 1L,  // <2>
## 	NOT_CONVERGED = 2L        // <3>
## };

## findroot_result findroot_bisect(
## 	const dfd& f,              // <1>
## 	double lower,              // <2>
## 	double upper,              // <3>
## 	const findroot_args& args  // <4>
## )
## 
## findroot_result findroot_bisect(
## 	const dfd& f,              // <1>
## 	double lower,              // <2>
## 	double upper               // <3>
## )

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/findroot-bisect.cpp")
out = findroot_bisect_ex(0, 10)
print(out)


## findroot_result findroot_brent(
## 	const dfd& f,              // <1>
## 	double lower,              // <2>
## 	double upper,              // <3>
## 	const findroot_args& args  // <4>
## )
## 
## findroot_result findroot_brent(
## 	const dfd& f,              // <1>
## 	double lower,              // <2>
## 	double upper               // <3>
## )

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/findroot-brent.cpp")
out = findroot_brent_ex(0, 10)
print(out)


## struct optimize_args
## {
## 	bool fnscale = 1;                          // <1>
## 	double tol = mach_eps_2r;                  // <2>
## 	unsigned int maxiter = 1000;               // <3>
## 	unsigned int report_period = uint_max;     // <4>
## 	error_action action = error_action::STOP;  // <5>
## 
## 	optimize_args() { };                       // <6>
## 	optimize_args(SEXP obj);                   // <7>
## 	operator SEXP() const;                     // <8>
## };
## 
## 

## struct optimize_result
## {
## 	double par;                   // <1>
## 	double value;                 // <2>
## 	unsigned int iter;            // <3>
## 	double tol;                   // <4>
## 	optimize_status status;       // <5>
## 	std::string message;          // <6>
## 
## 	operator SEXP() const;        // <7>
## };

## enum class optimize_status : unsigned int {
## 	OK = 0L,                  // <1>
## 	NUMERICAL_OVERFLOW = 1L,  // <2>
## 	NOT_CONVERGED = 2L        // <3>
## };

## optimize_result goldensection(
## 	const dfd& f,              // <1>
## 	double lower,              // <2>
## 	double upper,              // <3>
## 	const optimize_args& args  // <4>
## )
## 
## optimize_result goldensection(
## 	const dfd& f,              // <1>
## 	double lower,              // <2>
## 	double upper               // <3>
## )

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/goldensection.cpp")
out = goldensection_ex(0, 1)
print(out)


## optimize_result optimize_brent(
## 	const dfd& f,              // <1>
## 	double lower,              // <2>
## 	double upper,              // <3>
## 	const optimize_args& args  // <4>
## )
## 
## optimize_result optimize_brent(
## 	const dfd& f,              // <1>
## 	double lower,              // <2>
## 	double upper               // <3>
## )

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/optimize-brent.cpp")
out = optimize_brent_ex(0, 1)
print(out)


## neldermead_result neldermead(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const neldermead_args& args       // <3>
## )
## 
## neldermead_result neldermead(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f                      // <2>
## )

## struct neldermead_args
## {
## 	double alpha = 1.0;           // <1>
## 	double beta = 0.5;            // <2>
## 	double gamma = 2.0;           // <3>
## 	unsigned int trace = 0;       // <4>
## 	double abstol = R_NegInf;     // <5>
## 	double reltol = mach_eps_2r;  // <6>
## 	unsigned int maxit = 500;     // <7>
## 	double fnscale = 1.0;         // <8>
## 
## 	neldermead_args() { };        // <9>
## 	neldermead_args(SEXP obj);    // <10>
## 	operator SEXP() const;        // <11>
## };

## struct neldermead_result {
## 	std::vector<double> par;   // <1>
## 	double value;              // <2>
## 	neldermead_status status;  // <3>
## 	int fncount;               // <4>
## 
## 	operator SEXP() const;     // <5>
## };

## enum class neldermead_status : unsigned int {
## 	OK = 0L,                 // <1>
## 	NOT_CONVERGED = 1L,      // <2>
## 	SIMLEX_DEGENERACY = 10L  // <3>
## };

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/neldermead.cpp")
out = neldermead_ex(x0 = c(1, -1))
print(out)


## bfgs_result bfgs(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const vfv& g,                     // <3>
## 	const bfgs_args& args             // <4>
## )
## 
## bfgs_result bfgs(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const bfgs_args& args             // <4>
## )
## 
## bfgs_result bfgs(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const vfv& g                      // <3>
## )
## 
## bfgs_result bfgs(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f                      // <2>
## )

## struct bfgs_args {
## 	richardson_args deriv_args;    // <1>
## 	double parscale = 1;           // <2>
## 	int trace = 0;                 // <3>
## 	double fnscale = 1;            // <4>
## 	int maxit = 100;               // <5>
## 	int report = 10;               // <6>
## 	double abstol = R_NegInf;      // <7>
## 	double reltol = mach_eps_2r;   // <8>
## 
## 	bfgs_args() { };               // <9>
## 	bfgs_args(SEXP obj);           // <10>
## 	operator SEXP() const;         // <11>
## };

## struct bfgs_result {
## 	std::vector<double> par;  // <1>
## 	double value;             // <2>
## 	bfgs_status status;       // <3>
## 	int fncount;              // <4>
## 	int grcount;              // <5>
## 
## 	operator SEXP() const;    // <6>
## };

## enum class bfgs_status : unsigned int {
## 	OK = 0L,             // <1>
## 	NOT_CONVERGED = 1L   // <2>
## };

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/bfgs.cpp")
out = bfgs_ex(x0 = rep(1, 4))
print(out$numerical)
print(out$analytical)


## lbfgsb_result lbfgsb(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const vfv& g,                     // <3>
## 	const lbfgsb_args& args           // <4>
## )
## 
## lbfgsb_result lbfgsb(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const lbfgsb_args& args           // <4>
## )
## 
## lbfgsb_result lbfgsb(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const vfv& g                      // <3>
## )
## 
## lbfgsb_result lbfgsb(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f                      // <2>
## )

## struct lbfgsb_args {
## 	std::vector<double> lower;     // <1>
## 	std::vector<double> upper;     // <2>
## 	richardson_args deriv_args;    // <3>
## 	double parscale = 1;           // <4>
## 	int trace = 0;                 // <5>
## 	double fnscale = 1;            // <6>
## 	int lmm = 5;                   // <7>
## 	int maxit = 100;               // <8>
## 	int report = 10;               // <9>
## 	double factr = 1e7;            // <10>
## 	double pgtol = 0;              // <11>
## 
## 	lbfgsb_args() { };             // <12>
## 	lbfgsb_args(SEXP obj);         // <13>
## 	operator SEXP() const;         // <14>
## };

## struct lbfgsb_result {
## 	std::vector<double> par;  // <1>
## 	double value;             // <2>
## 	lbfgsb_status status;     // <3>
## 	int fncount;              // <4>
## 	int grcount;              // <5>
## 	std::string msg;          // <6>
## 
## 	operator SEXP() const;    // <7>
## };

## enum class lbfgsb_status : unsigned int {
## 	OK = 0L,             // <1>
## 	NOT_CONVERGED = 1L,  // <2>
## 	WARN = 51L,          // <3>
## 	ERROR = 52L,         // <4>
## };

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/lbfgsb.cpp")
out = lbfgsb_ex(x0 = rep(1, 4))
print(out$numerical)
print(out$analytical)


## cg_result bfgs(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const vfv& g,                     // <3>
## 	const cg_args& args               // <4>
## )
## 
## cg_result cg(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const cg_args& args               // <4>
## )
## 
## cg_result cg(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const vfv& g                      // <3>
## )
## 
## cg_result cg(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f                      // <2>
## )

## 
## struct cg_args
## {
## 	richardson_args deriv_args;    // <1>
## 	double parscale = 1;           // <2>
## 	double fnscale = 1;            // <3>
## 	double abstol = R_NegInf;      // <4>
## 	double reltol = mach_eps_2r;   // <5>
## 	int type = 1;                  // <6>
## 	int trace = 0;                 // <7>
## 	int maxit = 100;               // <8>
## 
## 	cg_args() { };                 // <9>
## 	cg_args(SEXP obj);             // <10>
## 	operator SEXP() const;         // <11>
## };

## struct cg_result {
## 	std::vector<double> par;  // <1>
## 	double value;             // <2>
## 	cg_status status;         // <3>
## 	int fncount;              // <4>
## 	int grcount;              // <5>
## 
## 	operator SEXP() const;    // <6>
## };

## enum class cg_status : unsigned int {
## 	OK = 0L,             // <1>
## 	NOT_CONVERGED = 1L   // <2>
## };

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/cg.cpp")
out = cg_ex(x0 = rep(1, 4))
print(out$numerical)
print(out$analytical)


## 
## nlm_result nlm(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const vfv& g,                     // <3>
## 	const mfv& h,                     // <4>
## 	const nlm_args& args              // <5>
## )
## 
## nlm_result nlm(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const vfv& g,                     // <3>
## 	const nlm_args& args              // <5>
## )
## 
## nlm_result nlm(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const nlm_args& args              // <5>
## )
## 
## nlm_result nlm(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const vfv& g,                     // <3>
## 	const mfv& h                      // <4>
## )
## 
## nlm_result nlm(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f,                     // <2>
## 	const vfv& g                      // <3>
## )
## 
## nlm_result nlm(
## 	const Rcpp::NumericVector& init,  // <1>
## 	const dfv& f                      // <2>
## )

## struct nlm_args
## {
## 	std::vector<double> typsize;   // <1>
## 	unsigned int print_level = 0;  // <2>
## 	double fscale = 1;             // <3>
## 	double fnscale = 1;            // <4>
## 	unsigned int ndigit = 12;      // <5>
## 	double gradtol = 1e-6;         // <6>
## 	double stepmax = R_PosInf;     // <7>
## 	double steptol = 1e-6;         // <8>
## 	int iterlim = 100;             // <9>
## 	unsigned int method = 1;       // <10>
## 	double trust_radius = 1.0;     // <11>
## 
## 	nlm_args() { };                // <12>
## 	nlm_args(SEXP obj);            // <13>
## 	operator SEXP() const;         // <14>
## };

## struct nlm_result
## {
## 	std::vector<double> par;      // <1>
## 	std::vector<double> grad;     // <2>
## 	double estimate;              // <3>
## 	int iterations;               // <4>
## 	nlm_status status;            // <5>
## 
## 	operator SEXP() const;        // <6>
## };

## enum class nlm_status : unsigned int {
## 	OK = 0L,                   // <1>
## 	GRADIENT_WITHIN_TOL = 1L,  // <2>
## 	ITERATES_WITH_TOL = 2L,    // <3>
## 	NO_LOWER_STEP = 3L,        // <4>
## 	ITERATION_MAX = 4L,        // <5>
## 	STEP_SIZE_EXCEEDED = 5L    // <6>
## };
## 

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/nlm.cpp")
out = nlm_ex(x0 = rep(1, 4))
nn = c("par", "grad")
print(out$res1[nn])
print(out$res2[nn])
print(out$res3[nn])


## -----------------------------------------------------------------------------
#| eval: false
# apply(X, c(1,2), f)
# apply(X, 1, f)
# apply(X, 2, f)


## template <typename T, int RTYPE>
## Rcpp::Vector<RTYPE> row_apply(
## 	const Rcpp::Matrix<RTYPE>& X,							// <1>
## 	const std::function<T(const Rcpp::Vector<RTYPE>&)>& f   // <2>
## )
## 
## template <typename T, int RTYPE>
## Rcpp::Vector<RTYPE> col_apply(
## 	const Rcpp::Matrix<RTYPE>& X,							// <1>
## 	const std::function<T(const Rcpp::Vector<RTYPE>&)>& f   // <2>
## )
## 
## template <typename T, int RTYPE>
## Rcpp::Matrix<RTYPE> mat_apply(
## 	const Rcpp::Matrix<RTYPE>& X,							// <1>
## 	const std::function<T(T)>& f                            // <2>
## )

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/apply.cpp")
X = matrix(1:12, 4, 3)
out = apply_ex(X)
print(out)


## Rcpp::NumericMatrix outer(
## 	const Rcpp::NumericMatrix& X,  // <1>
## 	const dfvv& f                  // <3>
## )
## 
## Rcpp::NumericMatrix outer(
## 	const Rcpp::NumericMatrix& X,  // <1>
## 	const Rcpp::NumericMatrix& Y,  // <2>
## 	const dfvv& f                  // <3>
## )
## 
## Rcpp::NumericVector outer_matvec(
## 	const Rcpp::NumericMatrix& X,  // <1>
## 	const dfvv& f,                 // <3>
## 	const Rcpp::NumericVector& a   // <4>
## )
## 
## Rcpp::NumericVector outer_matvec(
## 	const Rcpp::NumericMatrix& X,  // <1>
## 	const Rcpp::NumericMatrix& Y,  // <2>
## 	const dfvv& f,                 // <3>
## 	const Rcpp::NumericVector& a   // <4>
## )

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/outer.cpp")
m = 5; n = 3; d = 2
X = matrix(rnorm(10), m, d)
Y = matrix(rnorm(6), n, d)
a = rep(1, m)
b = rep(1, n)
out = outer_ex(X, Y, a, b)
print(out)


## cg_result solve_cg(
## 	const vfv& l,                     // <1>
## 	const Rcpp::NumericVector& b,     // <2>
## 	const Rcpp::NumericVector& init,  // <3>
## 	const cg_args& args               // <4>
## )
## 
## cg_result solve_cg(
## 	const vfv& l,                     // <1>
## 	const Rcpp::NumericVector& b,     // <2>
## 	const Rcpp::NumericVector& init   // <3>
## )
## 
## cg_result solve_cg(
## 	const vfv& l,                     // <1>
## 	const Rcpp::NumericVector& b      // <2>
## )

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/solve-cg.cpp")
b = rep(1, 10)
out = solve_cg_ex(b)
print(out)


## -----------------------------------------------------------------------------
A = matrix(0, 10, 10)
diag(A) = 2
A[cbind(1:9, 2:10)] = 1
A[cbind(2:10, 1:9)] = 1
solve(A, b)


## -----------------------------------------------------------------------------
#| eval: false
# which(f(X), arr.ind = TRUE)


## template <typename T, int RTYPE>
## Rcpp::IntegerMatrix which(
## 	const Rcpp::Matrix<RTYPE>& X,      // <1>
## 	const std::function<bool(T)>& f),  // <2>
## 	bool one_based = false             // <3>
## )

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/which.cpp")
x = runif(10, -1, 1)
X = matrix(x, 2, 5)
out = which_ex(X)
print(X)
print(out)


## -----------------------------------------------------------------------------
f = function(x) { x > 0 & x < 0.5 }
which(f(X), arr.ind = TRUE) - 1


## typedef std::function<double(double,bool)> density;       // <1>
## typedef std::function<double(double,bool,bool)> cdf;      // <2>
## typedef std::function<double(double,bool,bool)> quantile; // <3>

## double d_trunc(
## 	double x,                      // <1>
## 	double lo,                     // <2>
## 	double hi,                     // <3>
## 	const density& f,              // <4>
## 	const cdf& F,                  // <5>
## 	bool log = false               // <6>
## )
## 
## Rcpp::NumericVector d_trunc(
## 	const Rcpp::NumericVector& x,  // <1>
## 	const Rcpp::NumericVector& lo, // <2>
## 	const Rcpp::NumericVector& hi, // <3>
## 	const density& f,              // <4>
## 	const cdf& F,                  // <5>
## 	bool log = false               // <6>
## )

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/d_trunc.cpp")
x = seq(0, 1, length.out = 30)
d_trunc_ex(x, shape1 = 2, shape2 = 5, lo = 0.5, hi = 0.95)


## double p_trunc(
## 	double x,              // <1>
## 	double lo,             // <2>
## 	double hi,             // <3>
## 	const cdf& F,          // <4>
## 	bool lower = true,     // <5>
## 	bool log = false       // <6>
## )
## 
## Rcpp::NumericVector p_trunc(
## 	const Rcpp::NumericVector& x,  // <1>
## 	const Rcpp::NumericVector& lo, // <2>
## 	const Rcpp::NumericVector& hi, // <3>
## 	const cdf& F,                  // <4>
## 	bool lower = true,             // <5>
## 	bool log = false               // <6>
## )

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/p_trunc.cpp")
x = seq(0, 1, length.out = 30)
p_trunc_ex(x, shape1 = 2, shape2 = 5, lo = 0.5, hi = 0.95)


## double q_trunc(
## 	double p,              // <1>
## 	double lo,             // <2>
## 	double hi,             // <3>
## 	const cdf& F,          // <4>
## 	const quantile& Finv,  // <5>
## 	bool lower = true,     // <6>
## 	bool log = false       // <7>
## )
## 
## Rcpp::NumericVector q_trunc(
## 	const Rcpp::NumericVector& p,  // <1>
## 	const Rcpp::NumericVector& lo, // <2>
## 	const Rcpp::NumericVector& hi, // <3>
## 	const cdf& F,                  // <4>
## 	const quantile& Finv,          // <5>
## 	bool lower = true,             // <6>
## 	bool log = false               // <7>
## )

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/q_trunc.cpp")
p = seq(0, 1, length.out = 30)
q_trunc_ex(p, shape1 = 2, shape2 = 5, lo = 0.5, hi = 0.95)


## double r_trunc(
## 	double lo,             // <2>
## 	double hi,             // <3>
## 	const cdf& F,          // <4>
## 	const quantile& Finv   // <5>
## )
## 
## Rcpp::NumericVector r_trunc(
## 	unsigned int n,                // <1>
## 	const Rcpp::NumericVector& lo, // <2>
## 	const Rcpp::NumericVector& hi, // <3>
## 	const cdf& F,                  // <4>
## 	const quantile& Finv           // <5>
## )
## 

## -----------------------------------------------------------------------------
Rcpp::sourceCpp("examples/r_trunc.cpp")
r_trunc_ex(n = 20, shape1 = 2, shape2 = 5, lo = 0.5, hi = 0.95)

