// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include "SVJtrajectory.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
#include <climits>
#include <boost/random.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include "exact.h"
#include "bessel.h"
#include "utils.h"

//#include <ctime>

#ifndef PI
    #define PI M_PI
#endif

using namespace std;

typedef std::complex<double> Complex;

double quantileHestonProcessSVJ(double p1, double p2, double p3, double p4, int N, double lambda, double mu_s, double sigma_s, double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t) {

	//calculating mu_bar
	double mu_bar=exp(mu_s + 0.5*sigma_s*sigma_s) - 1.0;
	
	//calculating S_bar
	double vt = quantileSquareRootProcess(p1, kappa, theta, sigma, u, vu, t);
    double intv = quantileIntSquareRootProcess(p2, kappa, theta, sigma, u, vu, t, vt);
    double intsqrvdW = (vt - vu - kappa*theta*(t-u) + kappa * intv) / sigma;
    double m = ((r-lambda*mu_bar)*(t-u) - 0.5*intv + rho*intsqrvdW);//m would be different
    double s = sqrt((1-rho*rho)*intv);
    double S_bar = Su * exp(m + s*qnorm(p3));

	//got Y
	double ln_Y_mean = N * mu_s;
	double ln_Y_var  = N * sigma_s*sigma_s;
	double Y      = exp( ln_Y_mean + sqrt(ln_Y_var) * qnorm(p4)) ;
	
	return S_bar * Y;
}
