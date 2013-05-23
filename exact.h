// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#ifndef volmodels_exact_h
#define volmodels_exact_h

#include <complex>
#include <vector>
#include "utils.h"

using namespace std;

static double EPS = 10E-12;
static int STDEVS = 12;
static int DIGITS = 12;

typedef std::complex<double> Complex;

inline Complex gamma(const double& kappa, const double& sigma, const double& a);

Complex firstMomentIntegralSquareRootProcess(double kappa, double theta, double sigma, double u, double vu, double t, double vt);

Complex secondMomentIntegralSquareRootProcess(double kappa, double theta, double sigma, double u, double vu, double t, double vt);

double quantileSquareRootProcess(double p, double kappa, double theta, double sigma, double u, double vu, double t);

double quantileIntSquareRootProcess(double p, double kappa, double theta, double sigma, double u, double vu, double t, double vt);

double quantileHestonProcess(double p1, double p2, double p3, double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t);

#endif
