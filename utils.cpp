// 
// Dale Roberts <dale.o.roberts@gmail.com>
//

#define MATHLIB_STANDALONE
#include <Rmath.h>

//#undef qchisq
//#undef qpois
//#undef qgamma
//#undef qbinom
#undef qnorm
//#undef pnorm

#include "utils.h"
#include "bessel.h"

#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
//#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/special_functions/gamma.hpp>

using boost::math::tgamma;

double qchisq(double p, double df, double ncp) 
{
    boost::math::non_central_chi_squared dist(df, ncp);
    return quantile(dist, p);
}

double qnorm(double p) 
{
    boost::math::normal dist;
    return quantile(dist, p);
}

//double qpois(double p, double lambda) 
//{
//    boost::math::poisson dist(lambda);
//    return quantile(dist, p);
//}

double qgamma(double p, double shape, double scale)
{
    boost::math::gamma_distribution<double> dist(shape, scale);
    return quantile(dist, p);   
}

//double qbinom(double x, double n, double p)
//{
//    boost::math::binomial dist(n, p);
//    return quantile(dist, x);
//}

double pnorm(double z) 
{
    boost::math::normal dist;
    return cdf(dist, z);
}

//double pnorm(double q) {
//    return Rf_pnorm5(q, 0.0, 1.0, 1, 0);
//}
//
//double qnorm(double q) {
//    return Rf_qnorm5(q, 0.0, 1.0, 1, 0);
//}
//
//double qchisq(double p, double df, double ncp) {
//    return Rf_qnchisq(p, df, ncp, 1, 0);
//}
//
double qpois(double p, double lambda) {
    return qpois(p, lambda, 1, 0);
}
//
//double qgamma(double p, double shape, double scale) {
//    return Rf_qgamma(p, shape, scale, 1, 0);
//}

double qbinom(double p, double size, double prob) {
    return qbinom(p, size, prob, 1, 0);
}

int qbessel(double p, double nu, double z) 
{
    // see Glasserman-Kim, p.16
    
    int B = 0, count = 0;
    double probmass = pow(z/2, nu) / (besseli(z,nu).real()*gammafn(nu+1));
    double cprob = probmass;
    
    while (p > cprob & count < 1001) 
    {
        count++;
        probmass = z*z*probmass/(4*count*(count+nu));
        cprob = cprob + probmass;
    }
    
    B = count;
    
    return B;
}

double bscall(double S, double K, double t1, double t2, double r, double sigma) 
{
    double tau = t2-t1;
    double d1 = (log(S/K)+(r+0.5*(sigma*sigma))*tau)/(sigma*sqrt(tau));
    double d2 = d1-sigma*sqrt(tau);
    return S*pnorm(d1)-K*(exp(-r*tau))*pnorm(d2);
}

bool isPowerOfTwo(int x) 
{
    return ((x != 0) && ((x & (~x + 1)) == x));
}

int log2(int val) 
{
    int ret = -1;
    while (val != 0) {
        val >>= 1;
        ret++;
    }
    return ret;
}
