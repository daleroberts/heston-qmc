// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
#include <climits>
#include <boost/random.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

#include <boost/math/tools/toms748_solve.hpp>

#include "exact.h"
#include "bessel.h"
#include "utils.h"

#ifndef PI
    #define PI M_PI
#endif

using namespace std;

typedef std::complex<double> Complex;

struct CdfIntSquareRootProcessBridgeWithDeriv
{
    CdfIntSquareRootProcessBridgeWithDeriv(double kappa_, double theta_, double sigma_, double u_, double vu_, double t_, double vt_, double p_) 
        : kappa(kappa_), theta(theta_), sigma(sigma_), u(u_), vu(vu_), t(t_), vt(vt_), p(p_) {};
    
    boost::math::tuple<double, double> operator()(double const& x) {
        double eps = EPS;
        int m = STDEVS; // # of standard deviations to use
        
        int k = 1;
        int add2pi = 0;
     
        double delta = 4*theta*kappa/(sigma*sigma);
        Complex m1 = firstMomentIntegralSquareRootProcess(kappa, theta, sigma, u, vu, t, vt);
        Complex m2 = secondMomentIntegralSquareRootProcess(kappa, theta, sigma, u, vu, t, vt);

        double ueps = abs(m1) + m*abs(std::sqrt(m2-m1*m1));
        double h = 2*PI / (x+ueps);

        Complex b1 = besseli(4*kappa*sqrt(vu*vt)*exp(-0.5*kappa*(t-u)) / ((sigma*sigma)*(1-exp(-kappa*(t-u)))), 0.5*delta-1);

        Complex ga = gamma(kappa, sigma, h*k);
        Complex zc = sqrt(vu*vt)*4*ga*exp(-0.5*ga*(t-u))/((sigma*sigma)*(Complex(1,0)-exp(-ga*(t-u))));

        double angle = 0.0;
        double oldangle = atan2(zc.imag(), zc.real());

        Complex b2 = exp(add2pi*(0.5*delta-1)*2*PI*Complex(0,1))*besseli(zc, 0.5*delta-1);
        Complex t1 = ga*exp(-0.5*(ga-kappa)*(t-u))*(1-exp(-kappa*(t-u)))/(kappa*(Complex(1,0)-exp(-ga*(t-u))));
        Complex t2 = exp((vu+vt)/(sigma*sigma)*(kappa*(Complex(1,0)+exp(-kappa*(t-u)))/(1-exp(-kappa*(t-u))) - ga*(Complex(1,0)+exp(-ga*(t-u)))/(Complex(1,0)-exp(-ga*(t-u)))));
        
        double cc = (h*x/PI)+(2/PI)*sin(h*k*x)*(b2*t1*t2).real()/(b1.real()*k);
        double dd = (h/PI)+(2/PI)*h*cos(h*k*x)*(b2*t1*t2).real()/(b1.real()); // deriv wrt x        
        double tv = (b2*t1*t2).real()/b1.real();

        while (std::abs(tv)/k >= PI*eps/2) {
            
            if (k > 5000) {
                #ifndef _OPENMP
                cout << "k reached max: " << abs(tv)/k << " " << PI*eps/2 <<  endl;
                #endif
                break;
            }

            k = k + 1;
            
            ga = gamma(kappa, sigma, h*k);
            zc = sqrt(vu*vt)*4.0*ga*exp(-0.5*ga*(t-u))/((sigma*sigma)*(Complex(1,0)-exp(-ga*(t-u))));

            angle = atan2(zc.imag(), zc.real());
            if ((zc.real() < 0) && (zc.imag() < 0) && (oldangle * angle < 0))
                add2pi += 1;
            oldangle = angle;

            b2 = exp(add2pi*(0.5*delta-1)*2*PI*Complex(0,1))*besseli(zc, 0.5*delta-1);
            t1 = ga*exp(-0.5*(ga-kappa)*(t-u))*(1-exp(-kappa*(t-u)))/(kappa*(Complex(1,0)-exp(-ga*(t-u))));
            t2 = exp((vu+vt)/(sigma*sigma)*(kappa*(Complex(1,0)+exp(-kappa*(t-u)))/(1-exp(-kappa*(t-u))) - ga*(Complex(1,0)+exp(-ga*(t-u)))/(Complex(1,0)-exp(-ga*(t-u)))));

            tv = (b2*t1*t2).real()/b1.real();
            
            cc += (2/PI)*sin(h*k*x)*tv/k;
            dd += (2/PI)*h*cos(h*k*x)*tv;
        }
                
      return boost::math::make_tuple(cc-p, dd);
    }

private:
   double kappa, theta, sigma, u, vu, t, vt, p;
};

struct CdfIntSquareRootProcessBridge
{
    CdfIntSquareRootProcessBridge(double kappa_, double theta_, double sigma_, double u_, double vu_, double t_, double vt_, double p_) 
    : kappa(kappa_), theta(theta_), sigma(sigma_), u(u_), vu(vu_), t(t_), vt(vt_), p(p_) {};
    
    double operator()(double const& x) {
        double eps = EPS;
        int m = STDEVS; // # of standard deviations to use
        
        int k = 1;
        int add2pi = 0;
        
        double delta = 4*theta*kappa/(sigma*sigma);
        Complex m1 = firstMomentIntegralSquareRootProcess(kappa, theta, sigma, u, vu, t, vt);
        Complex m2 = secondMomentIntegralSquareRootProcess(kappa, theta, sigma, u, vu, t, vt);
        
        double ueps = abs(m1) + m*abs(std::sqrt(m2-m1*m1));
        double h = 2*PI / (x+ueps);
        
        Complex b1 = besseli(4*kappa*sqrt(vu*vt)*exp(-0.5*kappa*(t-u)) / ((sigma*sigma)*(1-exp(-kappa*(t-u)))), 0.5*delta-1);
        
        Complex ga = gamma(kappa, sigma, h*k);
        Complex zc = sqrt(vu*vt)*4*ga*exp(-0.5*ga*(t-u))/((sigma*sigma)*(Complex(1,0)-exp(-ga*(t-u))));
        
        double angle = 0.0;
        double oldangle = atan2(zc.imag(), zc.real());
        
        Complex b2 = exp(add2pi*(0.5*delta-1)*2*PI*Complex(0,1))*besseli(zc, 0.5*delta-1);
        Complex t1 = ga*exp(-0.5*(ga-kappa)*(t-u))*(1-exp(-kappa*(t-u)))/(kappa*(Complex(1,0)-exp(-ga*(t-u))));
        Complex t2 = exp((vu+vt)/(sigma*sigma)*(kappa*(Complex(1,0)+exp(-kappa*(t-u)))/(1-exp(-kappa*(t-u))) - ga*(Complex(1,0)+exp(-ga*(t-u)))/(Complex(1,0)-exp(-ga*(t-u)))));
        
        double cc = (h*x/PI)+(2/PI)*sin(h*k*x)*(b2*t1*t2).real()/(b1.real()*k);
        double dd = (h/PI)+(2/PI)*h*cos(h*k*x)*(b2*t1*t2).real()/(b1.real()); // deriv wrt x        
        double tv = (b2*t1*t2).real()/b1.real();
        
        while (std::abs(tv)/k >= PI*eps/2) {
            
            if (k > 5000) {
                #ifndef _OPENMP
                cout << "k reached max: " << abs(tv)/k << " " << PI*eps/2 <<  endl;
                #endif
                break;
            }
            
            k = k + 1;
            
            ga = gamma(kappa, sigma, h*k);
            zc = sqrt(vu*vt)*4.0*ga*exp(-0.5*ga*(t-u))/((sigma*sigma)*(Complex(1,0)-exp(-ga*(t-u))));
            
            angle = atan2(zc.imag(), zc.real());
            if ((zc.real() < 0) && (zc.imag() < 0) && (oldangle * angle < 0))
                add2pi += 1;
            oldangle = angle;
            
            b2 = exp(add2pi*(0.5*delta-1)*2*PI*Complex(0,1))*besseli(zc, 0.5*delta-1);
            t1 = ga*exp(-0.5*(ga-kappa)*(t-u))*(1-exp(-kappa*(t-u)))/(kappa*(Complex(1,0)-exp(-ga*(t-u))));
            t2 = exp((vu+vt)/(sigma*sigma)*(kappa*(Complex(1,0)+exp(-kappa*(t-u)))/(1-exp(-kappa*(t-u))) - ga*(Complex(1,0)+exp(-ga*(t-u)))/(Complex(1,0)-exp(-ga*(t-u)))));
            
            tv = (b2*t1*t2).real()/b1.real();
            
            cc += (2/PI)*sin(h*k*x)*tv/k;
            dd += (2/PI)*h*cos(h*k*x)*tv;
        }
        
        return cc-p;
    }
    
private:
    double kappa, theta, sigma, u, vu, t, vt, p;
};

inline Complex gamma(const double& kappa, const double& sigma, const double& a) {
    return sqrt(kappa*kappa-2*sigma*sigma*a*Complex(0,1));
}
 
Complex firstMomentIntegralSquareRootProcess(double kappa, double theta, double sigma, double u, double vu, double t, double vt) {
    using namespace std;

    double d = 4*theta*kappa/(sigma*sigma);

    Complex z = Complex(0,-1)*pow(kappa,-2)*pow(sigma,2) + Complex(0,0.5)*(t - u)*pow(kappa,-1)*pow(sigma,2) + Complex(0,1)*(t - u)*exp(kappa*(-t + u))*pow(kappa,-1)*pow(sigma,2)*pow(1 - exp(kappa*(-t + u)),-1) + (vt + vu)*pow(sigma,-2)*(Complex(0,-1)*(t - u)*exp(kappa*(-t + u))*(1 + exp(kappa*(-t + u))) * pow(sigma,2) * pow(1 - exp(kappa*(-t + u)),-2) - Complex(0,1)*(t - u) * exp(kappa*(-t + u))*pow(sigma,2)* pow(1 - exp(kappa*(-t + u)),-1) + Complex(0,1)*(1 + exp(kappa*(-t + u)))*pow(kappa,-1)* pow(sigma,2) * pow(1 - exp(kappa*(-t + u)),-1));

    Complex t1 = Complex(0,4)*(t - u)*exp((-3*kappa*(t - u))/2.)* pow(1 - exp(kappa*(-t + u)),-2) * sqrt(vt)*sqrt(vu) + Complex(0,2)*(t - u)*exp(-(kappa*(t - u))/2.)* pow(1 - exp(kappa*(-t + u)),-1) * sqrt(vt)*sqrt(vu) - Complex(0,4)*exp(-(kappa*(t - u))/2.)*pow(kappa,-1) * pow(1 - exp(kappa*(-t + u)),-1)*sqrt(vt)*sqrt(vu);
 
    double tt = 4*kappa*exp(-(kappa*(t - u))/2.)*pow(sigma,-2)* pow(1 - exp(kappa*(-t + u)),-1) * sqrt(vt)*sqrt(vu);
    
    Complex t2 = besseli(4*kappa*exp(-(kappa*(t - u))/2.)*pow(sigma,-2)* pow(1 - exp(kappa*(-t + u)),-1) * sqrt(vt)*sqrt(vu),-2 + d/2.) + besseli(4*kappa*exp(-(kappa*(t - u))/2.)*pow(sigma,-2)* pow(1 - exp(kappa*(-t + u)),-1) * sqrt(vt)*sqrt(vu), d/2.);
    
    Complex t3 = Complex(2,0)*besseli(4*kappa*exp(-(kappa*(t - u))/2.)*pow(sigma,-2)*pow(1 - exp(kappa*(-t + u)),-1) * sqrt(vt)*sqrt(vu),-1 + d/2.);

    z += t1*t2/t3;

    return z;
}

Complex secondMomentIntegralSquareRootProcess(double kappa, double theta, double sigma, double u, double vu, double t, double vt) {
    Complex z = (pow(sigma,4)*(t - u)*besseli((4*
          exp(-((t - u)*kappa)/2.)*sqrt(vt*vu)*
          kappa)/
        ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-(kappa*(t - u))))*
      exp(-((t - u)*kappa) - 
        ((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma)))/
    (pow(kappa,3)*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*
          sqrt(vt*vu))/((sigma*sigma)*(1 - exp(-(kappa*(t - u)))))
        ,-1 + (2*kappa*theta)/(sigma*sigma))*
      pow(1 - exp(-((t - u)*kappa)),2)) + 
   (pow(sigma,4)*besseli((4*exp(-((t - u)*kappa)/2.)*
          sqrt(vt*vu)*kappa)/
        ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-(kappa*(t - u))))*
      exp(-((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma)))/
    (kappa*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*sqrt(vt*vu))/
        ((sigma*sigma)*(1 - exp(-(kappa*(t - u))))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-((t - u)*kappa)))*
      pow(kappa*kappa,1.5)) + 
   (Complex(0,1)*(t - u)*besseli((4*
          exp(-((t - u)*kappa)/2.)*sqrt(vt*vu)*
          kappa)/
        ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-(kappa*(t - u))))*
      exp(-((t - u)*kappa) - 
        ((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma))*(sigma*sigma)*
      ((Complex(0,0.5)*(t - u)*(sigma*sigma))/kappa + 
        ((vt + vu)*((Complex(0,-1)*(t - u)*
                exp(-((t - u)*kappa))*(sigma*sigma))/
              (1 - exp(-((t - u)*kappa))) - 
             (Complex(0,1)*(t - u)*
                exp(-((t - u)*kappa))*
                (1 + exp(-((t - u)*kappa)))*
                (sigma*sigma))/
              pow(1 - exp(-((t - u)*kappa)),2) + 
             (Complex(0,1)*(1 + exp(-((t - u)*kappa)))*
                (sigma*sigma))/
              ((1 - exp(-((t - u)*kappa)))*
                kappa)))/(sigma*sigma)))/
    (kappa*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*sqrt(vt*vu))/
        ((sigma*sigma)*(1 - exp(-(kappa*(t - u))))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      pow(1 - exp(-((t - u)*kappa)),2)) - 
   (Complex(0,2)*besseli((4*exp(-((t - u)*kappa)/2.)*
          sqrt(vt*vu)*kappa)/
        ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-(kappa*(t - u))))*
      exp(-((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma))*(sigma*sigma)*
      ((Complex(0,0.5)*(t - u)*(sigma*sigma))/kappa + 
        ((vt + vu)*((Complex(0,-1)*(t - u)*
                exp(-((t - u)*kappa))*(sigma*sigma))/
              (1 - exp(-((t - u)*kappa))) - 
             (Complex(0,1)*(t - u)*
                exp(-((t - u)*kappa))*
                (1 + exp(-((t - u)*kappa)))*
                (sigma*sigma))/
              pow(1 - exp(-((t - u)*kappa)),2) + 
             (Complex(0,1)*(1 + exp(-((t - u)*kappa)))*
                (sigma*sigma))/
              ((1 - exp(-((t - u)*kappa)))*
                kappa)))/(sigma*sigma)))/
    (kappa*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*sqrt(vt*vu))/
        ((sigma*sigma)*(1 - exp(-(kappa*(t - u))))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-((t - u)*kappa)))*kappa) + 
   (Complex(0,1)*(t - u)*besseli((4*
          exp(-((t - u)*kappa)/2.)*sqrt(vt*vu)*
          kappa)/
        ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-(kappa*(t - u))))*
      exp(-((t - u)*kappa) - 
        ((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma))*(sigma*sigma)*
      ((Complex(0,1.5)*(t - u)*(sigma*sigma))/kappa + 
        ((vt + vu)*((Complex(0,-1)*(t - u)*
                exp(-((t - u)*kappa))*(sigma*sigma))/
              (1 - exp(-((t - u)*kappa))) - 
             (Complex(0,1)*(t - u)*
                exp(-((t - u)*kappa))*
                (1 + exp(-((t - u)*kappa)))*
                (sigma*sigma))/
              pow(1 - exp(-((t - u)*kappa)),2) + 
             (Complex(0,1)*(1 + exp(-((t - u)*kappa)))*
                (sigma*sigma))/
              ((1 - exp(-((t - u)*kappa)))*
                kappa)))/(sigma*sigma)))/
    (kappa*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*sqrt(vt*vu))/
        ((sigma*sigma)*(1 - exp(-(kappa*(t - u))))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      pow(1 - exp(-((t - u)*kappa)),2)) - 
   (2*pow(sigma,4)*besseli((4*
          exp(-((t - u)*kappa)/2.)*sqrt(vt*vu)*
          kappa)/
        ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-(kappa*(t - u))))*
      exp(-2*(t - u)*kappa - 
        ((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma))*((t - u)*(t - u)))/
    (kappa*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*sqrt(vt*vu))/
        ((sigma*sigma)*(1 - exp(-(kappa*(t - u))))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      pow(1 - exp(-((t - u)*kappa)),3)*
      kappa) + 
   (Complex(0,1)*(t - u)*(besseli((4*
            exp(-((t - u)*kappa)/2.)*sqrt(vt*vu)*
            kappa)/
          ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
         (2*kappa*theta)/(sigma*sigma)) + 
        besseli((4*exp(-((t - u)*kappa)/2.)*
            sqrt(vt*vu)*kappa)/
          ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
         -2 + (2*kappa*theta)/(sigma*sigma)))*
      (1 - exp(-(kappa*(t - u))))*
      exp(-((t - u)*kappa) - 
        ((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma))*(sigma*sigma)*
      ((Complex(0,4)*(t - u)*
           exp((-3*(t - u)*kappa)/2.)*sqrt(vt*vu))/
         pow(1 - exp(-((t - u)*kappa)),2) + 
        (Complex(0,2)*(t - u)*exp(-((t - u)*kappa)/2.)*
           sqrt(vt*vu))/(1 - exp(-((t - u)*kappa))) - 
        (Complex(0,4)*exp(-((t - u)*kappa)/2.)*
           sqrt(vt*vu))/
         ((1 - exp(-((t - u)*kappa)))*
           kappa)))/
    (kappa*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*sqrt(vt*vu))/
        ((sigma*sigma)*(1 - exp(-(kappa*(t - u))))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      pow(1 - exp(-((t - u)*kappa)),2)) - 
   (Complex(0,1)*(besseli((4*exp(-((t - u)*kappa)/2.)*
            sqrt(vt*vu)*kappa)/
          ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
         (2*kappa*theta)/(sigma*sigma)) + 
        besseli((4*exp(-((t - u)*kappa)/2.)*
            sqrt(vt*vu)*kappa)/
          ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
         -2 + (2*kappa*theta)/(sigma*sigma)))*
      (1 - exp(-(kappa*(t - u))))*
      exp(-((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma))*(sigma*sigma)*
      ((Complex(0,4)*(t - u)*
           exp((-3*(t - u)*kappa)/2.)*sqrt(vt*vu))/
         pow(1 - exp(-((t - u)*kappa)),2) + 
        (Complex(0,2)*(t - u)*exp(-((t - u)*kappa)/2.)*
           sqrt(vt*vu))/(1 - exp(-((t - u)*kappa))) - 
        (Complex(0,4)*exp(-((t - u)*kappa)/2.)*
           sqrt(vt*vu))/
         ((1 - exp(-((t - u)*kappa)))*
           kappa)))/
    (kappa*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*sqrt(vt*vu))/
        ((sigma*sigma)*(1 - exp(-(kappa*(t - u))))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-((t - u)*kappa)))*kappa) + 
   (besseli((4*exp(-((t - u)*kappa)/2.)*sqrt(vt*vu)*
          kappa)/
        ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-(kappa*(t - u))))*
      exp(-((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma))*(-(pow(sigma,4)*(t - u))/
         (2.*pow(kappa*kappa,1.5)) + 
        ((vt + vu)*(-((pow(sigma,4)*(t - u)*
                  exp(-((t - u)*kappa)))/
                (pow(kappa,2)*
                  (1 - exp(-((t - u)*kappa))))) - 
             (pow(sigma,4)*(t - u)*
                exp(-((t - u)*kappa))*
                (1 + exp(-((t - u)*kappa))))/
              (pow(kappa,2)*
                pow(1 - exp(-((t - u)*kappa)),2)) - 
             (pow(sigma,4)*
                (1 + exp(-((t - u)*kappa))))/
              ((1 - exp(-((t - u)*kappa)))*
                pow(kappa*kappa,1.5)) + 
             (2*pow(sigma,4)*exp(-2*(t - u)*kappa)*
                ((t - u)*(t - u)))/
              (pow(1 - exp(-((t - u)*kappa)),2)*
                kappa) + 
             (pow(sigma,4)*exp(-((t - u)*kappa))*
                ((t - u)*(t - u)))/
              ((1 - exp(-((t - u)*kappa)))*
                kappa) + 
             (2*pow(sigma,4)*exp(-2*(t - u)*kappa)*
                (1 + exp(-((t - u)*kappa)))*
                ((t - u)*(t - u)))/
              (pow(1 - exp(-((t - u)*kappa)),3)*
                kappa) + 
             (pow(sigma,4)*exp(-((t - u)*kappa))*
                (1 + exp(-((t - u)*kappa)))*
                ((t - u)*(t - u)))/
              (pow(1 - exp(-((t - u)*kappa)),2)*
                kappa)))/(sigma*sigma))*
      kappa)/
    (kappa*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*sqrt(vt*vu))/
        ((sigma*sigma)*(1 - exp(-(kappa*(t - u))))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-((t - u)*kappa)))) + 
   (besseli((4*exp(-((t - u)*kappa)/2.)*sqrt(vt*vu)*
          kappa)/
        ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-(kappa*(t - u))))*
      exp(-((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma))*(((Complex(0,0.5)*(t - u)*(sigma*sigma))/
           kappa + 
          ((vt + vu)*((Complex(0,-1)*
                  exp(-((t - u)*kappa))*(t - u)*
                  (sigma*sigma))/
                (1 - exp(-((t - u)*kappa))) - 
               (Complex(0,1)*exp(-((t - u)*kappa))*
                  (1 + exp(-((t - u)*kappa)))*(t - u)*
                  (sigma*sigma))/
                pow(1 - exp(-((t - u)*kappa)),2) + 
               (Complex(0,1)*
                  (1 + exp(-((t - u)*kappa)))*
                  (sigma*sigma))/
                ((1 - exp(-((t - u)*kappa)))*
                  kappa)))/(sigma*sigma))*
        ((Complex(0,0.5)*(t - u)*(sigma*sigma))/
           kappa + 
          ((vt + vu)*((Complex(0,-1)*
                  exp(-((t - u)*kappa))*(t - u)*
                  (sigma*sigma))/
                (1 - exp(-((t - u)*kappa))) - 
               (Complex(0,1)*exp(-((t - u)*kappa))*
                  (1 + exp(-((t - u)*kappa)))*(t - u)*
                  (sigma*sigma))/
                pow(1 - exp(-((t - u)*kappa)),2) + 
               (Complex(0,1)*
                  (1 + exp(-((t - u)*kappa)))*
                  (sigma*sigma))/
                ((1 - exp(-((t - u)*kappa)))*
                  kappa)))/(sigma*sigma)))*
      kappa)/
    (kappa*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*sqrt(vt*vu))/
        ((sigma*sigma)*(1 - exp(-(kappa*(t - u))))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-((t - u)*kappa)))) + 
   ((besseli((4*exp(-((t - u)*kappa)/2.)*sqrt(vt*vu)*
            kappa)/
          ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
         (2*kappa*theta)/(sigma*sigma)) + 
        besseli((4*exp(-((t - u)*kappa)/2.)*
            sqrt(vt*vu)*kappa)/
          ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
         -2 + (2*kappa*theta)/(sigma*sigma)))*
      (1 - exp(-(kappa*(t - u))))*
      exp(-((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma))*((Complex(0,0.5)*(t - u)*(sigma*sigma))/
         kappa + 
        ((vt + vu)*((Complex(0,-1)*(t - u)*
                exp(-((t - u)*kappa))*(sigma*sigma))/
              (1 - exp(-((t - u)*kappa))) - 
             (Complex(0,1)*(t - u)*
                exp(-((t - u)*kappa))*
                (1 + exp(-((t - u)*kappa)))*
                (sigma*sigma))/
              pow(1 - exp(-((t - u)*kappa)),2) + 
             (Complex(0,1)*(1 + exp(-((t - u)*kappa)))*
                (sigma*sigma))/
              ((1 - exp(-((t - u)*kappa)))*
                kappa)))/(sigma*sigma))*
      ((Complex(0,4)*(t - u)*
           exp((-3*(t - u)*kappa)/2.)*sqrt(vt*vu))/
         pow(1 - exp(-((t - u)*kappa)),2) + 
        (Complex(0,2)*(t - u)*exp(-((t - u)*kappa)/2.)*
           sqrt(vt*vu))/(1 - exp(-((t - u)*kappa))) - 
        (Complex(0,4)*exp(-((t - u)*kappa)/2.)*
           sqrt(vt*vu))/
         ((1 - exp(-((t - u)*kappa)))*
           kappa))*kappa)/
    (kappa*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*sqrt(vt*vu))/
        ((sigma*sigma)*(1 - exp(-(kappa*(t - u))))),
       -1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-((t - u)*kappa)))) + 
   ((besseli((4*exp(-((t - u)*kappa)/2.)*sqrt(vt*vu)*
            kappa)/
          ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
         (2*kappa*theta)/(sigma*sigma)) + 
        besseli((4*exp(-((t - u)*kappa)/2.)*
            sqrt(vt*vu)*kappa)/
          ((sigma*sigma)*(1 - exp(-((t - u)*kappa)))),
         -2 + (2*kappa*theta)/(sigma*sigma)))*
      (1 - exp(-(kappa*(t - u))))*
      exp(-((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma))*((4*(t - u)*
           exp((-3*(t - u)*kappa)/2.)*(sigma*sigma)*
           sqrt(vt*vu))/
         (pow(kappa,2)*pow(1 - 
             exp(-((t - u)*kappa)),2)) + 
        (2*(t - u)*exp(-((t - u)*kappa)/2.)*
           (sigma*sigma)*sqrt(vt*vu))/
         (pow(kappa,2)*(1 - exp(-((t - u)*kappa))))\
         + (4*exp(-((t - u)*kappa)/2.)*(sigma*sigma)*
           sqrt(vt*vu))/
         ((1 - exp(-((t - u)*kappa)))*
           pow(kappa*kappa,1.5)) - 
        (8*exp((-5*(t - u)*kappa)/2.)*(sigma*sigma)*
           ((t - u)*(t - u))*sqrt(vt*vu))/
         (pow(1 - exp(-((t - u)*kappa)),3)*
           kappa) - 
        (8*exp((-3*(t - u)*kappa)/2.)*(sigma*sigma)*
           ((t - u)*(t - u))*sqrt(vt*vu))/
         (pow(1 - exp(-((t - u)*kappa)),2)*
           kappa) - 
        (exp(-((t - u)*kappa)/2.)*(sigma*sigma)*
           ((t - u)*(t - u))*sqrt(vt*vu))/
         ((1 - exp(-((t - u)*kappa)))*
           kappa))*kappa)/
    (2.*kappa*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*
          sqrt(vt*vu))/((sigma*sigma)*(1 - exp(-(kappa*(t - u)))))
        ,-1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-((t - u)*kappa)))) + 
   ((1 - exp(-(kappa*(t - u))))*
      exp(-((t - u)*(-kappa + kappa))/2. + 
        ((vt + vu)*((kappa*(1 + exp(-(kappa*(t - u)))))/
              (1 - exp(-(kappa*(t - u)))) - 
             ((1 + exp(-((t - u)*kappa)))*
                kappa)/
              (1 - exp(-((t - u)*kappa)))))/
         (sigma*sigma))*((Complex(0,4)*(t - u)*
           exp((-3*(t - u)*kappa)/2.)*sqrt(vt*vu))/
         pow(1 - exp(-((t - u)*kappa)),2) + 
        (Complex(0,2)*(t - u)*exp(-((t - u)*kappa)/2.)*
           sqrt(vt*vu))/(1 - exp(-((t - u)*kappa))) - 
        (Complex(0,4)*exp(-((t - u)*kappa)/2.)*
           sqrt(vt*vu))/
         ((1 - exp(-((t - u)*kappa)))*
           kappa))*
      (((besseli((4*exp(-((t - u)*kappa)/2.)*
                 sqrt(vt*vu)*kappa)/
               ((sigma*sigma)*
                 (1 - exp(-((t - u)*kappa)))),
              -3 + (2*kappa*theta)/(sigma*sigma)) + 
             besseli((4*exp(-((t - u)*kappa)/2.)*
                 sqrt(vt*vu)*kappa)/
               ((sigma*sigma)*
                 (1 - exp(-((t - u)*kappa)))),
              -1 + (2*kappa*theta)/(sigma*sigma)))*
           ((Complex(0,4)*(t - u)*
                exp((-3*(t - u)*kappa)/2.)*sqrt(vt*vu))
               /pow(1 - exp(-((t - u)*kappa)),2) + 
             (Complex(0,2)*(t - u)*
                exp(-((t - u)*kappa)/2.)*sqrt(vt*vu))/
              (1 - exp(-((t - u)*kappa))) - 
             (Complex(0,4)*exp(-((t - u)*kappa)/2.)*
                sqrt(vt*vu))/
              ((1 - exp(-((t - u)*kappa)))*
                kappa)))/2. + 
        ((besseli((4*exp(-((t - u)*kappa)/2.)*
                 sqrt(vt*vu)*kappa)/
               ((sigma*sigma)*
                 (1 - exp(-((t - u)*kappa)))),
              -1 + (2*kappa*theta)/(sigma*sigma)) + 
             besseli((4*exp(-((t - u)*kappa)/2.)*
                 sqrt(vt*vu)*kappa)/
               ((sigma*sigma)*
                 (1 - exp(-((t - u)*kappa)))),
              1 + (2*kappa*theta)/(sigma*sigma)))*
           ((Complex(0,4)*(t - u)*
                exp((-3*(t - u)*kappa)/2.)*sqrt(vt*vu))
               /pow(1 - exp(-((t - u)*kappa)),2) + 
             (Complex(0,2)*(t - u)*
                exp(-((t - u)*kappa)/2.)*sqrt(vt*vu))/
              (1 - exp(-((t - u)*kappa))) - 
             (Complex(0,4)*exp(-((t - u)*kappa)/2.)*
                sqrt(vt*vu))/
              ((1 - exp(-((t - u)*kappa)))*
                kappa)))/2.)*kappa)/
    (2.*kappa*besseli((4*kappa*exp(-(kappa*(t - u))/2.)*
          sqrt(vt*vu))/((sigma*sigma)*(1 - exp(-(kappa*(t - u)))))
        ,-1 + (2*kappa*theta)/(sigma*sigma))*
      (1 - exp(-((t - u)*kappa))));

      return z;
}



double quantileSquareRootProcess(double p, double kappa, double theta, double sigma, double u, double vu, double t) {
  double df = 4*kappa*theta / (sigma * sigma);
  double ncp = vu * 4*kappa*exp(-kappa*(t-u))/(sigma * sigma *(1- exp(-kappa*(t-u))));
  double k = sigma * sigma * (1-exp(-kappa*(t-u))) / (4*kappa);
  return k * qchisq(p, df, ncp);
}


double quantileIntSquareRootProcess(double p, double kappa, double theta, double sigma, double u, double vu, double t, double vt) {
  // quantile function of \int_u^t v(s) ds given v(u)=vu and v(t)=vt
//    CdfIntSquareRootProcessBridgeWithDeriv f(kappa, theta, sigma, u, vu, t, vt, p);
    CdfIntSquareRootProcessBridge g(kappa, theta, sigma, u, vu, t, vt, p);

    double m1 = std::abs(firstMomentIntegralSquareRootProcess(kappa, theta, sigma, u, vu, t, vt));
    double m2 = std::abs(secondMomentIntegralSquareRootProcess(kappa, theta, sigma, u, vu, t, vt));
        
    double stddev = std::sqrt(m2 - m1*m1);
    
    double guess = m1;
    
    if (stddev > 0.05) {
        boost::math::normal dist(m1,stddev);
        guess = quantile(dist, p);
        if (guess < 0)
            guess = m1;
    }
    
    double factor = 2;    
    boost::uintmax_t max_iter = 1000;
    std::pair<double, double> r = bracket_and_solve_root(g, guess, factor, true, 
                                                         boost::math::tools::eps_tolerance<double>(std::numeric_limits<double>::digits), max_iter);    
    int digits = DIGITS;
    
//    double intVs = boost::math::tools::newton_raphson_iterate(f, guess, 0.0, 1.0, digits);
    double intVs = (r.first + r.second) / 2.0;
    
//    printf("p:%6.4f u:%6.4f vu:%6.4f t:%6.4f vt:%6.4f m1:%6.4f m2:%6.4f intVs:%6.4f guess: %6.4f r: %6.4f\n", p, u, vu, t, vt, m1, m2, intVs, guess, r.first);
    
    return intVs;
}

double quantileHestonProcess(double p1, double p2, double p3, double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t) {
    double vt = quantileSquareRootProcess(p1, kappa, theta, sigma, u, vu, t);
    double intv = quantileIntSquareRootProcess(p2, kappa, theta, sigma, u, vu, t, vt);
    double intsqrvdW = (vt - vu - kappa*theta*(t-u) + kappa * intv) / sigma;
    double m = (r*(t-u) - 0.5*intv + rho*intsqrvdW);
    double s = sqrt((1-rho*rho)*intv);
    return Su * exp(m + s*qnorm(p3));
}
