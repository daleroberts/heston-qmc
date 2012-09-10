#include <iostream>
#include <complex>

#include "exact.h"
#include "toms644.h"

typedef std::complex<double> Complex;
typedef double Real;

template<class T>
std::complex<T> besseli(const std::complex<T>& z, const T& nu) 
{
    double zr = z.real(), zi = z.imag(), fnu = nu, cyr = 0.0, cyi = 0.0;
    int kode = 1, n = 1, nz = 0, ierr = 0;
    Complex r;

    if (fnu < 0) {
	    //I(-nu,z) = I(nu,z) + (2/pi)*sin(pi*nu)*K(nu,z)
        fnu = -fnu;
        zbesi_(&zr, &zi, &fnu, &kode, &n, &cyr, &cyi, &nz, &ierr);
        r = Complex(cyr, cyi);
        zbesk_(&zr, &zi, &fnu, &kode, &n, &cyr, &cyi, &nz, &ierr);
        r = r + 2/PI*sin(PI*fnu)*Complex(cyr, cyi);
    } else {
        zbesi_(&zr, &zi, &fnu, &kode, &n, &cyr, &cyi, &nz, &ierr);
        r = Complex(cyr, cyi);
    }

    return r;
}


Complex besseli(const Real& z, const Real& nu) 
{
    double zr = z, zi = 0.0, fnu = nu, cyr = 0.0, cyi = 0.0;
    int kode = 1, n = 1, nz = 0, ierr = 0;
    Complex r;

    if (fnu < 0) {
	    //I(-nu,z) = I(nu,z) + (2/pi)*sin(pi*nu)*K(nu,z)
        fnu = -fnu;
        zbesi_(&zr, &zi, &fnu, &kode, &n, &cyr, &cyi, &nz, &ierr);
        r = Complex(cyr, cyi);
        zbesk_(&zr, &zi, &fnu, &kode, &n, &cyr, &cyi, &nz, &ierr);
        r = r + 2/PI*sin(PI*fnu)*Complex(cyr, cyi);
    } else {
        zbesi_(&zr, &zi, &fnu, &kode, &n, &cyr, &cyi, &nz, &ierr);
        r = Complex(cyr, cyi);
    }

    return r;
}


Complex gamma(const Real& kappa, const Real& sigma, const Real& a) 
{
    return sqrt(kappa*kappa-2*sigma*sigma*a*Complex(0,1));
}

Complex intsqrm1(Real kappa, Real theta, Real sigma, Real u, Real vu, Real t, Real vt) 
{
    using namespace std;

    Real d = 4*theta*kappa/(sigma*sigma);

    Complex z = Complex(0,-1)*pow(kappa,-2)*pow(sigma,2) 
        + Complex(0,0.5)*(t - u)*pow(kappa,-1)*pow(sigma,2) 
        + Complex(0,1)*(t - u)*exp(kappa*(-t + u))*pow(kappa,-1)
        * pow(sigma,2)*pow(1 - exp(kappa*(-t + u)),-1);
        
    z += (vt + vu)*pow(sigma,-2)*(Complex(0,-1)*(t - u)
        * exp(kappa*(-t + u)) * (1 + exp(kappa*(-t + u)))*pow(sigma,2)
        * pow(1 - exp(kappa*(-t + u)),-2) 
        - Complex(0,1)*(t - u)*exp(kappa*(-t + u))
        * pow(sigma,2)* pow(1 - exp(kappa*(-t + u)),-1) 
        + Complex(0,1)*(1 + exp(kappa*(-t + u)))
        * pow(kappa,-1)* pow(sigma,2)
        * pow(1 - exp(kappa*(-t + u)),-1));

    Complex t1 = Complex(0,4)*(t - u)*exp((-3*kappa*(t - u))/2.)
        * pow(1 - exp(kappa*(-t + u)),-2)*sqrt(vt)*sqrt(vu) 
        + Complex(0,2)*(t - u)*exp(-(kappa*(t - u))/2.) 
        * pow(1 - exp(kappa*(-t + u)),-1)*sqrt(vt)*sqrt(vu) 
        - Complex(0,4)*exp(-(kappa*(t - u))/2.)*pow(kappa,-1)
        * pow(1 - exp(kappa*(-t + u)),-1)*sqrt(vt)*sqrt(vu);

    Complex t2 = besseli(4*kappa*exp(-(kappa*(t - u))/2.)
        * pow(sigma,-2)* pow(1 - exp(kappa*(-t + u)),-1)
        * sqrt(vt)*sqrt(vu),-2 + d/2.) 
        + besseli(4*kappa*exp(-(kappa*(t - u))/2.)
        * pow(sigma,-2)* pow(1 - exp(kappa*(-t + u)),-1)
        * sqrt(vt)*sqrt(vu),d/2.);

    Complex t3 = Complex(2,0)*besseli(4*kappa*exp(-(kappa*(t - u))/2.)
        * pow(sigma,-2) * pow(1 - exp(kappa*(-t + u)),-1)
        * sqrt(vt)*sqrt(vu),-1 + d/2.);

    z += t1*t2/t3;

    return z;
}

Complex intsqrm2(Real kappa, Real theta, Real sigma, Real u, Real vu, Real t, Real vt) {
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

/*
RcppExport SEXP BesselI(SEXP z_, SEXP nu_) {
BEGIN_RCPP
    using namespace Rcpp;

    Rcomplex z = as<Rcomplex>(z_);
    double nu = as<double>(nu_);

    Complex x = besseli(Complex(z.r, z.i), nu);

    z.r = x.real();
    z.i = x.imag();

    return wrap(z);
END_RCPP
}
*/

SEXP pintsqr_c(SEXP kappa_, SEXP theta_, SEXP sigma_, SEXP u_, SEXP vu_, SEXP t_, SEXP vt_, SEXP x_) {
//BEGIN_RCPP
    using namespace Rcpp;
    using namespace std;

    Real kappa = as<Real>(kappa_);
    Real theta = as<Real>(theta_);
    Real sigma = as<Real>(sigma_);
    Real u     = as<Real>(u_);
    Real vu    = as<Real>(vu_);
    Real t     = as<Real>(t_);
    Real vt    = as<Real>(vt_);
    Real x     = as<Real>(x_);

    int m = 7;
    Real eps = 10E-7;
    int k = 1;
    int add2pi = 0;
 
    Real delta = 4*theta*kappa/(sigma*sigma);
    Complex m1 = intsqrm1(kappa, theta, sigma, u, vu, t, vt);
    Complex m2 = intsqrm2(kappa, theta, sigma, u, vu, t, vt);

    Real ueps = abs(m1) + m*abs(sqrt(m2-m1*m1));
    Real h = 2*PI / (x+ueps);

    Complex b1 = besseli(4*kappa*sqrt(vu*vt)*exp(-0.5*kappa*(t-u))/((sigma*sigma)*(1-exp(-kappa*(t-u)))), 0.5*delta-1);

    Complex ga = gamma(kappa, sigma, h*k);
    Complex zc = sqrt(vu*vt)*4*ga*exp(-0.5*ga*(t-u))/((sigma*sigma)*(Complex(1,0)-exp(-ga*(t-u))));

    Real angle = 0.0;
    Real oldangle = atan2(zc.imag(), zc.real());

    Complex b2 = exp(add2pi*(0.5*delta-1)*2*PI*Complex(0,1))*besseli(zc, 0.5*delta-1);
    Complex t1 = ga*exp(-0.5*(ga-kappa)*(t-u))*(1-exp(-kappa*(t-u)))/(kappa*(Complex(1,0)-exp(-ga*(t-u))));
    Complex t2 = exp((vu+vt)/(sigma*sigma)*(kappa*(Complex(1,0)+exp(-kappa*(t-u)))/(1-exp(-kappa*(t-u)))-ga*(Complex(1,0)+exp(-ga*(t-u)))/(Complex(1,0)-exp(-ga*(t-u))) ));
    Real    cc = (h*x/PI)+(2/PI)*sin(h*k*x)*(b2*t1*t2).real()/(b1.real()*k);

    Real tv = abs(b2*t1*t2/(b1*Complex(k,0)));

    while (abs(b2*t1*t2/(b1*Complex(k,0))) >= PI*eps/2) {
        if (k > 10000) {
            //cout << "k reached max" << endl;
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
        t2 = exp((vu+vt)/(sigma*sigma)*(kappa*(Complex(1,0)+exp(-kappa*(t-u)))/(1-exp(-kappa*(t-u)))-ga*(Complex(1,0)+exp(-ga*(t-u)))/(Complex(1,0)-exp(-ga*(t-u))) ));

        cc += (2/PI)*sin(h*k*x)*(b2*t1*t2).real()/(b1.real()*k);
    }

    return wrap(cc);
//END_RCPP
}

