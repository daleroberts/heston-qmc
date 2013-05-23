// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <iostream>
#include "ql_imps.h"

int main(int argc, char *argv[]) {
    double kappa = 6.21;
    double theta = 0.019;
    double sigma = 0.61;
    double u = 0;
    double vu = 0.010201;
    double t = 1.0;
    double r = 0.0319;
    double Su = 100;
    double K = 100;
    double rho = -0.7;

    double value = HestonCallAnalytic(r, rho, kappa, theta, sigma, vu, Su, (t-u), K);
    
    std::cout << value << std::endl;
}
