// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <iostream>
#include <iomanip>

#include "exact.h"
#include "euro.h"

using namespace std;

int main(int argc, char *argv[]) 
{
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
   
    // parse arguments
    
    long seed = argc>1 ? atol(argv[1]) : 1234;    
    int maxPowerTwo = argc>2 ? atoi(argv[2]) : 14;

    EPS    = argc>3 ? atof(argv[3]) : 1E-12;
    STDEVS = argc>4 ? atoi(argv[4]) : 12;
    DIGITS = argc>5 ? atoi(argv[5]) : 12;
        
    for (int powerTwo = 6; powerTwo <= maxPowerTwo; ++powerTwo)
    {
        long numberPaths = 1 << powerTwo;
        
        Result mcResult = HestonCallMCExact(r, rho, kappa, theta, sigma, u, vu, Su, t, K,  30 * numberPaths, seed);
        Result qmcResult = HestonCallQMCExact(r, rho, kappa, theta, sigma, u, vu, Su, t, K,  numberPaths, 30, seed);
        Result cqmcResult = HestonCallCondQMCExact(r, rho, kappa, theta, sigma, u, vu, Su, t, K,  numberPaths, 30, seed);
        
        cout.setf(ios_base::fixed,ios_base::floatfield);
        cout << setw(12) << numberPaths << " & " 
             << setw(18) << setprecision(6) << mcResult.mean << " (" << setprecision(6) << mcResult.stderr << ")" << " &"
             << setw(18) << setprecision(6) << qmcResult.mean << " (" << setprecision(6) << qmcResult.stderr << ")" << " &"
             << setw(18) << setprecision(6) << cqmcResult.mean << " (" << setprecision(6) << cqmcResult.stderr << ")" << " \\\\" << endl;
    }
}
