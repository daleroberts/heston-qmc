// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <stdio.h>

#include "exact.h"
#include "asian.h"

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
        
    EPS    = 1E-12;
    STDEVS = 12;
    DIGITS = 12;

    // parse arguments
    
    int seed = argc>1 ? atol(argv[1]) : 1234;    
    int powerTwo = argc>2 ? atoi(argv[2]) : 5;
    int max_tsteps = argc>3 ? atoi(argv[3]) : 4;
    
    int numberPaths = 1 << powerTwo;
    int numberBatches = 30;
    
    for (int tsteps = 4; tsteps <= max_tsteps; ++tsteps) 
    {
        double dt = t/tsteps;
        
        vector<double> times;
        for (int i = 0; i <= tsteps; i++)
            times.push_back(i*dt);

        Result naive = QMCAsian(r, rho, kappa, theta, sigma, vu, Su, times, K, numberPaths, numberBatches, seed);   
        Result bridge = QMCAsianBridge(r, rho, kappa, theta, sigma, vu, Su, times, K, numberPaths, numberBatches, seed);
        
        printf("%4i  &  %.6f (%.6f)  &   %.6f (%.6f)  \\\\\n", tsteps, naive.mean, naive.stderr, bridge.mean, bridge.stderr);
        
    }
}
