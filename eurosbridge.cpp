// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <stdio.h>

#include "exact.h"
#include "asian.h"
#include "ql_imps.h"

using namespace std;

int main(int argc, char *argv[]) 
{
    double lambda = 0.11;
	double mu_s = -0.1391;
	double sigma_s = 0.15;
    
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
    
    int tsteps = 4;
    int m = log2(tsteps);
    double dt = 1.0/tsteps;
    
    vector<double> times;
    for (int i = 0; i <= tsteps; i++)
        times.push_back(i*dt);
    
    // parse arguments
    
    long seed = argc>1 ? atol(argv[1]) : 4711;    
    int maxPowerTwo = argc>2 ? atoi(argv[2]) : 10;
    
    EPS    = argc>3 ? atof(argv[3]) : 1E-12;
    STDEVS = argc>4 ? atoi(argv[4]) : 12;
    DIGITS = argc>5 ? atoi(argv[5]) : 12;
    
    vector<double> value = {0.0, 2.67093, 4.25452, 5.6008, 6.8061};
    
    printf("        &");
    
    for (int i = 1; i < times.size(); ++i) {
        printf("  %6.4f           ", value[i]);
        if (i != 4) printf("&");
        else printf("\\\\\n");
    }
        
//    for (int powerTwo = 5; powerTwo <= maxPowerTwo; ++powerTwo)
//    {
//        long numberPaths = 1 << powerTwo;
//        int numberBatches = 30;
//                
//        QMCEurosBridge(r, rho, kappa, theta, sigma, vu, Su, times, K, numberPaths, numberBatches, seed);
//    }
    
    QMCEurosBridge(r, rho, kappa, theta, sigma, vu, Su, times, K, 64, 30, seed);
    
}
