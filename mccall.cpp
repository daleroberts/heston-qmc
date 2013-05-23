// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <stdio.h>
#include <ql/math/statistics/generalstatistics.hpp>
#include "ql_imps.h"
#include "utils.h"

using namespace std;
using QuantLib::GeneralStatistics;

Result run(double r, double rho, double kappa, double theta, double sigma, double v0, double S0, double T, double K, int discretization, int trials, int steps, int seed) {
    GeneralStatistics outerStats;
    
    for (int i = 0; i < 30; ++i) {
        double res = HestonCallMC(r, rho, kappa, theta, sigma, v0, S0, T, K, discretization, trials, steps, seed+i);
        outerStats.add(res);
    }
    
    double mean = outerStats.mean();
    double stdev = outerStats.standardDeviation();
    
    return Result(mean, stdev);
}

int main(int argc, char *argv[]) 
{
    double kappa = 6.21;
    double theta = 0.019;
    double sigma = 0.61;
    double v0 = 0.010201;
    double T = 1.0;
    double r = 0.0319;
    double S0 = 100;
    double K = 100;
    double rho = -0.7;
    
    // parse arguments
    
    long seed = argc>1 ? atol(argv[1]) : 1234;    
    int maxPowerTwo = argc>2 ? atoi(argv[2]) : 13;
    
    for (int powerTwo = 5; powerTwo <= maxPowerTwo; ++powerTwo)
    {
        long numberPaths = 1 << powerTwo;
        
        printf("%6i", numberPaths);
        for (int j = 1; j <= 6; ++j) {
            Result res = run(r, rho, kappa, theta, sigma, v0, S0, T, K, j, numberPaths, 30, seed);
            printf(" & %.5f (%.5f)", res.mean, res.stderr);
            if (j == 6) printf(" \\\\\n");
        }
        
    }
}
