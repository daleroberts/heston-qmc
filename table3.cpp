// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <iostream>
#include <iomanip>

#include "exact.h"
#include "asianSVJ.h"

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

    int tsteps = argc > 1 ? atoi(argv[1]) : 4;
//    int tsteps = 16;
    int m = log2(tsteps);
    double dt = 1.0/tsteps;
    
    vector<double> times;
    for (int i = 0; i <= tsteps; i++)
        times.push_back(i*dt);

    printf("tsteps: %i\n", tsteps);
    
    // parse arguments
    
    int seed = 123;
    int maxPowerTwo = argc>2 ? atoi(argv[2]) : 14;
    
    EPS    = argc>3 ? atof(argv[3]) : 1E-12;
    STDEVS = argc>4 ? atoi(argv[4]) : 12;
    DIGITS = argc>5 ? atoi(argv[5]) : 12;
    
    for (int powerTwo = 13; powerTwo <= maxPowerTwo; ++powerTwo)
    {
        long numberPaths = 1 << powerTwo;
        int numberBatches = 30;
        
        Result mcResult = MCAsianSVJ(lambda, mu_s, sigma_s, r, rho, kappa, theta, sigma, vu, Su, times, K, numberPaths, numberBatches, seed);
        Result qmcResult = QMCAsianSVJ(lambda, mu_s, sigma_s, r, rho, kappa, theta, sigma, vu, Su, times, K, numberPaths, numberBatches, seed);
        Result bqmcResult = QMCAsianBridgeSVJ(lambda, mu_s, sigma_s, r, rho, kappa, theta, sigma, vu, Su, times, K, numberPaths, numberBatches, seed);
        
        cout.setf(ios_base::fixed,ios_base::floatfield);
        cout << setw(12) << numberPaths << " & " 
        << setw(18) << setprecision(6) << mcResult.mean << " (" << setprecision(6) << mcResult.stderr << ")" << " &"
        << setw(18) << setprecision(6) << qmcResult.mean << " (" << setprecision(6) << qmcResult.stderr << ")" << " &"
        << setw(18) << setprecision(6) << bqmcResult.mean << " (" << setprecision(6) << bqmcResult.stderr << ")" << " \\\\" << endl;
    }
}
