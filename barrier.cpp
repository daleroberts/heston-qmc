// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <iostream>
#include <iomanip>
#include <vector>
#include <boost/random.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <ql/quantlib.hpp>
#include "euro.h"
#include "exact.h"
#include "sobolpoints.h"
#include "scramblepoints.h"
#include "utils.h"

using namespace std;
using QuantLib::GeneralStatistics;

typedef boost::numeric::ublas::matrix<bool> BooleanMatrix;
typedef boost::numeric::ublas::matrix<double> DoubleMatrix;

double QMCBarrierDownOutCallNaive(double r, double rho, double kappa, double theta, double sigma, double v0, double S0, double T, double K, double H, long numberPaths, long numberBatches, long prsSeed=123, long numberTimeMonitors = 32) 
{
    // Scrambled Sobol version - crappy parallel version using a critical section
    
    GeneralStatistics outerStats;
    
    double dt = T / numberTimeMonitors;
    
    srand(prsSeed);
    int m = 16, s = 1000;
    
    BooleanMatrix mat(m*s, m);
    getSobolpoints(mat);
    
   	int sact = 3 * numberTimeMonitors;
    int Kscramble=31;

	if(!isPowerOfTwo(numberPaths)) {
        cout << "numberPaths not a power of two!" << endl;
        return 0.0;
    }
    
    int mact = (int) log2(numberPaths);
        
    #ifdef _OPENMP

    #pragma omp parallel shared(outerStats, mat)
    {    
        #pragma omp for 
        for (long k=0; k < numberBatches; k++) {
            GeneralStatistics innerStats;
            
            //P will contain the scrambled points
           	DoubleMatrix P(sact, numberPaths); 
            
            // this is terrible, each thread has to lock and wait
            #pragma omp critical
            getscrambledpoints(P, mat, sact, mact, m, Kscramble, numberPaths);          
            
            for (long i=0; i < numberPaths; i++) 
            {
                double payoff;            
                bool alive = true;
                
                double t = dt, u = 0.0, vu = v0, Su = S0, St;
                double vt, intv, intsqrvdW;
                
                for (long j=0; j < numberTimeMonitors; j++) 
                {                    
                    vt = quantileSquareRootProcess(P(3*j,i), kappa, theta, sigma, u, vu, t);
                    intv = quantileIntSquareRootProcess(P(3*j+1,i), kappa, theta, sigma, u, vu, t, vt);
                    intsqrvdW = (vt - vu - kappa*theta*(t-u) + kappa * intv) / sigma;
                    St = Su * exp((r*(t-u) - 0.5*intv + rho*intsqrvdW) + sqrt((1-rho*rho)*intv)*qnorm(P(3*j+2,i)));

                    if (St <= H) 
                        alive = false;
                    
                    // save previous
                    Su = St;
                    vu = vt;
                    u = t;
                    
                    // increment time
                    t += dt;
                }
                
                payoff = St - K;
                payoff = (payoff>0 && alive ? payoff : 0) * exp(-r*(t-u));
                  
                innerStats.add(payoff);
            }
        
            #pragma omp critical
            outerStats.add(innerStats.mean());
        }
    }
    
    double mean = outerStats.mean();
    double sds = outerStats.errorEstimate();
    double stdev = outerStats.standardDeviation();
    
    cout.setf(ios_base::fixed,ios_base::floatfield);
    std::cout << setw(12) << setprecision(6) << mean << setw(12) << setprecision(6) << stdev;
    std::cout << std::endl;

    return mean;
    
    #else
    
    cout << "No OpenMP!" << endl;
    return 0.0;
    
    #endif
}
