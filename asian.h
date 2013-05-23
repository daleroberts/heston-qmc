// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#ifndef volmodels_asian_h
#define volmodels_asian_h

#include <stdio.h>
#include <vector>
#include <algorithm>

#include <boost/random.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "exact.h"
#include "sobolpoints.h"
#include "scramblepoints.h"
#include "sobol.h"
#include "utils.h"
#include "bridges.h"
#include "stats.h"
#include "matlab.h"

using namespace std;

Result MCAsian(double r, double rho, double kappa, double theta, double sigma, double v0, double S0, const vector<double>& times, double K, long numberPaths, long numberBatches, long prsSeed) 
{
    vector<double> outerStats;

    int tsteps = times.size();
    int tdim = tsteps-1;
    double tmax = *times.rbegin();
    
    #pragma omp parallel shared(outerStats)
    {    
        
        #pragma omp for 
        for (long k=0; k < numberBatches; k++) 
        {
            vector<double> innerStats;
            
            // new random number generator for each batch
            
            boost::mt19937 rng(prsSeed + k * numberPaths);
            boost::uniform_real<> uni_dist(0,1);
            boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);
                        
            for (unsigned long i=0; i < numberPaths; i++) 
            {
                vector<double> Ss(tsteps, S0);
                vector<double> Vs(tsteps, v0);
                
                double p;
                                
                for(long j = 1; j < tsteps; j++) 
                {
                    p = uni();
                    Ss[j] = quantileHestonProcess(p, uni(), uni(), r, rho, kappa, theta, sigma, times[j-1], Vs[j-1], Ss[j-1], times[j]);                    
                    Vs[j] = quantileSquareRootProcess(p, kappa, theta, sigma, times[j-1], Vs[j-1], times[j]);
                }
                
                double meanprice = 0.0;            
				for (std::vector<double>::iterator jj = Ss.begin() + 1; jj != Ss.end(); ++jj) {
					meanprice += *jj;
				}
                meanprice = meanprice / tdim;
                
                double payoff = meanprice - K;
                payoff = (payoff>0? payoff : 0) * exp(-r*tmax);
                
                innerStats.push_back(payoff);
            }
            
            double mean = Mean(innerStats);
            
            #pragma omp critical
            outerStats.push_back(mean);
        }
    }
    
    double mean = Mean(outerStats);
    double stderr = StdError(outerStats);
    
    return Result(mean, stderr);    
}


Result QMCAsian(double r, double rho, double kappa, double theta, double sigma, double v0, double S0, const vector<double>& times, double K, long numberPaths, long numberBatches, long prsSeed) 
{
    using boost::numeric::ublas::matrix;
    
    vector<double> outerStats;
        
    int tsteps = times.size();
    int tdim = tsteps - 1;
   	int sact = 3*tdim+1;
    double tmax = *times.rbegin();
    
    srand(prsSeed);
    int m = 16, s = 1000;
    matrix<bool> mat(m*s, m);
    getSobolpoints(mat);
    
    int Kscramble=31;
    
	if(!isPowerOfTwo(numberPaths)) {
        cout << "numberPaths not a power of two!" << endl;
        return Result(0., 0.);
    }
    int mact = log2(numberPaths);
    
    #pragma omp parallel shared(outerStats, mat)
    {    
        
        #pragma omp for 
        for (long k=0; k < numberBatches; k++) 
        {
            vector<double> innerStats;
            
            // P will contain the scrambled points
            
            matrix<double> P(sact, numberPaths); 
            
            #pragma omp critical
            getscrambledpoints(P, mat, sact, mact, m, Kscramble, numberPaths);
            
            for (long i=0; i < numberPaths; i++) 
            {
                vector<double> Ss(tsteps, S0);
                vector<double> Vs(tsteps, v0);
                
                for(long j = 1; j < tsteps; j++) 
                {					
                    Ss[j] = quantileHestonProcess(P(3*(j-1),i), P(3*(j-1)+1,i), P(3*(j-1)+2,i), r, rho, kappa, theta, sigma, times[j-1], Vs[j-1], Ss[j-1], times[j]);
                    Vs[j] = quantileSquareRootProcess(P(3*(j-1),i), kappa, theta, sigma, times[j-1], Vs[j-1], times[j]);
                    
                }
                
                double meanprice = 0.0;            
				for (std::vector<double>::iterator jj = Ss.begin() + 1; jj != Ss.end(); ++jj) {
					meanprice += *jj;
				}
                meanprice = meanprice / tdim;
                
                double payoff = meanprice - K;
                payoff = (payoff>0? payoff : 0) * exp(-r*tmax);
                
                innerStats.push_back(payoff);
            }
            
            double mean = Mean(innerStats);
            
            #pragma omp critical
            outerStats.push_back(mean);
        }
    }
    
    double mean = Mean(outerStats);
    double stderr = StdError(outerStats);
    
    return Result(mean, stderr);    
}

Result QMCAsianBridge(double r, double rho, double kappa, double theta, double sigma, double v0, double S0, const vector<double>& times, double K, long numberPaths, long numberBatches, long prsSeed) 
{
    using boost::numeric::ublas::matrix;
    using boost::numeric::ublas::matrix_column;
    
    vector<double> outerStats;
            
    int tsteps = times.size();
    int tdim = tsteps - 1;
   	int sact = 3*(tdim-1)+2*tdim+1;
    double tmax = *times.rbegin();
    
    srand(prsSeed);
    int m = 16, s = 1000;
    matrix<bool> mat(m*s, m);
    getSobolpoints(mat);
    
    int Kscramble=31;
    
	if(!isPowerOfTwo(numberPaths)) {
        cout << "numberPaths not a power of two!" << endl;
        return Result(0., 0.);
    }
    int mact = log2(numberPaths);

    #pragma omp parallel shared(outerStats, mat)
    {    
        
        #pragma omp for 
        for (int k=0; k < numberBatches; k++) 
        {
            vector<double> innerStats;
            
            // P will contain the scrambled points
           
            matrix<double> P(sact, numberPaths); 
            
            #pragma omp critical
            getscrambledpoints(P, mat, sact, mact, m, Kscramble, numberPaths);
        
            #ifdef DEBUG
                ostringstream filename;
                filename << "rep" << k+1 << ".dlm";
                ReadMatlabMatrix(filename.str(), P);
                P = boost::numeric::ublas::trans(P);
            
                matrix<double> result;
                ReadMatlabMatrix("result.dlm", result);
                result = boost::numeric::ublas::trans(result);
            #endif
            
            for (long i=0; i < numberPaths; i++) 
            {
                // extract a column from the scrambled points
                
                matrix_column<matrix<double>> Pc(P,i);
                
                // generate the value of the square root process at the time points
                
                std::vector<double> qs(Pc.begin(), Pc.begin()+3*(tdim-1)+1);
                std::vector<double> Vs = squareRootBridge(times, kappa, theta, sigma, v0, qs);
                
                // generate the integrals of V and the stochastic integrals of sqrt(V)
                
                qs.assign(Pc.begin()+3*(tdim-1)+1, Pc.begin()+3*(tdim-1)+tsteps+1);
                std::vector<double> intVs(tsteps, 0.0);
                std::vector<double> intsqrVdWs(tsteps, 0.0);
                
                for (int j=1; j < tsteps; ++j) {
                    intVs[j] = quantileIntSquareRootProcess(qs[j-1], kappa, theta, sigma, times[j-1], Vs[j-1], times[j], Vs[j]);
                    intsqrVdWs[j] = (Vs[j] - Vs[j-1] - kappa*theta*(times[j]-times[j-1]) + kappa * intVs[j]) / sigma;                    
                }

//                PIO(qs);
//                PIO(times); 
//                PIO(Vs);
//                PIO(intVs); 
//                PIO(intsqrVdWs);
                                                
                // generate the corrected drifts and vols
                
                std::vector<double> drifts(tsteps, log(S0));
                std::vector<double> vols(tsteps, 0.0);
                
                for (int j=1; j < tsteps; ++j) {
                    drifts[j] = drifts[j-1] + r*(times[j]-times[j-1]) - 0.5*intVs[j] + rho*intsqrVdWs[j];
                    vols[j] = vols[j-1] + (1-rho*rho)*intVs[j];
                }
                
                qs.assign(Pc.begin()+3*(tdim-1)+tsteps, Pc.end());
                std::vector<double> Ss = stockPriceBridge(times, drifts, vols, qs);
                
//                PIO(times);
//                PIO(drifts);
//                PIO(vols);
//                PIO(Ss);
                
                double meanprice = 0.0;            
				for (std::vector<double>::iterator jj = Ss.begin() + 1; jj != Ss.end(); ++jj) {
					meanprice += *jj;
				}
                meanprice = meanprice / tdim;
                
                double payoff = meanprice - K;
                payoff = (payoff>0? payoff : 0) * exp(-r*tmax);
                                
                innerStats.push_back(payoff);
                
//                PO(payoff);
            }
        
            double mean = Mean(innerStats);
                        
            #ifdef DEBUG
            printf("%.6f || %.6f\n", result(0,k), mean);
            #endif
            
            #pragma omp critical
            outerStats.push_back(mean);
        }
    }
    
    double mean = Mean(outerStats);
    double stderr = StdError(outerStats);
    
    return Result(mean, stderr);    
}

Result QMCAsianBridge2(double r, double rho, double kappa, double theta, double sigma, double v0, double S0, const vector<double>& times, double K, long numberPaths, long numberBatches, long prsSeed) 
{
    using boost::numeric::ublas::matrix;
    using boost::numeric::ublas::matrix_row;
    using boost::numeric::ublas::column_major;
    
    vector<double> outerStats;
        
    int tsteps = times.size();
    int tdim = tsteps - 1;
   	int sact = 3*(tdim-1)+2*tdim+1;
    double tmax = *times.rbegin();
    
    #pragma omp parallel shared(outerStats)
    {    
        
        #pragma omp for 
        for (int k=0; k < numberBatches; k++) 
        {            
            vector<double> innerStats;
            
            // P will contain the scrambled points
            
            matrix<double, column_major> P(numberPaths, sact); 
            
            #pragma omp critical
            ScrambleSobol(prsSeed + k*numberBatches, P);
            
            for (long i=0; i < numberPaths; i++) 
            {
                // extract a row from the scrambled points
                
                matrix_row<matrix<double, column_major>> Pc(P,i);
                                
                // generate the value of the square root process at the time points
                
                std::vector<double> qs(Pc.begin(), Pc.begin()+3*(tdim-1)+1);
                std::vector<double> Vs = squareRootBridge(times, kappa, theta, sigma, v0, qs);
                
                // generate the integrals of V and the stochastic integrals of sqrt(V)
                
                qs.assign(Pc.begin()+3*(tdim-1)+1, Pc.begin()+3*(tdim-1)+tsteps+1);
                
                std::vector<double> intVs(tsteps, 0.0);
                std::vector<double> intsqrVdWs(tsteps, 0.0);
                
                for (int j=1; j < tsteps; ++j) {
                    intVs[j] = quantileIntSquareRootProcess(qs[j-1], kappa, theta, sigma, times[j-1], Vs[j-1], times[j], Vs[j]);
                    intsqrVdWs[j] = (Vs[j] - Vs[j-1] - kappa*theta*(times[j]-times[j-1]) + kappa * intVs[j]) / sigma;                    
                }

                // generate the corrected drifts and vols
                
                std::vector<double> drifts(tsteps, log(S0));
                std::vector<double> vols(tsteps, 0.0);
                
                for (int j=1; j < tsteps; ++j) {
                    drifts[j] = drifts[j-1] + r*(times[j]-times[j-1]) - 0.5*intVs[j] + rho*intsqrVdWs[j];
                    vols[j] = vols[j-1] + (1-rho*rho)*intVs[j];
                }
                
                qs.assign(Pc.begin()+3*(tdim-1)+tsteps, Pc.end());
                std::vector<double> Ss = stockPriceBridge(times, drifts, vols, qs);
                
                double meanprice = 0.0;            
				for (std::vector<double>::iterator jj = Ss.begin() + 1; jj != Ss.end(); ++jj) {
					meanprice += *jj;
				}
                meanprice = meanprice / tdim;
                
                double payoff = meanprice - K;
                payoff = (payoff>0? payoff : 0) * exp(-r*tmax);
                
                innerStats.push_back(payoff);
            }
            
            double mean = Mean(innerStats);
                        
            #pragma omp critical
            outerStats.push_back(mean);
        }
    }
    
    double mean = Mean(outerStats);
    double stderr = StdError(outerStats);
    
    return Result(mean, stderr);    
}

Result QMCEurosBridge(double r, double rho, double kappa, double theta, double sigma, double v0, double S0, const vector<double>& times, double K, long numberPaths, long numberBatches, long prsSeed) 
{
    using boost::numeric::ublas::matrix;
    using boost::numeric::ublas::matrix_row;
    using boost::numeric::ublas::column_major;
    
    int tsteps = times.size();
    int tdim = tsteps - 1;
   	int sact = 3*(tdim-1)+2*tdim+1;

    vector<double> outerStats1;
    vector<double> outerStats2;
    vector<double> outerStats3;
    vector<double> outerStats4;
    
    vector<double> debugStats1;

    #pragma omp parallel shared(outerStats1,outerStats2,outerStats3,outerStats4)
    {    
        
        #pragma omp for 
        for (int k=0; k < numberBatches; k++) 
        {            
            vector<double> innerStats1;
            vector<double> innerStats2;
            vector<double> innerStats3;
            vector<double> innerStats4;
            
            // P will contain the scrambled points
            
            matrix<double, column_major> P(numberPaths, sact); 
            
            #pragma omp critical
            ScrambleSobol(prsSeed + k*numberBatches, P);
            
            for (long i=0; i < numberPaths; i++) 
            {
                // extract a row from the scrambled points
                
                matrix_row<matrix<double, column_major>> Pc(P,i);
//                vector<double> Pc = {0.9495, 0.4886, 0.1993, 0.4592, 0.9332, 0.8714, 0.7465, 0.5993, 0.9445, 0.9948, 0.9969, 0.3889, 0.0325, 0.0521, 0.2327, 0.3398, 0.4730, 0.5293 };
                
                // generate the value of the square root process at the time points
                
                std::vector<double> qs(Pc.begin(), Pc.begin()+3*(tdim-1)+1);
                
//                PIO(qs);
                
                std::vector<double> Vs = squareRootBridge(times, kappa, theta, sigma, v0, qs);
                
//                PIO(Vs);
                
                // generate the integrals of V and the stochastic integrals of sqrt(V)
                
                qs.assign(Pc.begin()+3*(tdim-1)+1, Pc.begin()+3*(tdim-1)+tsteps+1);
                
//                PIO(qs);
                
                std::vector<double> intVs(tsteps, 0.0);
                std::vector<double> intsqrVdWs(tsteps, 0.0);
                
                for (int j=1; j < tsteps; ++j) {
                    intVs[j] = quantileIntSquareRootProcess(qs[j-1], kappa, theta, sigma, times[j-1], Vs[j-1], times[j], Vs[j]);
                    intsqrVdWs[j] = (Vs[j] - Vs[j-1] - kappa*theta*(times[j]-times[j-1]) + kappa * intVs[j]) / sigma;                    
                }
                
//                PIO(intVs);
//                PIO(intsqrVdWs);

                // generate the corrected drifts and vols
                
                std::vector<double> drifts(tsteps, log(S0));
                std::vector<double> vols(tsteps, 0.0);
                
                for (int j=1; j < tsteps; ++j) {
                    drifts[j] = drifts[j-1] + r*(times[j]-times[j-1]) - 0.5*intVs[j] + rho*intsqrVdWs[j];
                    vols[j] = vols[j-1] + (1-rho*rho)*intVs[j];
                }
                                
                qs.assign(Pc.begin()+3*(tdim-1)+tsteps, Pc.end());
                
//                PIO(qs);
                
                std::vector<double> Ss = stockPriceBridge(times, drifts, vols, qs);
                
//                PIO(Ss);
                    
//                double payoff1 = Ss[1];
                double payoff1 = Ss[1] - K;            
                double payoff2 = Ss[2] - K;
                double payoff3 = Ss[3] - K;
                double payoff4 = Ss[4] - K;
                
                payoff1 = (payoff1>0? payoff1 : 0) * exp(-r*times[1]);
                payoff2 = (payoff2>0? payoff2 : 0) * exp(-r*times[2]);
                payoff3 = (payoff3>0? payoff3 : 0) * exp(-r*times[3]);
                payoff4 = (payoff4>0? payoff4 : 0) * exp(-r*times[4]);
                
                innerStats1.push_back(payoff1);
                innerStats2.push_back(payoff2);
                innerStats3.push_back(payoff3);
                innerStats4.push_back(payoff4);
                
//                if (payoff1 < 1) {
//                    printf("%6.4f : ", payoff1);
//                    for(int i = 0; i < Pc.size(); ++i) {
//                        printf("%6.4f ", Pc[i]);
//                    }
//                    printf("\n");
//                }

                double debug1 = *max_element(intVs.begin(), intVs.end());
                
                #pragma omp critical
                debugStats1.push_back(debug1);
            }
            
            double mean1 = Mean(innerStats1);
            double mean2 = Mean(innerStats2);
            double mean3 = Mean(innerStats3);
            double mean4 = Mean(innerStats4);
            

            #pragma omp critical
            {
                outerStats1.push_back(mean1);
                outerStats2.push_back(mean2);
                outerStats3.push_back(mean3);
                outerStats4.push_back(mean4);
                
            }
        }
    }
    
    PIO(debugStats1);
    
    printf("%6i  &  ", numberPaths);
    printf("%6.4f (%6.4f)  &  ", Mean(outerStats1), StdError(outerStats1));
    printf("%6.4f (%6.4f)  &  ", Mean(outerStats2), StdError(outerStats2));
    printf("%6.4f (%6.4f)  &  ", Mean(outerStats3), StdError(outerStats3));
    printf("%6.4f (%6.4f)   \\\\\n", Mean(outerStats4), StdError(outerStats4));
}


#endif
