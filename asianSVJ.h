// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#ifndef volmodels_asianSVJ_h
#define volmodels_asianSVJ_h

#include <iostream>
#include <vector>
#include <boost/random.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
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

vector<double> poissonPath(const vector<double>& times, double lambda, const vector<double>& qs)
{
    int tdim = times.size();
    
    vector<double> N(tdim, 0.0);
    
    for (int i=0; i < tdim-1; ++i) {
        N[i+1] = N[i] + qpois(qs[i],lambda*(times[i+1]-times[i]));
    }
    
    return N;
}

vector<double> sumPoissonPath(const vector<double>& poissonPath, double mu_s, double sigma_s, const vector<double> qs)
{
    int tdim = poissonPath.size();
    
    vector<double> S(tdim, 0.0);
    
    S[0] = poissonPath[0]*mu_s+sqrt(poissonPath[0])*sigma_s*qnorm(qs[0]);
    
    int incr, qindex = 1;
    
    for (int i = 1; i < tdim; ++i)
    {
        incr = poissonPath[i] - poissonPath[i-1];
        S[i] = S[i-1] + incr*mu_s + sqrt(incr)*sigma_s*qnorm(qs[qindex]);
        
        qindex++;
    }
    
    return S;
}

Result MCAsianSVJ(double lambda, double mu_s, double sigma_s, double r, double rho, double kappa, double theta, double sigma, double v0, double S0, const vector<double>& times, double K, long numberPaths, long numberBatches, long prsSeed) 
{
    vector<double> outerStats;
    
    int tsteps = times.size();
    int tdim = tsteps - 1;
    double tmax = *times.rbegin();
                
    #pragma omp parallel shared(outerStats)
    {    
        
        #pragma omp for
        for (long k=0; k < numberBatches; k++) 
        {
            vector<double> innerStats;
            
            // new prng for each batch
            
            boost::mt19937 rng(prsSeed  + k * numberPaths);
            boost::uniform_real<> uni_dist(0,1);
            boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);
            
            for (unsigned long i=0; i < numberPaths; i++) 
            {
                // generate a path of the square root process
                
                vector<double> qs(tdim, 0);
                for (int i = 0; i < tdim; ++i) {
                    qs[i] = uni();
                }
                
                vector<double> Vs = squareRootPath(times, kappa, theta, sigma, v0, qs);
                
                // generate the integrals of V and the stochastic integrals of sqrt(V)
                
                qs.resize(tdim);
                for (int i = 0; i < tdim; ++i) {
                    qs[i] = uni();
                }

                std::vector<double> intVs(tsteps, 0.0);
                std::vector<double> intsqrVdWs(tsteps, 0.0);
                
                for (int j=1; j < tsteps; ++j) {
                    intVs[j] = quantileIntSquareRootProcess(qs[j-1], kappa, theta, sigma, times[j-1], Vs[j-1], times[j], Vs[j]);
                    intsqrVdWs[j] = (Vs[j] - Vs[j-1] - kappa*theta*(times[j]-times[j-1]) + kappa * intVs[j]) / sigma;                    
                }
                
                qs.resize(tdim);
                for (int i = 0; i < tdim; ++i) {
                    qs[i] = uni();
                }

                double mubar = exp(mu_s + 0.5*sigma_s*sigma_s)-1;

                std::vector<double> Ss(tsteps, S0);
                
                double drift, vol;
                for (int j=1; j < tsteps; ++j) {
                    drift = log(Ss[j-1]) + (r-lambda*mubar)*(times[j]-times[j-1]) - 0.5*intVs[j] + rho*intsqrVdWs[j];
                    vol = sqrt((1-rho*rho)*intVs[j]);
                    Ss[j] = exp(drift+vol*qnorm(qs[j-1]));
                }
                
                qs.resize(tdim);
                for (int i = 0; i < tdim; ++i) {
                    qs[i] = uni();
                }

                std::vector<double> Ns = poissonPath(times, lambda, qs);
                
                qs.resize(tdim+1);
                for (int i = 0; i < tdim+1; ++i) {
                    qs[i] = uni();
                }
                
                std::vector<double> Ys = sumPoissonPath(Ns, mu_s, sigma_s, qs);                
                
                for (int j=1; j < tsteps; ++j) {
                    Ss[j] = Ss[j]*exp(Ys[j]);
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

Result QMCAsianSVJ(double lambda, double mu_s, double sigma_s, double r, double rho, double kappa, double theta, double sigma, double v0, double S0, const vector<double>& times, double K, long numberPaths, long numberBatches, long prsSeed) 
{
    using boost::numeric::ublas::matrix;
    using boost::numeric::ublas::matrix_column;

    vector<double> outerStats;
    
    srand(prsSeed);
    int m = 16, s = 1000;
    
    matrix<bool> mat(m*s, m);
    getSobolpoints(mat);
    
    int tsteps = times.size();
    int tdim = tsteps - 1;
    int tm = log2(tsteps-1);
   	int sact = 5*tdim + 1;
    double tmax = *times.rbegin();

    int Kscramble=31;
    
	if(!isPowerOfTwo(numberPaths)) {
        cout << "numberPaths not a power of two!" << endl;
        return Result(0.0, 0.0);
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
            
            // this is terrible, each thread has to lock and wait
            
            #pragma omp critical
            getscrambledpoints(P, mat, sact, mact, m, Kscramble, numberPaths);
            
            for (long i=0; i < numberPaths; i++) 
            {
                // extract a column from the scrambled points
                
                matrix_column<matrix<double>> Pc(P,i);
                
                // generate the value of the square root process at the time points
                
                std::vector<double> qs(Pc.begin(), Pc.begin()+tdim);
                
                vector<double> Vs = squareRootPath(times, kappa, theta, sigma, v0, qs);
                
                // generate the integrals of V and the stochastic integrals of sqrt(V)
                
                qs.assign(Pc.begin()+tdim, Pc.begin()+2*tdim);
                
                std::vector<double> intVs(tsteps, 0.0);
                std::vector<double> intsqrVdWs(tsteps, 0.0);
                
                for (int j=1; j < tsteps; ++j) {
                    intVs[j] = quantileIntSquareRootProcess(qs[j-1], kappa, theta, sigma, times[j-1], Vs[j-1], times[j], Vs[j]);
                    intsqrVdWs[j] = (Vs[j] - Vs[j-1] - kappa*theta*(times[j]-times[j-1]) + kappa * intVs[j]) / sigma;                    
                }

                // Share price
                
                qs.assign(Pc.begin()+2*tdim, Pc.begin()+3*tdim);
                
                double mubar = exp(mu_s + 0.5*sigma_s*sigma_s)-1;
                
                std::vector<double> Ss(tsteps, S0);
                
                double drift, vol;
                for (int j=1; j < tsteps; ++j) {
                    drift = log(Ss[j-1]) + (r-lambda*mubar)*(times[j]-times[j-1]) - 0.5*intVs[j] + rho*intsqrVdWs[j];
                    vol = sqrt((1-rho*rho)*intVs[j]);
                    Ss[j] = exp(drift+vol*qnorm(qs[j-1]));
                }

                // Poisson path
                
                qs.assign(Pc.begin()+3*tdim, Pc.begin()+4*tdim);
                
                std::vector<double> Ns = poissonPath(times, lambda, qs);
                
                // Compound path
                
                qs.assign(Pc.begin()+4*tdim, Pc.begin()+5*tdim+1);
                
                std::vector<double> Ys = sumPoissonPath(Ns, mu_s, sigma_s, qs);                
                
                for (int j=1; j < tsteps; ++j) {
                    Ss[j] = Ss[j]*exp(Ys[j]);
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

Result QMCAsianBridgeSVJ(double lambda, double mu_s, double sigma_s, double r, double rho, double kappa, double theta, double sigma, double v0, double S0, const vector<double>& times, double K, long numberPaths, long numberBatches, long prsSeed) 
{
    using boost::numeric::ublas::matrix;
    using boost::numeric::ublas::matrix_column;
        
    vector<double> outerStats;
    
    int tsteps = times.size();
    int tdim = tsteps - 1;
   	int sact = 1+3*(tdim-1)+4*tdim;
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
                        
            // this is terrible, each thread has to lock and wait
            
            #pragma omp critical
            getscrambledpoints(P, mat, sact, mact, m, Kscramble, numberPaths);
            
//            ostringstream filename;
//            filename << "rep" << k+1 << ".dlm";
//            ReadMatlabMatrix(filename.str(), P);
//            P = boost::numeric::ublas::trans(P);
            
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
                
                // generate the corrected drifts and vols
                
                std::vector<double> drifts(tsteps, log(S0));
                std::vector<double> vols(tsteps, 0.0);
                
                double mubar = exp(mu_s + 0.5*sigma_s*sigma_s)-1;
                
                for (int j=1; j < tsteps; ++j) {
                    drifts[j] = drifts[j-1] + (r-lambda*mubar)*(times[j]-times[j-1]) - 0.5*intVs[j] + rho*intsqrVdWs[j];
                    vols[j] = vols[j-1] + (1-rho*rho)*intVs[j];
                }
                
                qs.assign(Pc.begin()+3*(tdim-1)+tdim+1, Pc.begin()+(1+3*(tdim-1)+2*tdim));
                std::vector<double> Ss = stockPriceBridge(times, drifts, vols, qs);
                
                qs.assign(Pc.begin()+1+3*(tdim-1)+2*tdim, Pc.begin()+1+3*(tdim-1)+3*tdim);
                std::vector<double> Ns = poissonPathBridge(times, lambda, qs);
                                
                qs.assign(Pc.begin()+(1+3*(tdim-1)+3*tdim), Pc.end());
                std::vector<double> Js = sumPoissonPathBridge(Ns, mu_s, sigma_s, qs);
                
                for (int i = 0; i < Ss.size(); ++i) {
                    Ss[i] = Ss[i] * exp(Js[i]);
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
            
//            printf("%i %.6f\n", k+1, mean);
            
            #pragma omp critical
            outerStats.push_back(mean);
        }
    }
    
    double mean = Mean(outerStats);
    double stderr = StdError(outerStats);
    
    return Result(mean, stderr);    
}

#endif
