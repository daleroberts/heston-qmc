// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <iostream>
#include <iomanip>
#include <vector>

#include <boost/random.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

#include <ql/math/statistics/generalstatistics.hpp>

#include "euro.h"
#include "exact.h"
#include "sobolpoints.h"
#include "scramblepoints.h"
#include "utils.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

typedef boost::numeric::ublas::matrix<bool> BooleanMatrix;
typedef boost::numeric::ublas::matrix<double> DoubleMatrix;
typedef boost::numeric::ublas::vector<double> DoubleVector;

using namespace std;
using QuantLib::GeneralStatistics;

Result HestonCallMCExact(double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t, double Strike, long NumberOfPaths, long prsSeed=1) {
    
    double stock = 0.0;
    double mean = 0.0;
    double stddev = 0.0;
    double thisPayoff = 0.0;
    double runningSum = 0.0;
    double runningdlSum = 0.0;
    
    int nthreads = 1;
    long npaths = NumberOfPaths;
    long nOfPathsForThread = npaths;

    #pragma omp parallel reduction(+:runningSum) reduction(+:runningdlSum)
    {
        #ifdef _OPENMP
        long worker_seed = prsSeed + omp_get_thread_num() * NumberOfPaths;
        #else
        long worker_seed = prsSeed;
        #endif
        
        boost::mt19937 rng(worker_seed);

        boost::uniform_real<> uni_dist(0,1);
        boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);

        #pragma omp master
        {
            #ifdef _OPENMP
            nthreads = omp_get_num_threads();
            nOfPathsForThread = long(NumberOfPaths / nthreads);
            npaths = nthreads * nOfPathsForThread;
            #else
            nthreads = 1;
            #endif            
        }
        
        for (long i=0; i< nOfPathsForThread; i++) {
            
            stock = quantileHestonProcess(uni(), uni(), uni(), r, rho, kappa, theta, sigma, u, vu, Su, t);

            thisPayoff = stock - Strike;
            thisPayoff = (thisPayoff>0 ? thisPayoff : 0)*exp(-r*(t-u));
            runningSum += thisPayoff;
            runningdlSum += thisPayoff*thisPayoff;
        }
    }

    mean = runningSum / (double) npaths;
    stddev = sqrt((runningdlSum-mean*runningSum)/(npaths*(npaths-1)));
    
    return Result(mean, stddev);
}

Result HestonCallCondQMCExact(double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t, double K, 
                        long numberPaths, long numberBatches, long prsSeed=1234) {
    // Scrambled Sobol version - crappy parallel version using a critical section
    
    GeneralStatistics outerStats;
    
    srand(prsSeed);
    int m = 16, s = 1000;
    BooleanMatrix mat(m*s, m);
    getSobolpoints(mat);
    
   	int sact=2; //actual dimensionality is sact    
    int Kscramble=31; //number of digits scrambled

	if(!isPowerOfTwo(numberPaths)) {
        cout << "numberPaths not a power of two!" << endl;
        return Result(0.,0.);
    }
    
    int mact = (int) log2(numberPaths);
        
    #pragma omp parallel shared(outerStats, mat)
    {    
        #pragma omp for 
        for (long k=0; k < numberBatches; k++) {
            GeneralStatistics innerStats;
            vector<double> p;
            double payoff;
            
            //P will contain the scrambled points
           	DoubleMatrix P(sact, numberPaths); 
            
            // this is terrible, each thread has to lock and wait
            #pragma omp critical
            getscrambledpoints(P, mat, sact, mact, m, Kscramble, numberPaths);          
            
            for (long i=0; i < numberPaths; i++) {
                
                double vt = quantileSquareRootProcess(P(0,i), kappa, theta, sigma, u, vu, t);
                double intv = quantileIntSquareRootProcess(P(1,i), kappa, theta, sigma, u, vu, t, vt);
                double intsqrvdW = (vt - vu - kappa*theta*(t-u) + kappa * intv) / sigma;
                
                double xi = exp(-0.5*rho*rho*intv + rho*intsqrvdW);
                double sigmatilde = sqrt(intv/(t-u));
                                
                payoff = bscall(Su*xi, K, u, t, r, sigmatilde*sqrt(1-rho*rho));

                innerStats.add(payoff);
            }
        
            #pragma omp critical
            outerStats.add(innerStats.mean());
        }
    }
    
    double mean = outerStats.mean();
    double sds = outerStats.errorEstimate();
    double stdev = outerStats.standardDeviation();
    
    return Result(mean, stdev);
}

Result HestonCallQMCExact(double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t, double K, 
                              long numberPaths, long numberBatches, long prsSeed=1234) {
    // Scrambled Sobol version - crappy parallel version using a critical section
    
    GeneralStatistics outerStats;
    
    srand(prsSeed);
    int m = 16, s = 1000;
    BooleanMatrix mat(m*s, m);
    getSobolpoints(mat);
    
   	int sact=3; //actual dimensionality is sact    
    int Kscramble=31; //number of digits scrambled
    
	if(!isPowerOfTwo(numberPaths)) {
        cout << "numberPaths not a power of two!" << endl;
        return Result(0.,0.);
    }
    
    int mact = (int) log2(numberPaths);
    
    #pragma omp parallel shared(outerStats, mat)
    {    
        #pragma omp for 
        for (long k=0; k < numberBatches; k++) {
            GeneralStatistics innerStats;
            vector<double> p;
            double payoff;
            
            //P will contain the scrambled points
           	DoubleMatrix P(sact, numberPaths); 
            
            // this is terrible, each thread has to lock and wait
            #pragma omp critical
            getscrambledpoints(P, mat, sact, mact, m, Kscramble, numberPaths);          
            
            for (long i=0; i < numberPaths; i++) {
                payoff = quantileHestonProcess(P(0,i), P(1,i), P(2,i), r, rho, kappa, theta, sigma, u, vu, Su, t) - K;
                payoff = (payoff>0 ? payoff : 0)*exp(-r*(t-u));
                
                innerStats.add(payoff);
            }
            
            #pragma omp critical
            outerStats.add(innerStats.mean());
        }
    }
    
    double mean = outerStats.mean();
    double sds = outerStats.errorEstimate();
    double stdev = outerStats.standardDeviation();
    
    return Result(mean, stdev);
}

//double HestonCallMCExact(double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t, double K, long numberPaths, long numberBatches, long prsSeed=1) {
//    
//    boost::mt19937 rng(prsSeed);
//
//    boost::uniform_real<> uni_dist(0,1);
//    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);
//    
//    GeneralStatistics outerStats;
//    std::vector<double> p;
//    
//    for (long k=0; k < numberBatches; k++) {
//        GeneralStatistics innerStats;
//        double payoff;
//
//        for (long i=0; i < numberPaths; i++) {
//            payoff = quantileHestonProcess(uni(), uni(), uni(), r, rho, kappa, theta, sigma, u, vu, Su, t) - K;
//            payoff = (payoff>0 ? payoff : 0)*exp(-r*(t-u));
//
//            innerStats.add(payoff);
//        }
//        
//        outerStats.add(innerStats.mean());
//    }
//    
//    double mean = outerStats.mean();
//    double sds = outerStats.errorEstimate();
//    
//    std::cout << "mean: " << mean << " errorEstimate: " << sds << " paths: " << numberPaths << " batches: " << numberBatches << std::endl;
//
//    return mean;
//}

//double HestonCallQMCExact(double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t, double K, long numberPaths, long numberBatches, long ldsSeed=1, long prsSeed=1) {
//    
//    srand(prsSeed);
//    int m = 16, s = 1000;
//    BooleanMatrix mat(m*s, m);
//    getSobolpoints(mat);
//    
//  	int sact=3; //actual dimensionality is sact
//    
//    int Kscramble=31; //number of digits scrambled
//
//	if(!isPowerOfTwo(numberPaths)) {
//        std::cout << "numberPaths not a power of two!" <<std::endl;
//        return 0.0;
//    }
//    
//    int mact = (int) log2(numberPaths);
//    
//	DoubleMatrix P(sact, numberPaths); //P will contain the scrambled points	
//	DoubleVector est(numberBatches); //to contain the estimate for 1 repetition
//    
//    GeneralStatistics outerStats;
//    std::vector<double> p;
//    
//    for (long k=0; k < numberBatches; k++) {
//        GeneralStatistics innerStats;
//        double payoff;
//        
//  		getscrambledpoints(P, mat, sact, mact, m, Kscramble, numberPaths);
//
//        for (long i=0; i < numberPaths; i++) {
//            payoff = quantileHestonProcess(P(0,i), P(1,i), P(2,i), r, rho, kappa, theta, sigma, u, vu, Su, t) - K;
//            payoff = (payoff>0 ? payoff : 0)*exp(-r*(t-u));
//
//            innerStats.add(payoff);
//        }
//        
//        outerStats.add(innerStats.mean());        
//    }
//    
//    double mean = outerStats.mean();
//    double sds = outerStats.errorEstimate();
//    double stdev = outerStats.standardDeviation();
//    
//    cout.setf(ios_base::fixed,ios_base::floatfield);
////    std::cout << " Mean: " << setprecision(6) << mean << " stdDev: " << setprecision(6) << stdev;
//    std::cout << setw(12) << setprecision(6) << mean << setw(12) << setprecision(6) << stdev;
//
//    return mean;
//}
//
//double HestonCallQMCExact2(double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t, double K, long numberPaths, long numberBatches, long ldsSeed=1, long prsSeed=1) {
//    
//    srand(prsSeed);
//    int m = 16, s = 1000;
//    BooleanMatrix mat(m*s, m);
//    getSobolpoints(mat);
//    
//  	int sact=2; //actual dimensionality is sact
//    
//    int Kscramble=31; //number of digits scrambled
//
//	if(!isPowerOfTwo(numberPaths)) {
//        cout << "numberPaths not a power of two!" << endl;
//        return 0.0;
//    }
//    
//    int mact = (int) log2(numberPaths);
//    
//	DoubleMatrix P(sact, numberPaths); //P will contain the scrambled points	
//	DoubleVector est(numberBatches); //to contain the estimate for 1 repetition
//    
//    GeneralStatistics outerStats;
//    vector<double> p;
//    
//    for (long k=0; k < numberBatches; k++) {
//        GeneralStatistics innerStats;
//        double payoff;
//        
//  		getscrambledpoints(P, mat, sact, mact, m, Kscramble, numberPaths);
//
//        for (long i=0; i < numberPaths; i++) {
//            double vt = quantileSquareRootProcess(P(0,i), kappa, theta, sigma, u, vu, t);
//            double intv = quantileIntSquareRootProcess(P(1,i), kappa, theta, sigma, u, vu, t, vt);
//            double intsqrvdW = (vt - vu - kappa*theta*(t-u) + kappa * intv) / sigma;
//            
//            double xi = exp(-0.5*rho*rho*intv + rho*intsqrvdW);
//            double sigmatilde = sqrt(intv/(t-u));
//                            
//            payoff = bscall(Su*xi, K, u, t, r, sigmatilde*sqrt(1-rho*rho));
//
//            innerStats.add(payoff);
//        }
//        
//        outerStats.add(innerStats.mean());        
//    }
//    
//    double mean = outerStats.mean();
//    double sds = outerStats.errorEstimate();
//    double stdev = outerStats.standardDeviation();
//    
//    cout.setf(ios_base::fixed,ios_base::floatfield);
////    std::cout << " Mean: " << setprecision(6) << mean << " stdDev: " << setprecision(6) << stdev;
//    std::cout << setw(12) << setprecision(6) << mean << setw(12) << setprecision(6) << stdev;
//
//    return mean;
//}

//double HestonCallQMCExactOpenMP(double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t, double K, long numberPaths, long numberBatches, long ldsSeed=0, long prsSeed=1) {
//    
//    using namespace QuantLib;
//    
//    GeneralStatistics outerStats;
//
//    #ifdef _OPENMP
//    
//    #pragma omp parallel shared(outerStats)
//    {
//        #pragma omp master
//        {
//            std::cout << "Using OpenMP with " << omp_get_num_threads() << " threads" << std::endl;
//        }
//
//        long ldsSeedWorker = ldsSeed + omp_get_thread_num() * 2 * numberBatches;
//        long prsSeedWorker = prsSeed + omp_get_thread_num() * 2 * numberPaths;
//        
//        RandomizedLDS<SobolRsg, RandomSequenceGenerator<MersenneTwisterUniformRng> > rldsg(3, ldsSeedWorker, prsSeedWorker);        
//
//        #pragma omp for 
//        for (long k=0; k < numberBatches; k++) {
//            GeneralStatistics innerStats;
//            vector<double> p;
//            double payoff;
//            
//            for (long i=0; i < numberPaths; i++) {
//                p = rldsg.nextSequence().value;
//                payoff = quantileHestonProcess(p[0], p[1], p[2], r, rho, kappa, theta, sigma, u, vu, Su, t) - K;
//                payoff = (payoff>0 ? payoff : 0)*exp(-r*(t-u));
//
//                innerStats.add(payoff);
//            }
//        
//            #pragma omp critical
//            outerStats.add(innerStats.mean());
//            
//            rldsg.nextRandomizer();
//        }
//    }
//    
//    #else
//    cout << "No OpenMP!" << endl;
//    #endif
//    
//    double mean = outerStats.mean();
//    double sds = outerStats.errorEstimate();
//    
//    std::cout << "mean: " << mean << " errorEstimate: " << sds << " paths: " << numberPaths << " batches: " << numberBatches << std::endl;
//    
//    return mean;
//}

//double HestonCallQMCExactCond(double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t, double K, long numberPaths, long numberBatches, long ldsSeed=2, long prsSeed=3) {
//    
//    using namespace QuantLib;
//    
//    GeneralStatistics outerStats;
//    
//    #ifdef _OPENMP
//
//    #pragma omp parallel shared(outerStats)
//    {
//        #pragma omp master
//        {
//            std::cout << "Using OpenMP with " << omp_get_num_threads() << " threads" << std::endl;
//            std::cout << "QMC with Cond" << std::endl;
//
//        }
//
//        long ldsSeedWorker = ldsSeed + omp_get_thread_num() * 2 * numberBatches;
//        long prsSeedWorker = prsSeed + omp_get_thread_num() * 2 * numberPaths;
//        
//        RandomizedLDS<SobolRsg, RandomSequenceGenerator<MersenneTwisterUniformRng> > rldsg(2, ldsSeedWorker, prsSeedWorker);        
//
//        #pragma omp for 
//        for (long k=0; k < numberBatches; k++) {
//            GeneralStatistics innerStats;
//            vector<double> p;
//            double payoff;
//            
//            for (long i=0; i < numberPaths; i++) {
//                p = rldsg.nextSequence().value;
//                
//                double vt = quantileSquareRootProcess(p[0], kappa, theta, sigma, u, vu, t);
//                double intv = quantileIntSquareRootProcess(p[1], kappa, theta, sigma, u, vu, t, vt);
//                double intsqrvdW = (vt - vu - kappa*theta*(t-u) + kappa * intv) / sigma;
//                
//                double xi = exp(-0.5*rho*rho*intv + rho*intsqrvdW);
//                double sigmatilde = sqrt(intv/(t-u));
//                                
//                payoff = bscall(Su*xi, K, u, t, r, sigmatilde*sqrt(1-rho*rho));
//
//                innerStats.add(payoff);
//            }
//        
//            #pragma omp critical
//            outerStats.add(innerStats.mean());
//            
//            rldsg.nextRandomizer();
//        }
//    }
//    
//    #else
//    cout << "No OpenMP!" << endl;
//    #endif
//    
//    double mean = outerStats.mean();
//    double sds = outerStats.errorEstimate();
//    
//    std::cout << "mean: " << mean << " errorEstimate: " << sds << " paths: " << numberPaths << " batches: " << numberBatches << std::endl;
//    
//    return mean;
//}
