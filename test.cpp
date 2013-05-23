// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <stdio.h>
#include <iostream>
#include "asian.h"

using namespace std;

//void testPoissonPath()
//{
//    double lambda = 0.11;
//    
//    int tsteps = 4;
//    int m = log2(tsteps);
//    double dt = 1.0/tsteps;
//    
//    vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//    
//    vector<double> qs(tsteps, 0.999);
//    
//    PIO(times);
//    PIO(qs);
//
//    vector<double> v = poissonPath(times, lambda, qs);
//    
//    PIO(v);   
//}
//
//void testPoissonPathBridge()
//{
//    double lambda = 0.11;
//    
//    int tsteps = 4;
//    int m = log2(tsteps);
//    double dt = 1.0/tsteps;
//    
//    vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//    
//    vector<double> qs(tsteps, 0.999);
//    
//    vector<double> v = poissonPathBridge(times, lambda, qs);
//    
//    PIO(v);   
//}
//
//void testSumPoissonPath()
//{
//    double lambda = 0.11;
//    
//    int tsteps = 4;
//    int m = log2(tsteps);
//    double dt = 1.0/tsteps;
//    
//    vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//    
//    vector<double> qs(tsteps, 0.999);
//
//    vector<double> v = poissonPath(times, lambda, qs);
//    
//    double mubar = -0.12;
//    double sigma = 0.15;
//    double mu_s = std::log(1+mubar)-0.5*(sigma*sigma);
//    vector<double> qs2(tsteps+1, 0.999);
//    
//    vector<double> S = sumPoissonPath(v, mu_s, sigma, qs2);
//    
//    PIO(S);
//}
//
//void testSumPoissonPathBridge()
//{
//    double lambda = 0.11;
//    
//    int tsteps = 4;
//    int m = log2(tsteps);
//    double dt = 1.0/tsteps;
//    
//    vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//    
//    vector<double> qs(tsteps, 0.999);
//    
//    vector<double> v = poissonPath(times, lambda, qs);
//    
//    double mubar = -0.12;
//    double sigma = 0.15;
//    double mu_s = std::log(1+mubar)-0.5*(sigma*sigma);
//    vector<double> qs2(tsteps+1, 0.999);
//    
//    vector<double> S = sumPoissonPathBridge(v, mu_s, sigma, qs2);
//    
//    PIO(S);
//}
//
//void testIntSquareRootProcess()
//{
//    using namespace boost::numeric::ublas;
//
//    double kappa = 3.99;
//    double theta = 0.014;
//    double sigma = 0.27;
//    double u = 0.;
//    double vu = 0.008836;
//    double r = 0.0319;
//    double Su = 100.;
//    double K = 100.;
//    double rho = -0.79;
//    
//    double lambda = 0.11;
//    double mubar = -0.12;
//    double sigma_s = 0.15;
//    double mu_s=log(1+mubar)-0.5*sigma_s*sigma_s;
//
//    long numberPaths = 1;
//	long numberBatches = 1;
//
//    int tsteps = 4;
//    int m = log2(tsteps);
//    double dt = 5.0/tsteps;
//    
//    std::vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//    
//    PIO(times);
//    
//    double v0 = vu;
//    
//    tsteps = times.size();
//    int tm = log2(tsteps-1);
//   	int sact = 1+3*(2^(tm)-1)+4*(tsteps-1);
//    int Kscramble=31;
//                  
//    matrix<double> P(sact, numberPaths); 
//    
//    for (int i = 0; i < P.size1(); ++i)
//        for (int j = 0; j < P.size2(); ++j)
//            P(i,j) = 0.5;
//    
//    int i = 0;
//    
//    matrix_column<matrix<double>> Pc(P,i);
//    
//    // generate the value of the square root process at the time points
//
//    std::vector<double> qs(Pc.begin(), Pc.begin()+3*(tsteps-1));
//    std::vector<double> Vs = squareRootBridge(times, kappa, theta, sigma, v0, qs);
//    
//    PIO(Vs);
//    
//    // generate the integrals of V and the stochastic integrals of sqrt(V)
//    
//    qs.assign(Pc.begin()+3*(tsteps-2), Pc.begin()+3*(tsteps-2)+tsteps);
//    std::vector<double> intVs(tsteps, 0.0);
//    std::vector<double> intsqrVdWs(tsteps, 0.0);
//    
//    for (int j=1; j < tsteps; ++j) {
//        intVs[j] = quantileIntSquareRootProcess(qs[j-1], kappa, theta, sigma, times[j-1], Vs[j-1], times[j], Vs[j]);
//        intsqrVdWs[j] = (Vs[j] - Vs[j-1] - kappa*theta*(times[j]-times[j-1]) + kappa * intVs[j]) / sigma;                    
//    }
//    
//    PIO(intVs);
//    PIO(intsqrVdWs);
//
//}
//
//void testqbessel()
//{
//    double nu = 0.5325;
//    double z = 1.0093e-04;
//    double b = qbessel(0.5, nu, z);
//    PO(b);
//}
//
//void testqgamma()
//{
//    double a = 1.53251;
//    double b = 3.52428e-05;
//    double G = qgamma(0.5, a, 1/b);
//    PO(G);
//}
//
//int testMCAsian()
//{
//    double kappa = 6.21;
//    double theta = 0.019;
//    double sigma = 0.61;
//    double u = 0;
//    double vu = 0.010201;
//    double t = 1.0;
//    double r = 0.0319;
//    double Su = 100;
//    double K = 100;
//    double rho = -0.7;
//    
//    int tsteps = 8;
//    int m = log2(tsteps);
//    double dt = 1.0/tsteps;
//    
//    vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//    
//    PIO(times);
//    
//    long seed = 1;    
//    int maxPowerTwo = 1;
//    
//    for (int powerTwo = 1; powerTwo <= maxPowerTwo; ++powerTwo)
//    {
//        long numberPaths = 1 << powerTwo;
//        int numberBatches = 2;
//        
//        Result mcResult = MCAsian(r, rho, kappa, theta, sigma, vu, Su, times, K, numberPaths, numberBatches, seed);
//        
//        std::cout << mcResult.mean << " " << mcResult.stderr << endl;
//    }
//}
//
//int testQMCAsian()
//{
//    double kappa = 6.21;
//    double theta = 0.019;
//    double sigma = 0.61;
//    double u = 0;
//    double vu = 0.010201;
//    double t = 1.0;
//    double r = 0.0319;
//    double Su = 100;
//    double K = 100;
//    double rho = -0.7;
//    
//    int tsteps = 8;
//    int m = log2(tsteps);
//    double dt = 1.0/tsteps;
//    
//    vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//    
//    long seed = 1;    
//    int maxPowerTwo = 1;
//    
//    for (int powerTwo = 1; powerTwo <= maxPowerTwo; ++powerTwo)
//    {
//        long numberPaths = 1 << powerTwo;
//        int numberBatches = 2;
//        
//        Result mcResult = QMCAsian(r, rho, kappa, theta, sigma, vu, Su, times, K, numberPaths, numberBatches, seed);
//        
//        std::cout << mcResult.mean << " " << mcResult.stderr << endl;
//    }
//}
//

//void testsquareRootPath()
//{
//    double kappa=6.21;
//    double theta=0.019;
//    double sigma=0.61;
//    double v0=0.010201;
//
//    int tsteps = 8;
//    double dt = 1.0/tsteps;
//    
//    vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//
//    PIO(times);
//    
//    vector<double> qs(tsteps, 0.9);
//        
//    vector<double> path = squareRootPath(times, kappa, theta, sigma, v0, qs);
//    
//    PIO(path);
//}
//
//void testQuantileIntSqrProcess()
//{
//    double S0=100;
//    double K=100;
//    double V0=0.008836;
//
//    double kappa=3.99;
//    double theta=0.014;    
//    double sigma=0.27;
//    
//    double rho=-0.79;
//    double lambda=0.11;
//    double mubar=-0.12;
//    double sigma_s=0.15;
//    double mu_s=log(1+mubar)-0.5*(sigma_s*sigma_s);
//    double r=0.0319;
//    
//    double intVs = quantileIntSquareRootProcess(0.9, kappa, theta, sigma, 0.0, 0.008836, 1.25, 0.0289481);
//    
//    cout << setprecision(10) << intVs << endl;
//}
//
//
//
//int testMCAsianSVJ() 
//{
//    double S0=100;
//    double K=100;
//    double V0=0.008836;
//    double kappa=3.99;
//    double theta=0.014;
//    double sigma=0.27;
//    double rho=-0.79;
//    double lambda=0.11;
//    double mubar=-0.12;
//    double sigma_s=0.15;
//    double mu_s=log(1+mubar)-0.5*(sigma_s*sigma_s);
//    double r=0.0319;
//
//    int tsteps = 16;
//    int m = log2(tsteps);
//    double dt = 5.0/tsteps;
//    
//    vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//    
//    PIO(times);
//    
//    // parse arguments
//    
//    long seed = 1234;    
//    int maxPowerTwo = 1;
//    
//    Result mcResult = MCAsianSVJ(lambda, mu_s, sigma_s, r, rho, kappa, theta, sigma, V0, S0, times, K, 2, 30, seed);
//}
//
//
//
//int testQMCAsian() 
//{
//    double S0=100;
//    double K=100;
//    double V0=0.008836;
//    double kappa=3.99;
//    double theta=0.014;
//    double sigma=0.27;
//    double rho=-0.79;
//    double lambda=0.11;
//    double mubar=-0.12;
//    double sigma_s=0.15;
//    double mu_s=log(1+mubar)-0.5*(sigma_s*sigma_s);
//    double r=0.0319;
//    
//    int tsteps = 4;
//    int m = log2(tsteps);
//    double dt = 1.0/tsteps;
//    
//    vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//    
//    PIO(times);
//    
//    // parse arguments
//    
//    long seed = 4711;    
//    int maxPowerTwo = 1;
//    
//    Result mcResult = QMCAsian2(r, rho, kappa, theta, sigma, V0, S0, times, K, 2, 2, seed);
//    
//    PO(mcResult.mean);
//    PO(mcResult.stderr);
//}
//
//#include <iostream>
//#include <iomanip>
//
//#include "exact.h"
//#include "asian.h"
//#include "bridges.h"
//
//using namespace std;
//
int main(int argc, char *argv[]) 
{
    double kappa = 6.21;
    double theta = 0.019;
    double sigma = 0.61;
    double u = 0;
    double vu = 0.010201;
    double t = 5.0;
    double r = 0.0319;
    double Su = 100;
    double K = 100;
    double rho = -0.7;
        
    int tsteps = 16;
    int m = log2(tsteps);
    double dt = 1.0/tsteps;
    
    vector<double> times;
    for (int i = 0; i <= tsteps; i++)
        times.push_back(i*dt);
    
    // parse arguments
    
    int seed = argc>1 ? atol(argv[1]) : 4711;    
    int maxPowerTwo = argc>2 ? atoi(argv[2]) : 10;
    
    EPS    = argc>3 ? atof(argv[3]) : 1E-12;
    STDEVS = argc>4 ? atoi(argv[4]) : 12;
    DIGITS = argc>5 ? atoi(argv[5]) : 12;
    
    int numberPaths = 8;
    int numberBatches = 8;
    
    Result res = QMCAsianBridge(r, rho, kappa, theta, sigma, vu, Su, times, K, numberPaths, numberBatches, seed);
        
    printf("%.6f (%.6f)\n", res.mean, res.stderr);
}

//int main(int argc, char *argv[]) {
//    double val = qpois(0.45789, 1.086237055);
//    printf("%f\n", val);
//}

//int main(int argc, char *argv[])
//{
//    double kappa=6.21;
//    double theta=0.019;
//    double sigma=0.61;
//    double v0=0.010201;
//
//    int tsteps = 4;
//    double dt = 1.0/tsteps;
//    
//    vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//
//    vector<double> qs = {0.982065, 0.490717, 0.679456, 0.908689, 0.450825};
//    vector<double> Vs = {0.010201, 0.0215799, 0.00642798, 0.0020178, 0.00804512};
//
//    PIO(qs);
//    PIO(times);
//    PIO(Vs);
//            
//    std::vector<double> intVs(tsteps, 0.0);
//    std::vector<double> intsqrVdWs(tsteps, 0.0);
//    int qindex = 0;
//
//    qindex = 0;
//    for (int j=1; j < tsteps; ++j) {
//        intVs[j] = quantileIntSquareRootProcess(qs[qindex], kappa, theta, sigma, times[j-1], Vs[j-1], times[j], Vs[j]);
//        intsqrVdWs[j] = (Vs[j] - Vs[j-1] - kappa*theta*(times[j]-times[j-1]) + kappa * intVs[j]) / sigma;                    
//        qindex++;
//    }
//    
//    PIO(intVs); PIO(intsqrVdWs);
//    
//    qindex = 0;
//    for (int j = tsteps-1; j > 1; --j) {
//        intVs[j] = quantileIntSquareRootProcess(qs[qindex], kappa, theta, sigma, times[j-1], Vs[j-1], times[j], Vs[j]);
//        intsqrVdWs[j] = (Vs[j] - Vs[j-1] - kappa*theta*(times[j]-times[j-1]) + kappa * intVs[j]) / sigma;
//        qindex++;
//    }
//
//    PIO(intVs); PIO(intsqrVdWs);
//}


//#include <stdio.h>
//
//#include "exact.h"
//#include "asianSVJ.h"
//
//using namespace std;
//
//int main(int argc, char *argv[]) 
//{
////    double kappa = 6.21;
////    double theta = 0.019;
////    double sigma = 0.61;
////    double u = 0;
////    double vu = 0.010201;
////    double t = 1.0;
////    double r = 0.0319;
////    double Su = 100;
////    double K = 100;
////    double rho = -0.7;
////    
////    double lambda = 0.11;
////	double mu_s = -0.1391;
////	double sigma_s = 0.15;
//    
//    double S0=100;
//    double K=100;
//    double V0=0.008836;
//    double kappa=3.99;
//    double theta=0.014;
//    double sigma=0.27;
//    double rho=-0.79;
//    double lambda=0.11;
//    double mubar=-0.12;
//    double sigmas=0.15;
//    double mus=log(1+mubar)-0.5*(sigmas*sigmas);
//    double r=0.0319;
//    
//    int tsteps = 4;
//    int m = log2(tsteps);
//    double dt = 5.0/tsteps;
//    
//    vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//    
//    // parse arguments
//    
//    long seed = argc>1 ? atol(argv[1]) : 1234;    
//    int maxPowerTwo = argc>2 ? atoi(argv[2]) : 8;
//    
//    EPS    = argc>3 ? atof(argv[3]) : 1E-12;
//    STDEVS = argc>4 ? atoi(argv[4]) : 12;
//    DIGITS = argc>5 ? atoi(argv[5]) : 12;
//    
//    long numberPaths = 8;
//    int numberBatches = 30;
//            
//    Result res = QMCAsianBridgeSVJ(lambda, mus, sigmas, r, rho, kappa, theta, sigma, V0, S0, times, K, numberPaths, numberBatches, seed);
//}

//int main(int argc, char *argv[])
//{
//    double lambda = 0.11;
//    
//    int tsteps = 4;
//    int m = log2(tsteps);
//    double dt = 5.0/tsteps;
//    
//    vector<double> times;
//    for (int i = 0; i <= tsteps; i++)
//        times.push_back(i*dt);
//    
//    double B = qbinom(0.42254, 2, 0.5);
//    PO(B);
//    
//    vector<double> qs = {0.95425, 0.42254, 0.82708, 0.70041};
//    
//    vector<double> v = poissonPathBridge(times, lambda, qs);
//    
//    PO(lambda);
//    PIO(times);
//    PIO(qs);
//    PIO(v);   
//}
