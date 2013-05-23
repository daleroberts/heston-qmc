// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <iostream>
#include <iomanip>
#include "exact.h"
#include "utils.h"
#include "euro.h"
#include "barrier.h"
#include "asian.h"
#include "ql_imps.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

#define PO(e) std::cout << #e << ": " << e << std::endl
#define PIO(v) std::cout << #v << ": "; copy(v.begin(), v.end(), ostream_iterator<double>(cout, " ")); std::cout << " (size: " << v.size() << ")" << std::endl

using namespace std;

//void testHestonCallMCExact() 
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
//    long NumberOfPaths = 32 * 1024;
//
////    cout << "Enter Number of Paths" << endl;
////    cin >> NumberOfPaths;
//
//    #ifdef _OPENMP
//    std::cout << "Using OpenMP" << std::endl;
//    double initialTime = omp_get_wtime();
//    #else
//    clock_t initialTime = clock();
//    #endif
//
//    HestonCallMCExact(r, rho, kappa, theta, sigma, u, vu, Su, t, K,  NumberOfPaths);
//    
//    #ifdef _OPENMP
//    double timeSpent = omp_get_wtime() - initialTime;
//    #else 
//    double timeSpent = difftime(clock(), initialTime) / CLOCKS_PER_SEC;
//    #endif
//    
//    cout << "time: " << timeSpent << "s" << endl;
//
////    cin >> kappa;
//}

void testHestonCallQMCExact(long numberPaths = 512, long numberBatches = 32) 
{
    double kappa = 6.21;
    double theta = 0.019;
    double sigma = 0.61;
    double r = 0.0319;
    double rho = -0.7;
    double u = 0;
    double vu = 0.010201;
    double t = 1.0;
    double Su = 100;
    double K = 100;

    #ifdef _OPENMP
    double initialTime = omp_get_wtime();
    #else
    clock_t initialTime = clock();
    #endif
    
    HestonCallQMCExact3(r, rho, kappa, theta, sigma, u, vu, Su, t, K, numberPaths, numberBatches, 1234UL);
    
    #ifdef _OPENMP
    double timeSpent = omp_get_wtime() - initialTime;
    #else 
    double timeSpent = difftime(clock(), initialTime) / CLOCKS_PER_SEC;
    #endif
    
//    cout << " Time: " << timeSpent << "s" << endl;
    cout << setw(12) << setprecision(2) << timeSpent << endl;
}

void testHestonFdBarrier()
{
    double kappa = 6.21;
    double theta = 0.019;
    double sigma = 0.61;
    double r = 0.0319;
    double rho = -0.7;
    double u = 0;
    double vu = 0.010201;
    double t = 1.0;
    double Su = 100;
    double K = 100;
    double H = 95.;

    double npv = HestonFdBarrierDownOutCall(r, rho, kappa, theta, sigma, vu, Su, t-u, K, H); 

    std::cout << npv << std::endl;
}

void testHestonQMCBarrier()
{
    double kappa = 6.21;
    double theta = 0.019;
    double sigma = 0.61;
    double r = 0.0319;
    double rho = -0.7;
    double u = 0;
    double vu = 0.010201;
    double t = 1.0;
    double Su = 100;
    double K = 100;
    double H = 80;

    double npv = QMCBarrierDownOutCallNaive(r, rho, kappa, theta, sigma, vu, Su, t-u, K, H, 1024, 4, 123, 200); 

    std::cout << npv << std::endl;
}

void testHestonAsian() 
{
    double kappa = 6.21;
    double theta = 0.019;
    double sigma = 0.61;
    double u = 0.;
    double vu = 0.010201;
    double r = 0.0319;
    double Su = 100.;
    double K = 100.;
    double rho = -0.7;

    long numberPaths = 2;
	long numberBatches = 4;

	vector<double> t;

	t.push_back(0.25);
	t.push_back(0.50);
	t.push_back(0.75);
	t.push_back(1.00);

    cout << "Asian Option " 
		 << MCAsianNaive2(r, rho, kappa, theta, sigma, u, vu, Su, K, t, numberPaths, numberBatches) << endl;
}

void testQMCAsianNaive() 
{
    double kappa = 6.21;
    double theta = 0.019;
    double sigma = 0.61;
    double u = 0.;
    double vu = 0.010201;
    double r = 0.0319;
    double Su = 100.;
    double K = 100.;
    double rho = -0.7;

    long numberPaths = 256;
	long numberBatches = 4;

	double t = 1.0;
	
	//double QMCAsianNaive(double r, double rho, double kappa, double theta, double sigma, double v0, double S0, double T, double K, long numberPaths, long numberBatches, long prsSeed, long numberTimeMonitors) 
    cout << "Asian Option " 
		 << QMCAsianNaive(r, rho, kappa, theta, sigma, vu, Su, (t-u), K, numberPaths, numberBatches, 1L, 4) << endl;
}

void checkOpenMP()
{
    #ifdef _OPENMP

    #pragma omp parallel 
    {
       #pragma omp master
       std::cout << "OpenMP enabled. Threads: " << omp_get_num_threads() << std::endl;
	}

    #else

    std::cout << "OpenMP disabled." << std::endl;
    
    #endif            
}

void testStockPriceBridge()
{
    int tsteps = 4;
    double dt = 1.0/tsteps;
    
    vector<double> times;
    for (int i = 0; i <= tsteps; i++)
        times.push_back(i*dt);
    
//    PIO(times);
    
    vector<double> qs(tsteps, 0.5);
    
//    PIO(qs);

    double intsim1 = 0.5;
    double intVsdWs = 0.5;
    double r=0.0;
    double rho = -0.5;

    double S0 = 100;

    vector<double> drifts = {log(S0)};
    for (int i = 1; i <= tsteps; i++) {
        drifts.push_back(drifts[i-1] + r * (times[i] - times[i-1]) - 0.5 * intsim1 + rho * intVsdWs);
    }
    
//    PIO(drifts);
    
    vector<double> vols = {0.0};
    for (int i=1; i <= tsteps; i++) {
        vols.push_back(vols[i-1]+(1-rho*rho)*intsim1);
    }
    
//    PIO(vols);
            
    vector<double> Ss = stockPriceBridge(times, drifts, vols, qs);
    PIO(Ss);
}

void testSquareRootBridge()
{
    double S0=100;
    double K=100;
    double V0=0.010201;
    double kappa=6.21;
    double theta=0.019;
    double sigma=0.61;
    double rho=-0.7;
    double r=0.0319;
    
    int tsteps = 4;
    int m = log2(tsteps);
    double dt = 1.0/tsteps;
    
    vector<double> times;
    for (int i = 0; i <= tsteps; i++)
        times.push_back(i*dt);
    
    vector<double> qs(1+3*(tsteps-2), 0.5);

    vector<double> v = squareRootBridge(times, kappa, theta, sigma, V0, qs);
    
    PIO(v);
}

void testQMCAsianBridge() 
{
    double kappa = 6.21;
    double theta = 0.019;
    double sigma = 0.61;
    double u = 0.;
    double vu = 0.010201;
    double r = 0.0319;
    double Su = 100.;
    double K = 100.;
    double rho = -0.7;

    long numberPaths = 2;
	long numberBatches = 2;
    
    int tsteps = 4;
    int m = log2(tsteps);
    double dt = 1.0/tsteps;
    
    vector<double> times;
    for (int i = 0; i <= tsteps; i++)
        times.push_back(i*dt);
    
    cout << "Asian Option " 
		 << QMCAsianBridge(r, rho, kappa, theta, sigma, vu, Su, times, K, numberPaths, numberBatches, 1L) << endl;
}
                        
int main(int argc, char *argv[]) {

    // parse arguments
    EPS    = argc>1 ? atof(argv[1]) : 1E-9;
    STDEVS = argc>2 ? atoi(argv[2]) : 9;
    DIGITS = argc>3 ? atoi(argv[3]) : 9;
    
    long numberPaths = argc>4 ? atol(argv[4]) : 512;
    long numberBatches = argc>5 ? atol(argv[5]) : 4;
    
//    cout << setw(8) << scientific << setprecision(0) << EPS << setw(8) << STDEVS << setw(8) << DIGITS << setw(8) << numberPaths << setw(8) << numberBatches;
//    testHestonCallQMCExact(numberPaths, numberBatches);
//    testHestonFdBarrier();
//    testHestonQMCBarrier();
//    testHestonAsian();
	checkOpenMP();
//	testQMCAsianNaive();
//    testStockPriceBridge();
//    testSquareRootBridge();
//    testQMCAsianBridge();
//    testqgamma();

    #if defined _WIN32 || defined _WIN64
	char k;
	std::cin >> k;
    #endif
}
