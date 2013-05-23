// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#ifndef volmodels_barrier_h
#define volmodels_barrier_h

double QMCBarrierDownOutCallNaive(double r, double rho, double kappa, double theta, double sigma, double v0, double S0, double T, double K, double H, long numberPaths, long numberBatches, long prsSeed=123, long numberTimeMonitors = 32);

#endif
