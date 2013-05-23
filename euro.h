// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#ifndef volmodels_euro_h
#define volmodels_euro_h

#include "utils.h"

Result HestonCallMCExact(double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t, double Strike, long NumberOfPaths, long prsSeed);

Result HestonCallQMCExact(double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t, double K, long numberPaths, long numberBatches, long prsSeed);

Result HestonCallCondQMCExact(double r, double rho, double kappa, double theta, double sigma, double u, double vu, double Su, double t, double K, long numberPaths, long numberBatches, long prsSeed);

#endif
