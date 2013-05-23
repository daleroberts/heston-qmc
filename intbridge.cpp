// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <stdio.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "sobol.h"
#include "bridges.h"
#include "stats.h"
#include "exact.h"

int main(int argc, char *argv[]) {
    int number_batches = 30;
    int number_paths = 128;
    int seed = 31337;
    
    
    double kappa = 6.21;
    double theta = 0.019;
    double sigma = 0.61;
    double v0 = 0.010201;
    
    int tdim = 4;
    double dt = 1.0/tdim;
    
    vector<double> times;
    for (int i = 0; i <= tdim; i++)
        times.push_back(i*dt);
    
    int tsteps = times.size();
    
    int sact = 3*(tdim-1)+tsteps+1;

    using boost::numeric::ublas::matrix;
    using boost::numeric::ublas::matrix_row;
    using boost::numeric::ublas::column_major;

    vector<vector<double>> outerstats(tdim);

    #pragma omp parallel shared(outerstats)
    {
        #pragma omp for
        for (int batch_id = 0; batch_id < number_batches; ++batch_id) 
        {
            vector<vector<double>> innerstats(tdim);
            for (int i = 0; i < tdim; ++i)
                innerstats.reserve(number_paths);
            
            matrix<double, column_major> P(number_paths, sact);         
            
            #pragma omp critical
            ScrambleSobol(seed + batch_id*number_batches, P);
            
            for (int path_id = 0; path_id < number_paths; ++path_id) 
            {
                matrix_row<matrix<double, column_major>> Pc(P,path_id);
                
#if 0
                std::vector<double> qs(Pc.begin(), Pc.begin()+3*(tdim-1)+1);
                std::vector<double> Vs = squareRootBridge(times, kappa, theta, sigma, v0, qs);
                
                qs.assign(Pc.begin()+3*(tdim-1)+1, Pc.begin()+3*(tdim-1)+tsteps+1);
                std::vector<double> intVs(tsteps, 0.0);

                int k=0;
                for (int j = 1; j < tsteps; ++j) {
                    intVs[j] = quantileIntSquareRootProcess(qs.at(k), kappa, theta, sigma, times.at(j-1), Vs.at(j-1), times.at(j), Vs.at(j));
                    k++;
                }
#endif
                
#if 0 // diff
                std::vector<double> qs(Pc.begin(), Pc.begin()+3*(tdim-1)+1);
                std::vector<double> Vs = squareRootBridge(times, kappa, theta, sigma, v0, qs);
                
                qs.assign(Pc.begin()+3*(tdim-1)+1, Pc.begin()+3*(tdim-1)+tsteps+1);
                std::vector<double> z(tsteps, 0.0);
                
                int k=0;
                for (int j = 1; j < tsteps; ++j) {
                    z[j] = quantileIntSquareRootProcess(qs.at(k), kappa, theta, sigma, times[0], Vs[0], times.at(j), Vs.at(j));
                    k++;
                }
                
                std::vector<double> intVs(tsteps, 0.0);
                for (int i = tdim; i > 0; --i) {
                    intVs[i] = z[i] - z[i-1];
                }
#endif                
                
                
#if 0                
                std::vector<double> qs(Pc.begin(), Pc.begin()+3*(tdim-1)+1);
                std::vector<double> Vs = squareRootBridge(times, kappa, theta, sigma, v0, qs);
                
                qs.assign(Pc.begin()+3*(tdim-1)+1, Pc.begin()+3*(tdim-1)+tsteps+1);
                std::vector<double> intVs(tsteps, 0.0);

                // try 2: backwards
                int k=0;
                for (int j = tdim; j > 0; --j) {
                    intVs[j] = quantileIntSquareRootProcess(qs.at(k), kappa, theta, sigma, times.at(j-1), Vs.at(j-1), times.at(j), Vs.at(j));
                    k++;
                }
#endif                

#if 1
                std::vector<double> qs(Pc.begin(), Pc.end());

                int m = log2(times.size()-1);
                int h = 1 << m; // 2^m
                
                // the modified dates
                vector<double> s(h+1, 0.0);
                for (int i = 1; i < h+1; ++i) {
                    s[i] = sigma*sigma/(4*kappa) * (exp(kappa*times[i])-1);
                }
                
                // the dimension of the square-root process
                double delta = 4*kappa*theta/(sigma*sigma); 
                
                int jmax = 1;
                int iq = 0; // index of quantiles qs
                
                vector<double> x(h+1, 0.0); // squared Bessel bridge
                vector<double> v(h+1, v0);  // square-root process
                vector<double> z(h+1, 0.0); // int_0^{t_i} V ds
                
                x[0] = v0;
                x[h] = s[h] * qchisq(qs[iq], delta, x[0]/s[h]);
                
                iq++;
                
                v[h] = x[h] * exp(-kappa * times[h]);
                z[h] = quantileIntSquareRootProcess(qs[iq], kappa, theta, sigma, times[0], x[0], times[h], v[h]);

                iq++;
                
                for (int k = 1; k <= m; ++k) {
                    int imin = h/2, i = imin, l = 0, r = h; 
                    for (int j = 1; j <= jmax; ++j) {
                        double Z, P, lambda, a, b, G;
                        lambda=(1/(2*(s[r]-s[l])))*( (s[r]-s[i])*x[l]/(s[i]-s[l]) +(s[i]-s[l])*x[r]/(s[r]-s[i]));
                        P = qpois(qs[iq], lambda);
                        
                        iq++;
                        
                        Z = qbessel(qs[iq], delta/2-1, sqrt(x[l]*x[r])/(s[r]-s[l]));
                        
                        iq++;
                        
                        a = P + 2*Z + delta/2;
                        b = (s[r]-s[l])/(2*(s[i]-s[l])*(s[r]-s[i]));
                        G = qgamma(qs[iq], a, 1/b);
                        
                        iq++;

                        x[i] = G;
                        v[i] = x[i] * exp(-kappa * times[i]);
                        z[i] = quantileIntSquareRootProcess(qs[iq], kappa, theta, sigma, times[0], v[0], times[i], v[i]);

                        iq++;
                        
                        i = i+h;
                        l = l+h;
                        r = r+h;
                        
                    }
                    jmax = jmax*2;
                    h = imin;
                }
                
                std::vector<double> intVs(tsteps, 0.0);
                for (int i = tdim; i > 0; --i) {
                    intVs[i] = z[i] - z[i-1];
                }
#endif            

            

                for (int i = 0; i < tdim; ++i) {
                    innerstats[i].push_back(intVs[i+1]);
                }
            }            
            
            vector<double> means(tdim);
            for (int i = 0; i < tdim; ++i) {
                means[i] = Mean(innerstats[i]);
            }
          
            #pragma omp critical
            for (int i = 0; i < tdim; ++i) {
                outerstats[i].push_back(means[i]);
            }
        }
    }
    
    for (int i = 0; i < tdim; ++i) {
        printf("%.8f (%.8f)  ", Mean(outerstats[i]), StdError(outerstats[i]));
    }
    
    printf("\n");
}
