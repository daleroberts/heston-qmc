// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#ifndef volmodels_bridges_h
#define volmodels_bridges_h

#include "utils.h"

using namespace std;

vector<double> stockPricePath(const vector<double>& times, const vector<double>& drifts, const vector<double>& vols, const vector<double>& qs)
{
    int tdim = times.size();
    
    vector<double> s(tdim, 0);
    
    s[0] = exp(drifts[0]);
    for (int i = 1; i < tdim; ++i) {
        s[i] = exp(drifts[i] + sqrt(vols[i])*qnorm(qs[i-1]));
        
        #ifdef DEBUG
        qs.at(i-1);
        #endif
    }
    
    return s;
}

vector<double> stockPriceBridge(const vector<double>& times, const vector<double>& drifts, const vector<double>& vols, const vector<double>& qs)
{
    int tdim = times.size()-1;
    int m = log2(tdim);
    
    vector<double> s(times.size(), 0);
    vector<double> Z(tdim, 0);
    
    for (int i = 0; i < Z.size(); ++i) {
        Z[i] = qnorm(qs[i]);

        #ifdef DEBUG
        qs.at(i);
        #endif
    }

    int h = 1 << m; // 2^m
    int jmax = 1;
        
    s[h] = exp(drifts[h] + sqrt(vols[h])*Z[0]);
    s[0] = exp(drifts[0]);
    
    int qindex = 1;
    
    for (int k = 1; k <= m; ++k) {
        int l = 0, r = h, imin = h/2, i = imin;
        double a, b;
        for (int j = 1; j <= jmax; ++j) {
            a = drifts[i]-drifts[l]+log(s[l])+(vols[i]-vols[l])/(vols[r]-vols[l])*(log(s[r])-log(s[l])+drifts[l]-drifts[r]); 
            b = vols[i]-vols[l]-((vols[i]-vols[l])*(vols[i]-vols[l]))/(vols[r]-vols[l]);
            s[i] = exp(a + sqrt(b) * Z[qindex]);
            
            i = i + h;
            l = l + h;
            r = r + h;
            qindex++;
        }
        jmax = 2 * jmax;
        h = imin;
    }
        
    return s;
}

vector<double> squareRootPath(const vector<double>& times, double kappa, double theta, double sigma, double v0, const vector<double>& qs)
{
    int m = log2(times.size()-1);
    int h = 1 << m; // 2^m
    
    // the modified dates
    vector<double> s(h+1, 0.0);
    for (int i = 1; i < h+1; ++i) {
        s[i] = sigma*sigma/(4*kappa) * (exp(kappa*times[i])-1);
    }
    
    // the dimension of the square-root process
    double delta = 4*kappa*theta/(sigma*sigma); 
    
    // contains the squared Bessel bridge, which is finally modified to yield the square-root process 
    vector<double> x(h+1, v0);
    for (int i = 1; i < h+1; ++i) {
        x[i] = (s[i]-s[i-1])*qchisq(qs[i-1],delta,x[i-1]/(s[i]-s[i-1]));

        #ifdef DEBUG
        qs.at(i-1);
        #endif
    }
    
    vector<double> v(times.size(), v0);
    for (int i = 1; i < times.size(); ++i) {
        v[i] = x[i] * exp(-kappa * times[i]);
    }
    
    return v;
}

vector<double> squareRootBridge(const vector<double>& times, double kappa, double theta, double sigma, double v0, const vector<double>& qs)
{
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
        
    // contains the squared Bessel bridge, which is finally modified to yield the square-root process 
    vector<double> x(h+1, 0);
    x[0] = v0;
    x[h] = s[h] * qchisq(qs[iq], delta, x[0]/s[h]);
    
    #ifdef DEBUG
    qs.at(iq);
    #endif

    iq++;
        
    for (int k = 1; k <= m; ++k) {
        int imin = h/2, i = imin, l = 0, r = h; 
        for (int j = 1; j <= jmax; ++j) {
            double Z, P, lambda, a, b, G;
            lambda=(1/(2*(s[r]-s[l])))*( (s[r]-s[i])*x[l]/(s[i]-s[l]) +(s[i]-s[l])*x[r]/(s[r]-s[i]));
            P = qpois(qs[iq], lambda);
            
//            vector<double> vals = {qs[iq],lambda};
//            PIO(vals);
//            PO(P);
            
            #ifdef DEBUG
            qs.at(iq);
            #endif

            iq++;
            
//            PO(iq);
            
            Z = qbessel(qs[iq], delta/2-1, sqrt(x[l]*x[r])/(s[r]-s[l]));
            
//            PO(Z);

            #ifdef DEBUG
            qs.at(iq);
            #endif
            
            iq++;
            
            a = P + 2*Z + delta/2;
            b = (s[r]-s[l])/(2*(s[i]-s[l])*(s[r]-s[i]));
            G = qgamma(qs[iq], a, 1/b);
            
            #ifdef DEBUG
            qs.at(iq);
            #endif
            
            iq++;
            
//            PO(a);
//            PO(b);
//            PO(G);

            x[i] = G;

//            PIO(x);

            i = i+h;
            l = l+h;
            r = r+h;
            
        }
        jmax = jmax*2;
        h = imin;
    }
        
    vector<double> v(times.size(), v0);
    for (int i = 1; i < times.size(); ++i) {
        v[i] = x[i] * exp(-kappa * times[i]);
    }
    
    return v;
}

vector<double> poissonPathBridge(const vector<double>& times, double lambda, const vector<double>& qs)
{
    int tdim = times.size()-1;
    int m = log2(tdim);
    int h = 1 << m; // 2^m

    vector<double> N(tdim+1, 0.0);
    
    N[h] = qpois(qs[0], lambda*times[tdim]);
    
    #ifdef DEBUG
    qs.at(0);
    #endif
    
    int qindex = 1;
    int jmax = 1;

    for (int k = 1; k <= m; ++k) 
    {
        int l = 0, r = h, imin = h/2, i = imin;
        int a;
        double b, B;

        for (int j = 1; j <= jmax; ++j) 
        {
            a = N[r] - N[l]; 
            b = (times[i]-times[l])/(times[r]-times[l]);
            B = qbinom(qs[qindex],a,b);
            
            N[i] = N[l] + B;
            
            #ifdef DEBUG
            qs.at(qindex);
            #endif

            i = i + h;
            l = l + h;
            r = r + h;
            
            qindex++;
        }
        
        jmax = 2 * jmax;
        h = imin;
    }
    
    return N;
}

vector<double> sumPoissonPathBridge(const vector<double>& poissonPath, double mu_s, double sigma, const vector<double> qs)
{
    int tdim = poissonPath.size();
    
    vector<double> mus(tdim, 0.0);
    vector<double> sigmasqs(tdim, 0.0);
    
    for (int i = 1; i < tdim; ++i) 
    {
        mus[i] = poissonPath[i] * mu_s;
        sigmasqs[i] = poissonPath[i] * sigma * sigma;
    }
        
    int m = log2(tdim-1);
    int h = 1 << m; // 2^m
    
    vector<double> S(tdim, 0.0);

    S[h] = mus[h] + sqrt(sigmasqs[h]) * qnorm(qs[0]);
    
    #ifdef DEBUG
    qs.at(0);
    #endif

    
    int qindex = 1;
    int jmax = 1;

    for (int k = 1; k <= m; ++k) 
    {
        int l = 0, r = h, imin = h/2, i = imin;
        double a, b;

        for (int j = 1; j <= jmax; ++j) 
        {
            if (poissonPath[l] == poissonPath[r])
            {
                S[i] = S[l];
            }
            else
            {
                a = S[l]+mus[i]-mus[l] + (sigmasqs[i]-sigmasqs[l])*(S[r]-S[l]+mus[l]-mus[r])/(sigmasqs[r]-sigmasqs[l]);
                b = sigmasqs[i]-sigmasqs[l]-((sigmasqs[i]-sigmasqs[l])*(sigmasqs[i]-sigmasqs[l]))/(sigmasqs[r]-sigmasqs[l]);
 
                S[i] = a+sqrt(b)*qnorm(qs[qindex]);
                
                #ifdef DEBUG
                qs.at(qindex);
                #endif
            }

            i = i + h;
            l = l + h;
            r = r + h;
            
            qindex++;
        }
        
        jmax = 2 * jmax;
        h = imin;
    }
    
    return S;
}


#endif
