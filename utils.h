// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#pragma once

#ifdef DEBUG
    #define PO(e) std::cout << #e << ": " << std::setprecision(10) << e << std::endl
    #define PIO(v) std::cout << #v << ": "; copy(v.begin(), v.end(), ostream_iterator<double>(cout, " ")); std::cout << " (size: " << v.size() << ")" << std::endl
#else
    #define PO(e) {}
    #define PIO(v) {}
#endif

double pnorm(double z);
double qnorm(double p);
double qchisq(double p, double df, double ncp = 0.0);
double qpois(double p, double lambda);
int qbessel(double p, double nu, double z);
double qgamma(double p, double shape, double scale);
double qbinom(double p, double n, double prob);
double bscall(double S, double K, double t1, double t2, double r, double sigma);
bool isPowerOfTwo(int x);
int log2(int val);

struct Dummy {
	double _value;
    
	Dummy(double value) : _value(value) {};
    
	double operator()(void) {
		return _value;
	}
};

struct Result {
    double mean;
    double stderr;
    
    Result(double mean_, double stderr_): mean(mean_), stderr(stderr_) {}
};
