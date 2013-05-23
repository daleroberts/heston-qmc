// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#ifndef Scramble_H
#define Scramble_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

bool randomintzo(void);
boost::numeric::ublas::vector <bool> getdigits(int n, int m);
void getmatrixdimj(boost::numeric::ublas::matrix<int>& mj, boost::numeric::ublas::matrix<bool>& mat,int d,int mact,int m);
void createscrmatdimj(boost::numeric::ublas::matrix<int>&Ltilde,boost::numeric::ublas::matrix<int>& mj,int d,int mact,int K);
void getejvec(boost::numeric::ublas::vector <bool>& ej,int sact,int K);
void getscrambledpoints(boost::numeric::ublas::matrix<double>& P,boost::numeric::ublas::matrix<bool>& mat,int sact,int mact,int m,int K,int n);


#endif
