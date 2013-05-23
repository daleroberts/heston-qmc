// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#pragma once

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

using boost::numeric::ublas::column_major;
using boost::numeric::ublas::matrix;
using std::vector;

extern "C" {
void sobol_(double *QN, int *n, int *dimen, double *quasi, int *ll, int *count, int *SV, int *iflag, int *iseed, int *init, int *transform);
}

void InitSobol(matrix<double, column_major>& QN) {
    int n = QN.size1(), dimen = QN.size2(), maxbit = 30;
    int iflag = 0, transform = 0, iseed = 1, init = 1, ll;
    matrix<int, column_major> SV(dimen, maxbit);
    vector<double> quasi(dimen);
    int count;

    sobol_(&QN.data()[0], &n, &dimen, &quasi[0], &ll, &count, &SV.data()[0], &iflag, &iseed, &init, &transform);
}

void ScrambleSobol(const int seed, matrix<double, column_major>& QN) {
    int n = QN.size1(), dimen = QN.size2(), maxbit = 30;
    int iflag = 1, transform = 0, iseed = seed, init = 1, ll;
    matrix<int, column_major> SV(dimen, maxbit);
    vector<double> quasi(dimen);
    int count = 0;

    sobol_(&QN.data()[0], &n, &dimen, &quasi[0], &ll, &count, &SV.data()[0], &iflag, &iseed, &init, &transform);
}
