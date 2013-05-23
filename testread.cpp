// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <stdio.h>

#include <boost/numeric/ublas/matrix.hpp>

#include "matlab.h"

using namespace std;

int main(int argc, char *argv[]) {
    
    string filename = "rep2.dlm";
    boost::numeric::ublas::matrix<double> m;

    ReadMatlabMatrix(filename, m);
    
    for (int i = 0; i < m.size1(); i++) {
        for (int j = 0; j < m.size2(); j++) {
            printf("%.6f ", m(i,j));
        }
        printf("\n");
    }
}
