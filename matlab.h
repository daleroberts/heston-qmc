// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#ifndef volmodels_matlab_h
#define volmodels_matlab_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>

bool ReadMatlabMatrix(const std::string &filename, boost::numeric::ublas::matrix<double> &m )
{
    using namespace std;
    
    ifstream infile(filename);
    vector<double> values;
    
    if ( infile.bad() ) 
    {
        std::cerr << "Error:" << __FUNCTION__ << ":"
                  << "Cannot open file : " << filename << "\n";
        return false;
    }
    
    string line;
    int linecount = 0;

    while (getline(infile, line)) 
    {
        istringstream linestream(line);
        string item;
        double value;
    
        while (getline(linestream, item, ',')) 
        {
            istringstream itemstream(item);
            itemstream >> value;
        
            values.push_back(value);
        }
        
        linecount++;
    }
    
    int rows = linecount;
    int cols = values.size() / rows;
    
    m.resize(rows, cols);
    
    int k = 0;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            m(i,j) = values[k++];
    
    return true;
}

#endif
