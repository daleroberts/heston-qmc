// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#pragma once

using namespace std;

double Mean(vector<double>& vec) {
    double mean = 0.0;
    for (vector<double>::iterator iter = vec.begin(); iter != vec.end(); ++iter)
        mean += *iter;
    mean = mean / vec.size();
    
    return mean;
}

double StdError(vector<double>& vec) {
    double mean = Mean(vec);
    
    double var = 0.0;
    for (vector<double>::iterator iter = vec.begin(); iter != vec.end(); ++iter)
        var += (*iter - mean)*(*iter - mean);
    var = var / (vec.size() * (vec.size()-1));
    
    return sqrt(var);
}
