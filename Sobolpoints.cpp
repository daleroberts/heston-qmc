// 
// Dale Roberts <dale.o.roberts@gmail.com>
//


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <fstream>
using std::ifstream;

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

void getSobolpoints(boost::numeric::ublas::matrix<bool>& mat)
{
	   ifstream indata; // indata is like cin
	   double num; // variable for input value

	   int m= 16, s = 1000; //to be changed TOGETHER with the filename! hardcoded, as it corresponds to the filename
	   int ccount = 0, rcount = 0;
	   
	   //boost::numeric::ublas::matrix < bool > mat(m*s, m);

      indata.open("sobmatm16s1000vers1.dat"); // opens the file
       if(!indata) { // file couldn't be opened
          cerr << "Error: file could not be opened" << endl;
          exit(1);
       }

      indata >> num;
       while ( !indata.eof() ) 
	   { // keep reading until end-of-file
          //cout << "The next number is " << num << endl;
		   if ( num == 1.0 )
			mat(rcount, ccount) = 1;
		   else
			   mat(rcount, ccount) = 0;
		   ccount++;

		   if ( ccount >= m )
		   {
			   ccount = 0;
			   rcount++;
			   //std::cout<<rcount;
		   }

          indata >> num; // sets EOF flag if no value found
       }
       indata.close();
}; 
