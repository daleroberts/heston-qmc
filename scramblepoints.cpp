// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

bool randomintzo(void)
{
	double u=(rand()+1.0)/static_cast<double>(RAND_MAX+2.0);
	//double u=0;
	if(u<=0.5)
		return(0);
	else
		return(1);
};


boost::numeric::ublas::vector <bool> getdigits(int n, int m)
{
	double nt;
	bool ndig;
	boost::numeric::ublas::vector <bool> ndigits(m);

	nt=((double)n)/(pow((double) 2.0,(double)m));

	for(int i=1;i<(m+1);i++)
	{
		//cout << i << endl;
		//cout << 2.0*nt << endl;
		ndig= floor(2.0*nt);
		//cout << ndig << endl;
		nt=2.0*nt-ndig;
		//cout << nt << endl;
		ndigits[m-i]=ndig;
		//cout << ndigits[m-i] << endl;
	};
	return(ndigits);
};


void getmatrixdimj(boost::numeric::ublas::matrix<int>& mj, boost::numeric::ublas::matrix<bool>& mat,int d,int mact,int m)
{
	for(int i=0;i<mact;i++)
	{
		for(int k=0;k<mact;k++)
		{
			mj(i,k)=mat(d*m+i,k);
		};
	};
};


void createscrmatdimj(boost::numeric::ublas::matrix<int>&Ltilde,boost::numeric::ublas::matrix<int>& mj,int d,int mact,int K)
{
	//multiply L_j C_j(which is mj), where L_j is K*mact, lower triangular and filled with random 0 ones
	boost::numeric::ublas::matrix<int> Ltilde2(K,mact);
	
	for(int i=0;i<mact;i++)
	{
		for(int k=0;k<i;k++)
		{
			Ltilde2(k,i)=0;
		};

		Ltilde2(i,i)=1;
		for(int k=i+1;k<K;k++)
		{
			Ltilde2(k,i)=randomintzo();
		};
	};

	//cout << Ltilde2 << endl;


	//now multiply Ltilde with C_j(=mj) and save teh result in Ltilde
	for(int i=0;i<K;i++)
	{
		for(int j=0;j<mact;j++)
		{
			Ltilde(i+d*K,j)=0;
			for(int k=0;k<mact;k++)
			{
				Ltilde(i+d*K,j)+=Ltilde2(i,k)*mj(k,j);
			};
		};
	};


};

void getejvec(boost::numeric::ublas::vector <bool>& ej,int sact,int K)
{
	for(int i=0;i<(sact*K);i++)
	{
		ej(i)=randomintzo();
	};
	//cout << ej << endl;
};

void getscrambledpoints(boost::numeric::ublas::matrix<double>& P,boost::numeric::ublas::matrix<bool>& mat,int sact,int mact,int m,int K,int n)
{
	boost::numeric::ublas::vector <bool> ndigits(mact); // to contain the result of the matrix vector multiplication
	boost::numeric::ublas::vector <int> xdigits(K); // to contain the digits of the point

	boost::numeric::ublas::matrix <int> mj(mact,mact); //used to store the matrix of the original point set for a particular dimension
	boost::numeric::ublas::vector <bool> ej(sact*K); //used to shift the digits

	boost::numeric::ublas::matrix <int> Ltilde(sact*K,mact); //to contain the products 

	for(int d=0;d<sact;d++)
	{
		getmatrixdimj(mj, mat,d,mact,m); //mj contains C_j
		createscrmatdimj(Ltilde,mj,d,mact,K);
			   
	};
	getejvec(ej,sact,K); //get the whole vector ej(1:K) for dim 1, ej(K+1:2*k) for dim 2 etc



	for(int ni=0;ni<n;ni++)
	{
		   
		ndigits=getdigits(ni,mact); //contains the b-adic expansion of ni
		  
		for(int d=0;d<sact;d++) //for dimension d
		{

			for(int i=0;i<K;i++) //obtain the K digits for point P(ni,d) 
			{
				xdigits[i]=0;
				for(int j=0;j<mact;j++)
				{
							//cout << mat(d*m+i,j) << endl;
							//cout << ndigits[j] << endl;
					xdigits[i]+=(Ltilde(d*K+i,j)*ndigits(j));
				};
				   //cout << xdigits[i] << endl;
				xdigits[i]=(xdigits(i)+ej(d*K+i))%2;
				   //cout << xdigits[i] << endl;
				   //cout << endl;
				   //cout << xdigits[i] << endl;
			};
			P(d,ni)=0;
			for(int i=0;i<K;i++)
			{
				   //cout << (double) pow(float(2),i+1) << endl;
				   //cout << (double) xdigits[i] << endl;
				P(d,ni)+=(double(xdigits(i))/((double) pow(float(2),i+1)));
			};

					
		};
				
	};



};
