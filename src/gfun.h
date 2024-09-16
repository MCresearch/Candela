#ifndef GFUN_H
#define GFUN_H

#include <cstdlib>
#include <unistd.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <cassert>
#include <ctime>

#ifdef __MPI
#include <mpi.h>
#endif

#include "vec3.h"

using namespace std;

extern ofstream ofs_running;
extern int NPROC;
extern int RANK;


//==========================================================
// Constants 
//==========================================================
extern double PI;
extern double KB; // Boltzmann constant
extern double EPSILON0; // vacuum permitivity
extern double BOHR; // atomic unit of length

//==========================================================
// GLOBAL FUNCTION :
//==========================================================
template <class T>
static void READ_VALUE(ifstream &ifs, T &v)
{
    ifs >> v;
    ifs.ignore(20000, '\n');
    return;
}

bool SCAN_BEGIN(ifstream &ifs, const string &TargetName, const bool restart=1);
// Mohan warning : the last term can't be written as const bool &restart,
// I don't know why.



void TITLE(const string &class_name,const string &function_name);

void QUIT(const string &name);

template <class T>
inline void ZEROS(T a, const int &size)
{
	assert(size>0);
	for(int i=0; i<size; ++i) a[i]=0;
}

double shortest(const double &tmp_pos, const double &pos, const double &celldm);

double Polynomial_Interpolation
(
    const double *table,
    const int &table_length,
    const double &table_interval,
    const double &x				// input value
);

double Polynomial_Interpolation_xy
(
    const double *xpoint,
    const double *ypoint,
    const int table_length,
    const double &x             // input value
);

double distance
(
	const Vector3<double> &pos1, 
	const Vector3<double> &pos2, //mohan add 2016-10-21
	const double &a1,
	const double &a2,
	const double &a3
);
//qianrui add 2020-1-6
//check if x==y
#define ifnecheckv(x,y)\
{\
	if(!ifne(x,y))\
	{\
		cout<<"please check variation "<<#x<<" & "<<#y<<" in "<<__FILE__<<':'<<__LINE__<<'!'<<endl;\
		exit(0);\
	}\
} 
#define ifnechecke(x,y)\
{\
	if(!ifelt(x,y))\
	{\
		cout<<"please check variation "<<#x<<" & "<<#y<<" in "<<__FILE__<<':'<<__LINE__<<'!'<<endl;\
		exit(0);\
	}\
} 
#define ifnecheckp(x,y,n)\
{\
	if(!ifne(x,y,n))\
	{\
		cout<<"please check pointer "<<#x<<" & "<<#y<<" in "<<__FILE__<<':'<<__LINE__<<'!'<<endl;\
		exit(0);\
	}\
} 
template<class T>//qianrui add 2020-1-6
bool ifne(const T x, const T y)
{
	if(x!=y)
	{
		cout<<x<<" != "<<y<<endl;
		cout<<"Error in read data!"<<endl;
		return false;
	}
	return true;

}

template<class T>//qianrui add 2020-1-6
bool ifelt(const T x, const T y)
{
	double error=abs(double(x)-double(y));
	if(x>0&&y>0)
	{
		double re1=error/x;
		double re2=error/y;
		if(re1<=0.0001&&re2<=0.0001)
			return true;
		else
		{
			cout<<x<<" != "<<y<<endl;
			cout<<"Error in read data!"<<endl;
			return false;
		}
	}
	else if(x<0&&y<0)
	{
		double re1=-1*error/x;
		double re2=-1*error/y;
		if(re1<=0.0001&&re2<=0.0001)
			return true;
		else
		{
			cout<<x<<" != "<<y<<endl;
			cout<<"Error in read data!"<<endl;
			return false;
		}
	}
	else if(double(x)==0.0&&double(y)==0.0)
		return true;
	else
	{
		cout<<x<<" != "<<y<<endl;
		cout<<"Error in read data!"<<endl;
		return false;
	}
	return true;
}
template<class T>//qianrui add 2020-1-6
bool ifne(T*x,T*y,int n)
{
	for(int i=0;i<n;i++)
	{
		if(x[i]!=y[i])
		{
			cout<<i<<" : "<<x[i]<<" != "<<y[i]<<endl;
			cout<<"Error in read data!"<<endl;
			return false;
		}
	}
	return true;
}

//return the conjugate of one complex number.
template<class T> //qianrui add 2020-1-6
complex<T> conjugate(complex<T> & a)
{
	complex<T> y;
	T im=-1*a.imag();
	y.real(a.real());
	y.imag(im);
	return y;
}
//convert string to int//qianrui add on 2020-2-4
int str2int(const string str);
string int2str(const int y);
double str2dou(const string str);
string dou2str(const double y);
template<class T>
void pclean(T* &p)
{
	if(p!=nullptr)
	{
		delete[] p;
		p=nullptr;
	}
	return;
}

void searchead(ifstream &,string &,const string,int);
double mysecond(long int& time);
#endif
