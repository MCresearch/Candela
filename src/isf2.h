#ifndef ISF2_H
#define ISF2_H
#include <sstream>
#include <iostream>
#include <fstream>
#include "input.h"
#include "vec3.h"
//A program to calculate intermediate scattering function.
//By qianrui on 2019-1-9
class ISF2
{
	public:
	ISF2(){};
	~ISF2(){};
	
	void Routine();
	
	private:
	//calculate the intermediate scattering function. 
	void cal();
	
	private:
	

};
template<class T>
bool getrhok(int,Vector3<T>*,int,double*,int&);
bool targetk(int,int,int);
double gettime(int);
int getindex(int,int,int,int);
int getigeo(int);
void geoignore();
void writeisf(double *,int *);
#endif

/*
INPUT para:
calculation	isf2
geo_in_type
geo_directory
isf_outfile

geo_1,2( make sure: interval*(nt+nT)<=geo_2-geo_1)
geo_interval

natom
isf_nt1
isf_nt2
dt_snapshots (ps)

isf_target_q (A-1)
isf_ngx,y,z
isf_dgx,y,z (A-1)

*/
