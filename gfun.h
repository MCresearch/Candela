#ifndef GFUN_H
#define GFUN_H

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

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
    ifs.ignore(150, '\n');
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

#endif
