#ifndef WAVEFUNC_H
#define WAVEFUNC_H

#include <complex>
#include "gfun.h"
using namespace std;
//Class to store wavefunction.
//By qianrui 2020-1-6
class Wavefunc
{
	public: 
		Wavefunc()
		{
			Wavegg=nullptr;
			gkk_x=gkk_y=gkk_z=nullptr;
			occ=eigE=nullptr;
			kpoint_x=kpoint_y=kpoint_z=0;
			wk=0;
		};
		~Wavefunc();
		Wavefunc(Wavefunc &);
		void print(int );
		void checknorm(int,int);
	public:
		//For ABACUS, QE:
		//    <Wavegg|Wavegg> = 1
		//For PWmat:
		//    <Wavegg|Wavegg>*volume = 1 
		complex <double> * Wavegg = nullptr; //in bohr^-3/2 or 1
		double factor;// in unit 1 or bohr^3
		//Wavegg^2*factor in unit 1
		double *gkk_x = nullptr,*gkk_y = nullptr,*gkk_z = nullptr;//in bohr^-1
		double kpoint_x,kpoint_y,kpoint_z;//in bohr^-1
		double *occ = nullptr; // normalized to 1
		double *eigE = nullptr; //in eV
		int nband,ngtot;
		double wk;//weight of kpoint, normalized to 1
	public:
		int ig0 = 0; //ig when gkk = (0,0,0)
};
//add fftw in future.

#endif
