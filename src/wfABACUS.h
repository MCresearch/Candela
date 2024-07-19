#ifndef WFABACUS_H
#define WFABACUS_H
#include "wavefunc.h"
#include "gfun.h"
#include "binfstream.h"
#include <fstream>
//read wave function in ABACUS
class WfABACUS
{
	public:
	WfABACUS(){
	wk=nullptr;
	energy=nullptr;
	occ=nullptr;
	};
	~WfABACUS(){
		pclean(wk);
		this->clean();
	};
	void readWF(Wavefunc &,  int&);
	void readOUT(Wavefunc &);
	void readOCC(Wavefunc &, int&);
	void clean();
	private:
	binfstream rwswf;
	ifstream ifsocc;
	ifstream ifskwt;
	
	private:
	Vector3<double> axes[3];
	int nband;
	double *energy,*occ;
	int nkpoint;
	int ngtot;
	int num_ele;
	double * wk;
	double alat;
};
#endif
