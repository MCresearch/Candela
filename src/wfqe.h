#ifndef WFQE_H
#define WFQE_H
#include "wavefunc.h"
#include "gfun.h"
#include "binfstream.h"
#include <fstream>
//read wave function in QE
class WfQE
{
	public:
	WfQE(){
	wk=nullptr;
	};
	~WfQE(){
	pclean(wk);
	};
	void readGKK(Wavefunc &, int&);
	void readWF(Wavefunc &,  int&);
	void readOUT(Wavefunc &);
	void readOCC(Wavefunc &, int&);
	//qianrui 2020-2-18
	void readOCC2(Wavefunc &, int&);
	void readWF2(Wavefunc &,  int&);
	void readOUT2(Wavefunc &);
	void clean();
	private:
	binfstream rwsgkk;
	binfstream rwswf;
	ifstream ifsocc;
	ifstream ifskwt;
	
	private:
	Vector3<double> axes[3];
	int nband;
	int nkpoint;
	int ngtot;
	int num_ele;
	double * wk;
	double alat;
};
#endif
