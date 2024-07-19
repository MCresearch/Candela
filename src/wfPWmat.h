#ifndef WFPWMAT_H
#define WFPWMAT_H
#include "wavefunc.h"
#include "gfun.h"
#include "binfstream.h"
#include <fstream>
//Read wave funciton in PWmat
class WfPWmat
{
	public:
	WfPWmat();
	~WfPWmat();
	void Init(Wavefunc &);//before ik loop
	void readOCC(Wavefunc &,  int&);//under ik loop
	void readGKK(Wavefunc &,  int&);//under ik loop
	void readWF(Wavefunc &,  int&);//under ik loop
	void clean();
	private:
	binfstream rwsgkk;
	binfstream rwswf;
	ifstream ifsocc;
	
	private:
	int *ngtotnod_9,*ngtotnod_9_t;
	double *kpoint_x, *kpoint_y, *kpoint_z;
	double *energy,*occ;
	Vector3<double> axes[3];
	int n1,n2,n3;
	int nband;
	int nnodes;
	int nkpoint;
	int num_ele;
	double vol;
	double *wk;//weight of kpoint
	bool init;
	int ngtot;
	double alat;//crystal parameter
	int mg_nx;//max g points
};
#endif
