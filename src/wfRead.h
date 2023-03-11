#ifndef WFREAD_H
#define WFREAD_H
#include "input.h"
#include "wfqe.h"
#include "wfPWmat.h"
#include "wavefunc.h"
#include "wfABACUS.h"

//created by qianrui on 2020-2-5
//A class to operate some functions about WF
class WfRead  
{
	public:
	WfRead(){
	wfpt=nullptr;
	};
	~WfRead(){};
	WfRead(Wavefunc &wf);
	void setval(Wavefunc &wf);
	void Init();
	void readWF(int ik);
	void readOCC(int ik);
	void readvmatrix(const int ik);
	void calvmatrix();
	void clean();
	void cleanclass();
	void ignore(int );
	void ignoreOCC(int ik);
	Wavefunc *wfpt;
	WfQE wfqe;
	WfPWmat wfpwmat;
	WfABACUS wfabacus;
};
#endif
