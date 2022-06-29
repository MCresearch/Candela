#ifndef ISF_H
#define ISF_H

#include "cell.h"

//--------------------------------------------------
// Intermediate scattering function 
// F(q,t) = 1/N * <
// ( \sum_{j=1,N} e^{-iqR_j(t)}) *
// ( \sum_{i=1,N} e^{ iqR_i{0}})
// >
//--------------------------------------------------
class ISF
{
	public: 
	
	ISF();
	~ISF();

	void Routine();

	private:

	void selectq();
	void cal_Fqt(float *Fqt, float *qx, float *qy, float *qz);
	void output_Fqt(const string &out_file, const float *Fqt)const;

	float cal_eiqR(float *qx, float *qy, float *qz, const Cell &cel1, const Cell &cel2);
	
	private:

	int nconfig;
	int ncorrelation;
	int nat;
	int nat1; // selected atoms, for parallel purpose
	int nat2;

	// for the q vectors
	float target_q;
	int count_q;
	float* qsaved_x;
	float* qsaved_y;
	float* qsaved_z;


};


#endif

