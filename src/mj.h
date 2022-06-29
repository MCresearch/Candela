#ifndef MJ_H
#define MJ_H

#include "cell.h"
#include "water.h"

class ProtonTransfer
{
	public: 

	ProtonTransfer();
	~ProtonTransfer();

	bool use_pt;
	void setup_PT();
	string which_type_pt(const int &snapshot_i);

	int npt; // number of proton transfer events

	// data set 1
	int* snapshot_index; 
	double* snapshot_time;
	int* ion_index;
	int* nacc; // number of accepted HB at snapshot 'snapshot_index'
	int* ndon;

	// data set 2
	int* snapshot_index_pt; 
	double* snapshot_time_pt;
	int* ion_index_pt;
	int* nacc_pt; // number of accepted HB in the final snapshot of this ion
	int* ndon_pt;

	// data set 3
	int* indexH;

	string* type_pt;
	int* ds_pt; // how many snapshots it takes to transfer a proton
	double* dt_pt; // how much time it takes to transfer a proton 
	int* ds_min; // range (minimum snapshot) to analyze the compression of water wire
	int* ds_max; // range (maximum snapshot) to analyze the compression of water wire


};

// multiple jump 
class MJ
{
	public: 
	
	MJ();
	~MJ();

	void Routine();
	int ndim;
	double dt_max;

	int *npt; // number of proton transfer
	int *npt_hyper; // number of proton transfer via hyper
	int *npt_tetra; // number of proton transfer via tetra

};


#endif
