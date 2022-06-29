#ifndef TETRA_ORDER_H
#define TETRA_ORDER_H

#include "cell.h"
#include "water.h"

// tetrahedral order parameter
class TOP
{
	public: 
	
	TOP(){tetra_order_parameter=0.0;ave_adj=0.0;count_adj=0;};
	~TOP(){};

	void Routine();

	private:

	void cal_tetra_order_para(const Cell &cell);

	double ave_adj;
	int count_adj;
	int count_adj2;
	double tetra_order_parameter;

	double* dis_tetra;
	int* dis_count;
	double dr;
	double rcut;
	int nmesh;

	// TOP in terms of water density
	double* den_tetra;
	int* den_count;
	double dr_den;
	double dcut;
	int nden;
	double density; // density of this snapshot

};


#endif
