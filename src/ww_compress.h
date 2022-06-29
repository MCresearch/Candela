#ifndef WW_COMPRESS_H
#define WW_COMPRESS_H

#include "cell.h"
#include "water.h"
#include "mj.h"


// Compression of water wires 
class WW_Compress 
{
	public: 
	
	WW_Compress();
	~WW_Compress();

	void Routine();

	private:

	void count_compression(const Cell &cel, const int &igeo);
	void neighbour1(const Cell &cel, Water *water, const int &ito, const int &ith);

	int count_geometry_number;
	double count_min;

	double* ww;
	double* ww1;
	double* ww2;
	double* ww_h2o;
	double* ww1_h2o;
	double* ww2_h2o;
	double dr;
	double rcut;
	int nmesh;

	double* ww_all;
	double* ww_h2o_all;
};


#endif
