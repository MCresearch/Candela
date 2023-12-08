#ifndef MSD_MULTIPLE_H
#define MSD_MULTIPLE_H

#include "cell.h"
#include "water.h"
#include "mj.h"

class MSD_Single
{
	public:
	
	MSD_Single();
	~MSD_Single();
	
	void allocate(const double &t0);
	void deallocate();
	
	double t0; // starting time
	double t1; // ending time

	int count_msd; // count msd

	double* t; // time
	double* msd; // mean square displacement
	int* mmm; // counts per dt
	int ndim; // dimension of msd

	int natom; // number of atoms that need MSD 
	Vector3<double>* pre_wpos; // Positions of water molecules in previous step, has boundary
	Vector3<double>* wpos;     // Positions of water molecules, no boundary
	Vector3<double>* wpos0;    // Positions of water molecules at time 0, no boundary

	int saved_ion; // saved ion index
	double sx, sy, sz; // saved x, y, z
	double mx, my, mz; // accumulation of jump distances of proton transfer

	double xmc0, ymc0, zmc0; // mass center recored in the first step
};

// Multiple Mean Square Displacement
class MSD_Multiple
{
	public: 
	
	MSD_Multiple();
	~MSD_Multiple();

	void Routine();

	private:

	void compute_msd(const Cell &cel, const int &igeo);
	void each_msd(const Cell &cel, const int &ito, const int &it_select, MSD_Single &ms, const Water* water);
	int round(double r);
	MSD_Single* ms; // mean square displacement array	
};


#endif
