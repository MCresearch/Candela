#ifndef MSD_H
#define MSD_H

#include "cell.h"
#include "water.h"
#include "mj.h"


// Mean Square Displacement
class MSD
{
	public: 
	
	MSD();
	~MSD();

	void Routine();

	private:

	void compute_msd(const Cell &cel, const int &igeo, ofstream &ofs_msd);

	int count_geometry_number;
	double* msd;
	double* msd2;
	int nmsd;

	Vector3<double> pre_pos;  // the position of ion in previous step, has boundary
	Vector3<double> ion_pos; // the position of ion, no boundary
	Vector3<double> ion_pos0; // the position of ion at time 0, no boundary

	Vector3<double>* pre_wpos; // the positions of water molecules in previous step, has boundary
	Vector3<double>* wpos; // the positions of water molecules, no boundary
	Vector3<double>* wpos0; // the positions of water molecules at time 0, no boundary

	Vector3<double>* pre_wpos2; // the positions of water molecules in previous step, has boundary
	Vector3<double>* wpos2; // the positions of water molecules, no boundary
	Vector3<double>* wpos02; // the positions of water molecules at time 0, no boundary
	bool allocate_water;

	int count_msd;
	double start_time;

	ProtonTransfer PT; // mohan added 2017-09-06
};


#endif
