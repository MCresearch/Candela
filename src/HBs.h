#ifndef HBS_H
#define HBS_H

#include "cell.h"
#include "water.h"

// hydrogen bonds 
class HBs
{
	public: 
	
	HBs();
	~HBs();

	void Routine();

	static double angle(const Cell &cel, Vector3<double> &pos1, Vector3<double> &pos2, Vector3<double> &pos3);

	static bool accepted(const Cell &cel, Water *water, const int &ito, 
	const int &ith, const int &ia, const int &ia2, double &angle0, int &H_id, double &dis);

	static bool donate(const Cell &cel, Water *water, const int &ito, 
	const int &ith, const int &ia, const int &ia2, double &angle0, int &H_id, double &dis);

	static void setup_water(const Cell &cell, Water *water);

	private:

	ofstream ofs_bond;
	ofstream ofs_hbcase;
	
	ofstream ofs_en; // record "en"vironments of oxygen atoms, mohan added 2018-04-27

	ofstream ofs_trans;
	ofstream ofs_wire;
	ofstream ofs_zundel; // mohan add 2017-12-10
	int old_ion_index;

	void search_hydrogen_bonds(const Cell &cell, const int &igeo);

	double avg_aion; // average accepted HBs for ions
	double avg_dion; // average donating HBs for ions
	int avg_count_ion; // average counting for ions

	double avg_awater; // average accepted HBs for water molecules
	double avg_dwater; // average donating HBs for water molecules
	int avg_count_water; // average counting for water molecules

	// angle distribution
	double* angle_ion;
	double* angle_water; // for all water molecules
	
	double* accept_dis;
	double* donate_dis;

	int* accept_count;
	int* donate_count;
	int* HB_count;
	int** HB_count_acc;

	int nacc; // number of accepted HBs of current ion
	int ndon; // number of donating HBs of current ion
	int last_nacc; // number of accepted HBs of the last snapshot of previous ion
	int last_ndon; // number of donating HBs of the last snapshot of previous ion
	int snapshot_index; // snapshot index of current ion
	double snapshot_time; // snapshot time of current ion
	int last_snapshot_index; // snapshot index of the last snapshot of previous ion
	double last_snapshot_time; // snapshot time of the last snapshot of previous ion
	int last_nH;
	int* last_Hindex;

	int count_geometry_number;

	// condiering broken HBs along z direction
	int nz;
	double* z1d;
	double* z2d;



};


#endif
