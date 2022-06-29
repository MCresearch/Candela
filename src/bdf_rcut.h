#ifndef BDF_RCUT_H
#define BDF_RCUT_H

#include "cell.h"

//----------------------------------
// Bond angle distribution function
//----------------------------------
class BDF_Rcut
{
	public: 
	
	BDF_Rcut(){};
	~BDF_Rcut(){};

	void Routine();
	
	private:

	void cal();

	void get_adjacent_atom_positions(
		const Cell &cel,
		const double &dv,
		const int npoints,
		double* bond_df,
		double** bond_detail);

	void bond_angle(
	    const int &ntype,
		const int &iat,
	    const int &it,
	    int *itindex,
		const double &x1,
		const double &y1,
		const double &z1,
		const int &shortest_natom,
		const double &dv,
		const bool* bonded, // mohan added 2017-02-23
		const double* posx,
		const double* posy,
		const double* posz,
		double* bond_df,
		double** bond_detail);

	private:

	double cal_angle(
		const double &a,
		const double &b,
		const double &c);

	double cal_norm(
		const double &x1,
		const double &y1,
		const double &z1,
		const double &x2,
		const double &y2,
		const double &z2);

	int tot_adj;
	int tot_count;

	double** coord_xy;
	double** aabb;
	double* ccc;
	int nmesh;

};


#endif
