#ifndef BDF_H
#define BDF_H

#include "cell.h"

//----------------------------------
// Bond angle distribution function
//----------------------------------
class BDF
{
	public: 
	
	BDF(){};
	~BDF(){};

	void Routine();
	
	private:

	void cal();

	void get_adjacent_atom_positions(
		const int &shortest_natom,
		const Cell &cel,
		const double &dv,
		const int npoints,
		int* bond_df,
		int** bond_detail,
		int* bond_type);

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
		const double* posx,
		const double* posy,
		const double* posz,
		int* bond_df,
		int** bond_detail,
		int* bond_type);

	void shift2back(const int &pairs,
		const int &ip_start,
		double* dist,
		double* posx,
		double* posy,
		double* posz,
		int* itindex,
		int* iaindex);

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

	bool is_ion(const int &it, const int &ia);

};


#endif
