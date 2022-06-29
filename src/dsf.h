#ifndef DSF_H
#define DSF_H

#include "cell.h"

//-------------------------------
// Dynamics Structure Factor (SSF)
//-------------------------------
class DSF
{
	public: 
	
	DSF(){};
	~DSF(){};

	void Routine();
	
	private:

	// Calculate all the ssf issues inside this function. 
	void cal();

	int count_geo();

	// Calculate the static structure factor.
	void dsf_3D( 
		const Cell &cel1,
		const Cell &cel2,
		float *sum_exp, 
		float *sf ) const;

	// Sum up the terms involing phase (k \dot x)
	void sumup( float *sum_exp, const float dr[3] ) const;

	// Calculate how many |G|.
	int cal_diff_norm(
		int *nG_1D, 
		int *which_norm
	) const;

	// Calculate the 1D SSF in sf_1D.
	void dsf_1D( 
		float *G_1D, 
		float **sf_1D, 
		const int *norm_index,
		const int *nG_1D,
		float **sf	
	) const;

	// Output the static structure factor.
	void write_dsf( 
		const float *G_1D, 
		float **sf_1D 
	) const;

	private:

	float dg; // dG between nearest G points.
	int ngx, ngy, ngz; // Number of G in each direction.
	int ngtot; // Total number of G.
	int diff_norm; // Number of differnt |G|.

	double timetot;
	int nt;

};


#endif
