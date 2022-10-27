#ifndef SSF_H
#define SSF_H

#include "cell.h"

//-------------------------------
// Static Structure Factor (SSF)
//-------------------------------
class SSF
{
	public: 
	
	SSF(){};
	~SSF(){};

	void Routine();
	
	private:

	// check if the ssf file exists.
	void check_file_exist( const string &name );

	// Calculate all the ssf issues inside this function. 
	void cal();

	// Calculate the static structure factor.
	void ssf_3D( const Cell &cel, float *sf ) const;

	// Sum up the terms involing phase (k \dot x)
	void sumup( float *, float *, const float dr[3] ) const; //qianrui modify

	// Calculate how many |G|.
	int cal_diff_norm(
		int *nG_1D, 
		int *which_norm
	) const;

	// Calculate the 1D SSF in sf_1D.
	void ssf_1D( 
		float *G_1D, 
		float *sf_1D, 
		const int *norm_index,
		const int *nG_1D,
		const float *sf	
	) const;

	// Output the static structure factor.
	void rank_ssf(
		const int diff_norm, 
		float *G_1D, 
		float *sf_1D,
		int *nG_1D 
	) const;
	void write_ssf( 
		const int diff_norm,
		const float *G_1D, 
		const float *sf_1D 
	) const;
	void write_smoothssf(const float *G_1D,const float *sf, const int *nG_1D, const int diff_norm, const double dG, const string filename) const;

	private:

	float dgx; // dG between nearest G points. 2pi/a1
    float dgy; // 2pi / a2
	float dgz; // 2pi / a3

	int ngx, ngy, ngz; // Number of G in each direction.
	int ngtot; // Total number of G.
	int diff_norm; // Number of differnt |G|.


};


#endif
