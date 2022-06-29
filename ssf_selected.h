#ifndef SSF_SELECTED_H
#define SSF_SELECTED_H

#include "cell.h"

//-------------------------------
// Static Structure Factor (SSF)
//-------------------------------
class SSF_Selected
{
	public: 
	
	SSF_Selected(){};
	~SSF_Selected(){};

	void Routine();
	
	private:

	void generate_input();
	void select_g(ifstream &ifs2);
	void cal(ifstream &ifs);

	// Calculate the static structure factor.
	void ssf_3D( const Cell &cel, float *sf ) const;
	// Sum up the terms involing phase (k \dot x)
	void sumup( float *sum_exp, const float dr[3] ) const;

	private:

	float dgx; // dG between nearest G points. 2pi/a1
    float dgy; // 2pi / a2
	float dgz; // 2pi / a3

	int ngx, ngy, ngz; // Number of G in each direction.
	int ngtot; // Total number of G.
	int diff_norm; // Number of differnt |G|.

	float* gx;
	float* gy;
	float* gz;
	float* norm_value;

};


#endif
