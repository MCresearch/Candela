#ifndef MDP_H
#define MDP_H

#include "cell.h"

// Copmute mean density profile (MDP) based on
// instantaneous liquid interfaces (ILI)
class MDP
{
	public: 
	
	MDP();
	~MDP();

	void Routine();

	static void which_surface(int &six, int &siy,
		const double &norm1, const double &norm2, 
		const double &posx, const double &posy, const double &posz,
		double& drx, double& dry, double& drz, double** interface_in);

	static void read_ili(ifstream &ifs, double **interface, double ***gradient);
	
	private:

	double** joint_cd; // joint conditional distribution
	double** record_water; // record the information for water molecules
	double** record_ion; // record the information for ion

	double dr;
	double rcut;
	int nmesh;
	double du1;
	double du2;

	double** interface;
	double*** gradient;

	void proximity(ifstream &ifs, const Cell &cel, ofstream &ofs3);
	
	void joint_conditional_distribution(const double &aaa, const int &ith, const Cell &cel, 
		const double& posx, const double& posy, const double &posz, const double *grad);

	double ave_accept;
	double ave_donate;
	int ave_count;

	double* ili_dis; // distances of atoms to ILI 

};


#endif
