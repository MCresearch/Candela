#ifndef DIELECTRIC_H
#define DIELECTRIC_H

#include "cell.h"

class Dielectric
{
	public:
    
	Dielectric();
    ~Dielectric();

	void Routine(); // compute dipole moments

	private:
	
	int natom;
	double* dipole_x;
	double* dipole_y;
	double* dipole_z;
	double* dipole_tot;
	double epsilon_water;
	int count_epsilon;

	double dd; // delta of dielectric constant
	double dcut; // cutoff for dielectric constant
	double* dis_epsilon;


	// to compute dielectric radial distribution functions
	private:

	void dipole_rdf();
	void pairs( const Cell &cel, double *gr);
	void compute_GkR(double *gr);

	// setup parameters for computing dipole RDF
	double dr; // delta r in discretizing dipole RDF
	double rcut; // rcut of dipole RDF
	int nmesh;	// number of real-space points of dipole RDF
	int count_geometry_number;
	double rho_ion;    

	// macroscopic moment
	double M2_intra;
	double M2_inter;
	double M2;
	double Mx;
	double My;
	double Mz;

	// parameters to compute GkR
	double *GkR; 
	double volume;

	// for 2D plot
	double** coord_xy;

};

#endif //Dielectric

