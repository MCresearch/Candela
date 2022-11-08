#ifndef CELL_H
#define CELL_H

#include "gfun.h"
#include "vec3.h"
#include "atoms.h"
#include "water.h"
#include "input.h"
// each file containing a set of type of atoms.
class Cell
{
	public: 
	
	Cell();
	~Cell();

	public:

	Atoms *atom = nullptr; // atom class
	string coordinate; //which type of coordinate
	Vector3<double> a1,a2,a3; // lattice vectors
	double volume; // volume of cell
	string file_name;
	string system_name;
	int nat; // total atom number
	int nbonds; // toal number of bonds, only used for LAMMPS
	int nbtype; // total number of bond types
	int natype; // total number of angle types
	int nangles; // total number of angles, only used for LAMMPS
	int ntype; //needed if we add new species to it.
	bool read_and_used;
	int snapshot_index;
	double snapshot_time;

	int nbands; // number of bands, used for analyzing Wannier functions
	Vector3<double> *wan_centers = nullptr; // wannier centers
	double* eig = nullptr;

	public:

	void clean();

	// transform the coordinates from direct to cartesian
	void direct2cartesian(const int &it, const int &i);

	// transform the coordinates from cartesian to direct 
	void cartesian2direct(const int &it, const int &i);	

	void add_cell_length(const int &it, const int &ia, const int &i, const int &j, const int &k, 
		double &x, double &y, double &z) const ;

	void add_cell_length(const int &it, const int &ia, const int &i, const int &j, const int &k, 
		float &x, float &y, float &z) const ;

	void Print_water(ofstream &ofs, const int &ito, const int &ith, 
		const double &norm1, const double &norm2, const double &norm3, const double &rcut_oh, const int &output_type);

	void read_wannier_centers(ifstream &ifs_wan, const int &nbands);

	void read_eig(ifstream &ifs_wan, const int &nbands);

	void read_pos_ili(ifstream &ifs_pos_ili, const int &it);

	void atom_mass();
	void init_cel(Input& Inp);

};


#endif
