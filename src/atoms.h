#ifndef ATOMS_H
#define ATOMS_H

#include "gfun.h"
#include "vec3.h"

// each element containing a set of data
class Atoms
{
	public:

	Atoms();
	~Atoms();

	public:

	string id;//atom id.
	string pot_file;//pseudopotential file name.
	int na; //number of atoms.
	Vector3<double>* pos; //atom positions in Cartesian coordinates.
	Vector3<double>* posd; //atom positions in Direct coordinates.
	Vector3<double>* vel; //atom velocities.
	double mass; // mass of element in a.u., added by mohan on 2016-10-26
	double charge; // charge of element, only used in LAMMPS, 2016-10-26
	double* pos_ili; // distance to ili, mohan added 2017-04-08

	bool allocate_pos;
	bool allocate_posd;
	bool allocate_vel;
	bool allocate_pos_ili; // mohan added 2017-04-08

	public: 

	// read the 'fractional' coordinates of atoms
	// or the 'cartesian' coordinates of atoms.
	void read_pos(ifstream &ifs,bool frac);
	void read_pos_2(ifstream &ifs,bool frac);
	void read_pos_3(ifstream &ifs); // mohan add 2015-07-24 for ABACUS
	void read_pos_4(ifstream &ifs, Vector3<double> &a1, Vector3<double> &a2, Vector3<double> &a3); // mohan add 2016-08-02 for QE 
	void read_pos_5(ifstream &ifs, Vector3<double> &a1, Vector3<double> &a2, Vector3<double> &a3, double &lat_const);
	void read_vel(ifstream &ifs);

};

#endif
