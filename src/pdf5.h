#ifndef PDF5_H
#define PDF5_H

#include "cell.h"
#include "Honeycutt.h"

class PDF5
{
	public: 
	
	PDF5(){};
	~PDF5(){};

	void Routine();
	
	private:

	// calculate the pair distribution function.
	void cal();

	// calculate the paris under periodic boudary conditions.
	void periodic_pairs( const Cell &cell, double** gr, const int &option);

	private:

	int nn; //number of neighbours needed

	double dr; // delta r
	double rcut; // radius cutoff
	int nmesh; // number of mesh points.
	int count_geometry_number;

	double rho_ion;
	double dg;
	int ng;

};


#endif
