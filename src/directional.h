#ifndef DIRECTIONAL_H
#define DIRECTIONAL_H

#include "cell.h"
#include "water.h"

// Define directional H-bonds
class Directional
{
	public: 
	
	Directional();
	~Directional();

	void Routine();

	private:

	void compute_dist2D(const Cell &cel);
	void locate(const double &angle0, const double &dis);
	int count_geometry_number;
	double** dist2D;

	double drx;
	double dry;

	double dux;
	double duy;
	double duz;

	double upper_natom;
	double lower_natom;
};


#endif
