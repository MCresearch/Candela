#ifndef DIST_H
#define DIST_H

#include "cell.h"
#include "water.h"

// 3D distribution
class Dist
{
	public: 
	
	Dist();
	~Dist();

	void Routine();

	private:

	void compute_dist2D(const Cell &cel);
	void locate(const double &angle0, const double &dis);
	void compute_dist3D();
	int count_geometry_number;
	double** dist2D;
	double*** dist3D;

	double drx;
	double dry;

	double dux;
	double duy;
	double duz;

	double upper_natom;
	double lower_natom;
};


#endif
