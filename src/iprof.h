#ifndef Iprofile_H
#define Iprofile_H

#include "cell.h"

class Iprofile
{
	public: 
	
	Iprofile(){};
	~Iprofile(){};

	void Routine();
	
	private:

	double dr;
	double surface_area;

	// calculate the ionic density profile
	// using gaussian smearing
	void cal_gauss( const Cell &cel, double* r, double** iprof, double** iprof_unit);

};


#endif
