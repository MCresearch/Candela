#ifndef XY_PROFILE_H
#define XY_PROFILE_H

#include "cell.h"
#include "water.h"

class XY_Profile 
{
	public: 
	
	XY_Profile();
	~XY_Profile();

	void Routine();

	private:

	void compute();
	void xy_coordinate(const Cell &cel, const int &igeo);

	int count_geometry_number;
	double** coord_xy;

};


#endif
