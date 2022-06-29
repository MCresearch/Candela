#ifndef DENSITY2D_H
#define DENSITY2D_H

#include "cell.h"

class Density2D
{
public:
    
	Density2D();
    ~Density2D();
    
	void Routine(); 

	private:

	void covalent(const Cell &cel, const int &igeo);

	double** coord_xy;
	
	int nx;
	int ny;
	double x0;
	double y0;
	double dx;
	double dy;

};

#endif //Density2D
