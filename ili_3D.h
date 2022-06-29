#ifndef ILI3D_H
#define ILI3D_H

#include "cell.h"

// plot the 3D instantaneous liquid interface
class ILI_3D
{
	public: 
	
	ILI_3D();
	~ILI_3D();

	void Routine();
	
	private:

	void plot(ifstream &ifs, const Cell &cel);
	int count_geometry_number;

	double** interface;
	int** zindex;
	double*** gradient;

};


#endif
